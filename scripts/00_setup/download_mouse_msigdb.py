"""Download and verify the frozen Mouse MSigDB pathway collections."""

from __future__ import annotations

import csv
import hashlib
import urllib.request
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
VERSION = "2026.1.Mm"
BASE_URL = f"https://data.broadinstitute.org/gsea-msigdb/msigdb/release/{VERSION}"
OUTPUT_DIR = PROJECT_ROOT / "external_data" / "msigdb" / VERSION
MANIFEST_CANDIDATES = (
    PROJECT_ROOT / "metadata" / "mouse_msigdb_sha256.tsv",
    PROJECT_ROOT / "01_manifest" / "mouse_msigdb_sha256.tsv",
)

FILES = (
    f"mh.all.v{VERSION}.symbols.gmt",
    f"m2.cp.reactome.v{VERSION}.symbols.gmt",
    f"m5.go.bp.v{VERSION}.symbols.gmt",
)


def load_expected_manifest() -> dict[str, dict[str, str]]:
    manifest = next((path for path in MANIFEST_CANDIDATES if path.exists()), None)
    if manifest is None:
        raise FileNotFoundError(
            "Frozen Mouse MSigDB checksum manifest is missing; expected one of "
            + ", ".join(str(path) for path in MANIFEST_CANDIDATES)
        )
    with manifest.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    expected = {row["filename"]: row for row in rows}
    missing = sorted(set(FILES) - set(expected))
    if missing:
        raise RuntimeError(
            f"Checksum manifest {manifest} is missing entries for: {missing}"
        )
    return expected


def main() -> None:
    expected = load_expected_manifest()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for filename in FILES:
        url = f"{BASE_URL}/{filename}"
        destination = OUTPUT_DIR / filename
        if not destination.exists() or destination.stat().st_size == 0:
            request = urllib.request.Request(
                url,
                headers={"User-Agent": "BPD-PH-reanalysis/1.0"},
            )
            with urllib.request.urlopen(request, timeout=300) as response:
                destination.write_bytes(response.read())
        payload = destination.read_bytes()
        if b"\t" not in payload:
            raise RuntimeError(f"Downloaded file does not look like GMT: {destination}")
        observed_hash = hashlib.sha256(payload).hexdigest()
        expected_hash = expected[filename]["sha256"].lower()
        expected_bytes = int(expected[filename]["bytes"])
        if len(payload) != expected_bytes or observed_hash != expected_hash:
            raise RuntimeError(
                f"Checksum mismatch for {filename}: expected "
                f"{expected_bytes} bytes/{expected_hash}, observed "
                f"{len(payload)} bytes/{observed_hash}"
            )
        print(
            f"{filename}: {len(payload):,} bytes; "
            f"SHA-256 verified ({observed_hash})"
        )


if __name__ == "__main__":
    main()
