"""Download the three GEO datasets and fixed external resources."""
from __future__ import annotations

import argparse
import csv
import hashlib
import tarfile
import urllib.request
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    root = args.project_root.resolve()
    manifest_path = root / "config" / "download_manifest.tsv"
    with manifest_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        destination = root / row["relative_path"]
        destination.parent.mkdir(parents=True, exist_ok=True)
        if args.overwrite or not destination.exists() or destination.stat().st_size == 0:
            request = urllib.request.Request(
                row["source_url"],
                headers={"User-Agent": "BPD-PH-reanalysis/1.0"},
            )
            with urllib.request.urlopen(request, timeout=600) as response:
                destination.write_bytes(response.read())
        digest = hashlib.sha256(destination.read_bytes()).hexdigest()
        expected = row["sha256"].strip().lower()
        if not expected:
            raise RuntimeError(
                f"Missing frozen SHA-256 for {destination.relative_to(root)}"
            )
        if digest != expected:
            raise RuntimeError(
                f"SHA-256 mismatch for {destination.relative_to(root)}: "
                f"expected {expected}, observed {digest}"
            )
        print(
            f"{row['accession']}: {destination.relative_to(root)} "
            f"(SHA-256 verified: {digest})"
        )
        if row["extract_tar"].lower() == "true":
            extraction_dir = root / row["extract_to"]
            extraction_dir.mkdir(parents=True, exist_ok=True)
            with tarfile.open(destination, "r") as archive:
                archive.extractall(extraction_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
