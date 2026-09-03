#!/usr/bin/env python3
"""Restore a semantic R10 release archive to canonical project-relative paths."""

from __future__ import annotations

import argparse
import hashlib
import shutil
import tempfile
import zipfile
from pathlib import Path

import pandas as pd


MANIFEST = "release_restore_manifest.tsv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def safe_target(root: Path, relative: str) -> Path:
    target = (root / relative).resolve()
    try:
        target.relative_to(root.resolve())
    except ValueError as exc:
        raise ValueError(f"Restore target escapes project root: {relative}") from exc
    return target


def restore_from_directory(source: Path, root: Path) -> int:
    manifest_path = source / MANIFEST
    if not manifest_path.exists():
        raise FileNotFoundError(f"Missing {manifest_path}")
    manifest = pd.read_csv(manifest_path, sep="\t", dtype=str)
    restored = 0
    for row in manifest.itertuples(index=False):
        origin = (source / row.archive_path).resolve()
        if not origin.exists():
            raise FileNotFoundError(f"Release asset missing: {row.archive_path}")
        observed = sha256(origin)
        if observed.lower() != str(row.sha256).lower():
            raise ValueError(
                f"Release SHA-256 mismatch for {row.archive_path}: {observed}"
            )
        target = safe_target(root, row.project_path)
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists() and sha256(target).lower() == observed.lower():
            continue
        shutil.copy2(origin, target)
        if sha256(target).lower() != observed.lower():
            raise IOError(f"Restored SHA-256 mismatch for {row.project_path}")
        restored += 1
    print(f"Verified {len(manifest)} release assets; restored {restored} files.")
    return len(manifest)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--data-archive", type=Path)
    parser.add_argument("--source-root", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    if args.data_archive:
        archive = args.data_archive.resolve()
        if not archive.exists():
            raise FileNotFoundError(archive)
        with tempfile.TemporaryDirectory(prefix="r10_restore_") as temporary:
            source = Path(temporary)
            with zipfile.ZipFile(archive) as handle:
                names = set(handle.namelist())
                if MANIFEST not in names:
                    raise FileNotFoundError(f"{MANIFEST} is absent from {archive}")
                handle.extractall(source)
            restore_from_directory(source, root)
    else:
        source = (args.source_root or root).resolve()
        restore_from_directory(source, root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
