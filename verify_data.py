"""Verify your locally-built REAL data extracts against data_manifest.json.

NOTE: nothing in data/ ships in the repo (it is gitignored). These extracts only
exist after you build them locally (see data/README.md); this script checks that
your local build matches the manifest's recorded sha256/rows.

Usage:
    python verify_data.py            # check local extracts' sha256 vs manifest
    python verify_data.py --update   # recompute + write hashes/rows for local extracts

Only 'derived_extract' entries (the small locally-built CSVs) are hashed. External
caches (_tmp_fetch/*) and live APIs are documented in the manifest but not hashed
here — they are reproduced via their build scripts / endpoints.
"""
from __future__ import annotations
import hashlib
import json
import sys
from pathlib import Path

PKG = Path(__file__).resolve().parent
MANIFEST = PKG / "data_manifest.json"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def rows(path: Path) -> int:
    with open(path, "r", encoding="utf-8") as fh:
        return sum(1 for _ in fh) - 1  # minus header


def main() -> int:
    update = "--update" in sys.argv
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    datasets = manifest["datasets"]
    problems = 0

    for name, d in datasets.items():
        if d.get("kind") != "derived_extract":
            continue
        fp = PKG / d["local_file"]
        if not fp.exists():
            print(f"MISSING  {name}: {d['local_file']} (build it locally — see data/README.md)")
            problems += 1
            continue
        digest, n = sha256(fp), rows(fp)
        if update:
            d["sha256"], d["rows"] = digest, n
            print(f"UPDATED  {name}: sha256={digest[:12]}… rows={n}")
        else:
            ok_hash = digest == d.get("sha256")
            ok_rows = n == d.get("rows")
            if ok_hash and ok_rows:
                print(f"OK       {name}: {d['local_file']} (rows={n})")
            else:
                print(f"DRIFT    {name}: hash {'ok' if ok_hash else 'CHANGED'}, "
                      f"rows {n} vs manifest {d.get('rows')}")
                problems += 1

    if update:
        MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        print("manifest updated.")
        return 0

    print(f"\n{'all local extracts verified.' if problems == 0 else str(problems) + ' problem(s).'}")
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
