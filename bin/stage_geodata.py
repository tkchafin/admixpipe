#!/usr/bin/env python3
"""
stage_geodata.py  ─ stage vector layers, rewrite JSON paths

Example
-------
python stage_geodata.py \
    --json      assets/streams_overlay.json \
    --outdir    geo_data_files \
    --out       geo_data_files/geo_data_layers.json \
    --base_path $PWD        # optional
"""

from __future__ import annotations
import argparse, glob, json, shutil
from pathlib import Path


# ───────────────────────── helpers ─────────────────────────
def copy_layer(src: Path, dest_dir: Path) -> Path:
    """
    Copy *src* and any sibling files that share the same stem (foo.*)
    into *dest_dir*.  Returns the Path to the copied primary file.
    """
    dest_dir.mkdir(parents=True, exist_ok=True)

    stem = src.with_suffix("")
    siblings = glob.glob(str(stem) + ".*")
    files = siblings or [str(src)]      # fall back to single file

    for f in files:
        shutil.copy2(f, dest_dir)

    return dest_dir / src.name


def stage_layers(src_json: Path, out_dir: Path, base: Path | None) -> list[dict]:
    layers = json.loads(src_json.read_text())
    if isinstance(layers, dict):
        layers = [layers]

    default_base = src_json.parent.resolve()

    for layer in layers:
        raw_path = Path(layer["path"])
        src      = (
            raw_path
            if raw_path.is_absolute()
            else (base or default_base) / raw_path
        )
        src = src.expanduser().resolve()

        if not src.exists():
            raise FileNotFoundError(src)

        staged = copy_layer(src, out_dir)
        layer["path"] = str(staged)     # rewrite path

    return layers


# ───────────────────────── main ────────────────────────────
def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--json", required=True,  type=Path, help="Input JSON with layer specs")
    p.add_argument("--outdir", required=True, type=Path, help="Directory where files are copied")
    p.add_argument("--out", required=True,    type=Path, help="Path for rewritten JSON")
    p.add_argument("--base_path", type=Path, help="Base path for resolving relative layer paths")
    args = p.parse_args(argv)

    args.outdir.mkdir(parents=True, exist_ok=True)

    updated = stage_layers(args.json, args.outdir, args.base_path)
    args.out.write_text(json.dumps(updated, indent=2))

    print(f"✅  Staged {len(updated)} layer(s) into '{args.outdir}'")
    print(f"✅  Rewritten JSON → {args.out}'")


if __name__ == "__main__":
    main()
