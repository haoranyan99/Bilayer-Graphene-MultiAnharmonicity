from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from config import load_config, section
from susceptibility import (
    calculate_susceptibility,
    default_chi_output_path,
    default_patch_output_path,
    load_patch,
    save_chi_csv,
)


def main() -> int:
    parser = argparse.ArgumentParser(description="Calculate BLG electron susceptibility from one patch .npz file.")
    parser.add_argument("-c", "--config", default="configs/bilayer_sk.json", help="JSON config path")
    parser.add_argument("-p", "--patch", default=None, help="input patch .npz path; overrides config")
    parser.add_argument("-o", "--output", default=None, help="output CSV path; overrides config/default")
    args = parser.parse_args()

    config = load_config(args.config)
    sus_cfg = dict(section(config, "susceptibility"))
    patch_cfg = dict(config.get("hamiltonian", {}))
    patch_cfg.update(config.get("patch", {}))

    patch_path = Path(args.patch or sus_cfg.get("patch_file", config.get("patch", {}).get("output", default_patch_output_path(config, patch_cfg))))
    patch = load_patch(patch_path)
    rows, meta = calculate_susceptibility(patch, sus_cfg)

    out_path = Path(args.output or sus_cfg.get("output", default_chi_output_path(config, patch_path, sus_cfg)))
    save_chi_csv(out_path, rows, meta, patch.meta)

    print(f"Read: {patch_path}")
    print(f"Wrote: {out_path}")
    print(f"q_count: {len(rows)}")
    print(f"chi_re span: [{rows[:, 4].min():.8g}, {rows[:, 4].max():.8g}]")
    print(f"chi_im span: [{rows[:, 5].min():.8g}, {rows[:, 5].max():.8g}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
