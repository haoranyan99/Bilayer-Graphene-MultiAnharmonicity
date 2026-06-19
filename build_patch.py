from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from config import load_config, section
from model_factory import build_model, build_structure
from susceptibility import build_patch_data, default_patch_output_path, save_patch


def main() -> int:
    parser = argparse.ArgumentParser(description="Build one BLG k-patch file with energies, occupations, and eigenvectors.")
    parser.add_argument("-c", "--config", default="configs/bilayer_sk.json", help="JSON config path")
    parser.add_argument("-o", "--output", default=None, help="output .npz path; overrides config/default")
    args = parser.parse_args()

    config = load_config(args.config)
    patch_cfg = dict(section(config, "patch"))
    model_cfg = dict(config.get("hamiltonian", {}))
    model_cfg.update(patch_cfg)

    st = build_structure(config)
    model = build_model(config, model_cfg, st)
    patch = build_patch_data(config, model_cfg, st, model)

    out_path = Path(args.output or patch_cfg.get("output", default_patch_output_path(config, model_cfg)))
    save_patch(out_path, patch)

    print(f"Wrote: {out_path}")
    print(f"model: {patch.meta['model']}")
    print(f"nk_total: {patch.meta['nk_total']}")
    print(f"nbands: {patch.meta['nbands']}")
    print(f"kcell_area_Ainv2: {patch.meta['kcell_area_Ainv2']:.12g}")
    print(f"integration_weight: {patch.meta['integration_weight']:.12g}")
    print(f"energy_zero: {patch.meta['energy_zero']}")
    if patch.mu_values is not None:
        print(f"mu_values: {len(patch.mu_values)} points from {patch.mu_values[0]:.8g} to {patch.mu_values[-1]:.8g} eV")
    if patch.occ_mu is not None:
        print(f"occ_mu shape: {patch.occ_mu.shape}")
    print(f"energy span: [{patch.evals.min():.8g}, {patch.evals.max():.8g}] eV")
    print(f"occupation sum: {patch.occ.sum():.8g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
