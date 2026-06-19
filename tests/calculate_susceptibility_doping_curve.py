from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from config import load_config, read_range, section
from susceptibility import (
    PatchData,
    calculate_susceptibility,
    compact_float_tag,
    default_patch_output_path,
    fermi_occupation,
    load_patch,
)


def default_output_path(config: dict, patch_path: Path, curve_cfg: dict) -> Path:
    data_dir = Path(str(config.get("io", {}).get("data_dir", "data")))
    q = curve_cfg.get("q_index", [0, 0])
    temp = compact_float_tag(float(curve_cfg.get("temperature_K", config.get("patch", {}).get("temperature_K", 0.0))))
    return data_dir / f"chi_doping_curve_{patch_path.stem}_q{int(q[0])}_{int(q[1])}_T{temp}.csv"


def patch_with_mu(patch: PatchData, mu_eV: float, temperature_K: float, occ: np.ndarray | None = None) -> PatchData:
    if occ is None:
        occ = fermi_occupation(patch.evals, mu_eV, temperature_K)
    meta = dict(patch.meta)
    meta["mu_eV"] = float(mu_eV)
    meta["temperature_K"] = float(temperature_K)
    return PatchData(
        meta=meta,
        iq=patch.iq,
        jq=patch.jq,
        kx=patch.kx,
        ky=patch.ky,
        evals=patch.evals,
        occ=occ,
        evecs=patch.evecs,
        dk1=patch.dk1,
        dk2=patch.dk2,
        center_k=patch.center_k,
        mu_values=patch.mu_values,
        occ_mu=patch.occ_mu,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Calculate chi(q) versus doping from one BLG patch .npz file.")
    parser.add_argument("-c", "--config", default="configs/bilayer_sk.json", help="JSON config path")
    parser.add_argument("-p", "--patch", default=None, help="input patch .npz path; overrides config")
    parser.add_argument("-o", "--output", default=None, help="output CSV path; overrides config/default")
    args = parser.parse_args()

    config = load_config(args.config)
    curve_cfg = dict(section(config, "susceptibility_doping_curve"))
    sus_base = dict(config.get("susceptibility", {}))

    patch_cfg = dict(config.get("hamiltonian", {}))
    patch_cfg.update(config.get("patch", {}))
    patch_path = Path(args.patch or curve_cfg.get("patch_file", config.get("patch", {}).get("output", default_patch_output_path(config, patch_cfg))))
    patch = load_patch(patch_path)

    if "mu_eV" in curve_cfg:
        mu_values = read_range(curve_cfg["mu_eV"])
        occ_values = [None] * len(mu_values)
    elif patch.mu_values is not None and len(patch.mu_values) > 0:
        mu_values = [float(x) for x in patch.mu_values]
        if patch.occ_mu is not None and len(patch.occ_mu) == len(mu_values):
            occ_values = [patch.occ_mu[i] for i in range(len(mu_values))]
        else:
            occ_values = [None] * len(mu_values)
    else:
        mu_values = read_range({"min": -0.05, "max": 0.05, "num": 101})
        occ_values = [None] * len(mu_values)
    temperature_K = float(curve_cfg.get("temperature_K", config.get("patch", {}).get("temperature_K", 0.0)))
    q = curve_cfg.get("q_index", [0, 0])
    q_i = int(q[0])
    q_j = int(q[1])

    sus_cfg = dict(sus_base)
    sus_cfg.update(
        {
            "q_list": [{"iq": q_i, "jq": q_j}],
            "eta_eV": float(curve_cfg.get("eta_eV", sus_base.get("eta_eV", 1e-4))),
            "use_form_factor": bool(curve_cfg.get("use_form_factor", sus_base.get("use_form_factor", True))),
            "diagonal_band_only": bool(curve_cfg.get("diagonal_band_only", sus_base.get("diagonal_band_only", False))),
            "band_indices": curve_cfg.get("band_indices", sus_base.get("band_indices", [])),
        }
    )

    rows = []
    for mu, occ in zip(mu_values, occ_values, strict=True):
        patch_mu = patch_with_mu(patch, float(mu), temperature_K, occ)
        chi_rows, _ = calculate_susceptibility(patch_mu, sus_cfg)
        row = chi_rows[0]
        rows.append([float(mu), row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7]])

    rows_arr = np.asarray(rows, dtype=float)
    out_path = Path(args.output or curve_cfg.get("output", default_output_path(config, patch_path, curve_cfg)))
    out_path.parent.mkdir(parents=True, exist_ok=True)

    meta = {
        "format": "blg_susceptibility_doping_curve_csv_v1",
        "patch_file": str(patch_path),
        "q_index": [q_i, q_j],
        "temperature_K": temperature_K,
        "susceptibility": sus_cfg,
        "columns": ["mu_eV", "iq", "jq", "qx_Ainv", "qy_Ainv", "chi_re", "chi_im", "nKpair", "nK"],
    }
    header = "\n".join(
        [
            "metadata = " + json.dumps(meta, sort_keys=True),
            "patch_metadata = " + json.dumps(patch.meta, sort_keys=True),
            "mu_eV,iq,jq,qx_Ainv,qy_Ainv,chi_re,chi_im,nKpair,nK",
        ]
    )
    np.savetxt(
        out_path,
        rows_arr,
        delimiter=",",
        header=header,
        comments="# ",
        fmt=["%.12e", "%d", "%d", "%.12e", "%.12e", "%.12e", "%.12e", "%d", "%d"],
    )

    print(f"Read: {patch_path}")
    print(f"Wrote: {out_path}")
    print(f"mu_count: {len(mu_values)}")
    print(f"q_index: [{q_i}, {q_j}]")
    print(f"chi_re span: [{rows_arr[:, 5].min():.8g}, {rows_arr[:, 5].max():.8g}]")
    print(f"chi_im span: [{rows_arr[:, 6].min():.8g}, {rows_arr[:, 6].max():.8g}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
