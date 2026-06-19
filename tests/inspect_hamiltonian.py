from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np

from config import load_config, section
from model_factory import build_model, build_structure, choose_kpoint


def main() -> int:
    parser = argparse.ArgumentParser(description="Inspect the AB bilayer graphene SK Hamiltonian.")
    parser.add_argument("-c", "--config", default="configs/bilayer_sk.json", help="JSON config path")
    parser.add_argument("--k-point", default=None, help="Gamma, M, K, or Kp; overrides config")
    parser.add_argument("--precision", type=int, default=6, help="printed numeric precision")
    args = parser.parse_args()

    config = load_config(args.config)
    h_cfg = section(config, "hamiltonian")
    st = build_structure(config)
    model = build_model(config, h_cfg, st)

    k_name = args.k_point or str(h_cfg.get("k_point", "K"))
    k = choose_kpoint(st, k_name)
    hk = model.build_hk(k)
    evals = np.linalg.eigvalsh(hk)

    np.set_printoptions(precision=args.precision, suppress=False, linewidth=160)
    print("=== inspect_hamiltonian ===")
    print(f"config    = {args.config}")
    print(f"model     = {h_cfg.get('model', 'sk')}")
    print(f"k_point   = {k_name}")
    print(f"k_Ainv    = {k}")
    print(f"onsite_layer_gap_eV = {model.onsite_layer_gap_eV:g}")
    print(f"atoms     = {len(st.atoms)}")
    for i, site in enumerate(st.atoms):
        print(f"  {i}: layer={site.layer} sub={site.sublattice} r_A={site.r}")
    if hasattr(model, "hoppings"):
        print(f"hoppings  = {len(model.hoppings)}")
    elif hasattr(model, "params"):
        p = model.params
        print(
            "tb_params = "
            f"gamma0={p.gamma0:g}, gamma1={p.gamma1:g}, "
            f"gamma3={p.gamma3:g}, gamma4={p.gamma4:g}, "
            f"delta_prime={p.delta_prime:g}"
        )
        print(f"valley    = {model.valley:+d}")
    print("H(k) eV =")
    print(hk)
    print("eigvalsh(H) eV =")
    print(evals)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
