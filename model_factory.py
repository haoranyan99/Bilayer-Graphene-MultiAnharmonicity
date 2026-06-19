from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from sk_model import BilayerSKModel, SKParams
from structure import BilayerGrapheneStructure
from tb_model import BilayerTBModel, TBParams


FloatArray = NDArray[np.float64]


def build_structure(config: dict) -> BilayerGrapheneStructure:
    lattice = config.get("lattice", {})
    return BilayerGrapheneStructure(
        a0_A=float(lattice.get("a0_A", 2.46)),
        d0_A=float(lattice.get("d0_A", 3.35)),
        pressure_GPa=float(lattice.get("pressure_GPa", lattice.get("pressure", 0.0))),
        vacuum_A=float(lattice.get("vacuum_A", 10.0)),
    )


def read_onsite_layer_gap_eV(task_cfg: dict) -> float:
    if "onsite_layer_gap_eV" in task_cfg:
        return float(task_cfg["onsite_layer_gap_eV"])
    if "Dfield_eV" in task_cfg:
        return float(task_cfg["Dfield_eV"])
    return 0.0


def choose_kpoint(st: BilayerGrapheneStructure, name: str) -> FloatArray:
    key = name.strip().lower()
    if key in {"gamma", "g"}:
        return st.Gamma
    if key == "m":
        return st.M
    if key == "k":
        return st.K
    if key in {"kp", "k'", "kprime"}:
        return st.Kp
    raise ValueError(f"unknown k point {name!r}; use Gamma, M, K, or Kp")


def build_model(config: dict, task_cfg: dict, st: BilayerGrapheneStructure):
    model_name = str(task_cfg.get("model", "sk")).lower()
    onsite_gap = read_onsite_layer_gap_eV(task_cfg)
    if model_name == "sk":
        return BilayerSKModel(st, SKParams.from_dict(config.get("sk_params", {})), onsite_layer_gap_eV=onsite_gap)
    if model_name == "tb":
        return BilayerTBModel(
            st,
            TBParams.from_dict(config.get("tb_params", {})),
            onsite_layer_gap_eV=onsite_gap,
            valley=int(task_cfg.get("valley", +1)),
        )
    raise ValueError("model must be 'sk' or 'tb'")
