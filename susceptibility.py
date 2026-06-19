from __future__ import annotations

import json
from dataclasses import dataclass
from math import pi
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

from model_factory import choose_kpoint, read_onsite_layer_gap_eV
from structure import BilayerGrapheneStructure


FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]
ComplexArray = NDArray[np.complex128]

KB_EV_PER_K = 8.617333262145e-5


@dataclass(frozen=True)
class PatchData:
    meta: dict[str, Any]
    iq: IntArray
    jq: IntArray
    kx: FloatArray
    ky: FloatArray
    evals: FloatArray
    occ: FloatArray
    evecs: ComplexArray
    dk1: FloatArray
    dk2: FloatArray
    center_k: FloatArray
    mu_values: FloatArray | None = None
    occ_mu: FloatArray | None = None


def compact_float_tag(x: float) -> str:
    if abs(x) < 5e-13:
        x = 0.0
    return f"{x:.10g}".replace("-", "m").replace(".", "p")


def fermi_occupation(evals: FloatArray, mu_eV: float, temperature_K: float) -> FloatArray:
    if temperature_K <= 0.0:
        occ = np.zeros_like(evals, dtype=float)
        occ[evals < mu_eV] = 1.0
        occ[np.isclose(evals, mu_eV, atol=1e-12, rtol=0.0)] = 0.5
        return occ

    beta_arg = (evals - mu_eV) / (KB_EV_PER_K * temperature_K)
    beta_arg = np.clip(beta_arg, -700.0, 700.0)
    return 1.0 / (np.exp(beta_arg) + 1.0)


def _cell_area_2d(v1: FloatArray, v2: FloatArray) -> float:
    return float(abs(v1[0] * v2[1] - v1[1] * v2[0]))


def read_scan_values(spec: dict[str, float | int] | list[float]) -> list[float]:
    if isinstance(spec, list):
        if not spec:
            raise ValueError("scan list must not be empty")
        return [float(x) for x in spec]
    xmin = float(spec["min"])
    xmax = float(spec["max"])
    num = int(spec["num"])
    if num < 1:
        raise ValueError("scan num must be >= 1")
    if num == 1:
        return [xmin]
    return [xmin + (xmax - xmin) * i / (num - 1) for i in range(num)]


def read_patch_mu_values(patch_cfg: dict[str, Any]) -> FloatArray:
    spec = patch_cfg.get("mu_eV_list", patch_cfg.get("mu_list_eV", [float(patch_cfg.get("mu_eV", 0.0))]))
    return np.asarray(read_scan_values(spec), dtype=float)


def energy_zero_shift_eV(patch_cfg: dict[str, Any], iq: IntArray, jq: IntArray, evals: FloatArray) -> tuple[float, dict[str, Any]]:
    if "energy_shift_to_midgap" in patch_cfg:
        do_shift = bool(patch_cfg["energy_shift_to_midgap"])
    else:
        do_shift = str(dict(patch_cfg.get("energy_zero", {"mode": "middle_at_center"})).get("mode", "middle_at_center")).lower() not in {"none", "raw", "off"}

    if not do_shift:
        return 0.0, {"energy_shift_to_midgap": False, "shift_eV": 0.0}

    nbands = evals.shape[1]
    if nbands < 2:
        raise ValueError("energy_shift_to_midgap requires at least two bands")

    center_matches = np.flatnonzero((iq == 0) & (jq == 0))
    if len(center_matches) != 1:
        raise ValueError("energy_shift_to_midgap requires exactly one patch center point with iq=0,jq=0")
    idx = int(center_matches[0])
    b0 = nbands // 2 - 1
    b1 = nbands // 2
    shift = 0.5 * float(evals[idx, b0] + evals[idx, b1])
    return shift, {
        "energy_shift_to_midgap": True,
        "shift_eV": shift,
        "reference_iq": 0,
        "reference_jq": 0,
        "bands": [b0, b1],
    }


def build_patch_grid(st: BilayerGrapheneStructure, patch_cfg: dict[str, Any]) -> tuple[IntArray, IntArray, FloatArray, FloatArray, FloatArray]:
    mesh_cfg = dict(patch_cfg.get("mesh", {}))
    nx = int(patch_cfg.get("Nx", mesh_cfg.get("Nx", 41)))
    ny = int(patch_cfg.get("Ny", mesh_cfg.get("Ny", nx)))
    if nx < 1 or ny < 1:
        raise ValueError("patch Nx and Ny must be >= 1")
    if nx % 2 != 1 or ny % 2 != 1:
        raise ValueError("patch Nx and Ny must be odd so that iq=0,jq=0 is included")

    center_name = str(patch_cfg.get("center", mesh_cfg.get("center", "K")))
    center_k = choose_kpoint(st, center_name)

    dx = float(patch_cfg.get("dx_Ainv", mesh_cfg.get("dx_Ainv", mesh_cfg.get("dk_Ainv", 0.002))))
    dy = float(patch_cfg.get("dy_Ainv", mesh_cfg.get("dy_Ainv", mesh_cfg.get("dk_Ainv", dx))))
    if dx <= 0.0 or dy <= 0.0:
        raise ValueError("patch dx_Ainv and dy_Ainv must be positive")
    dk1 = np.array([dx, 0.0], dtype=float)
    dk2 = np.array([0.0, dy], dtype=float)

    ix = np.arange(-(nx // 2), nx // 2 + 1, dtype=np.int64)
    iy = np.arange(-(ny // 2), ny // 2 + 1, dtype=np.int64)
    iq_grid, jq_grid = np.meshgrid(ix, iy, indexing="ij")
    iq = iq_grid.ravel()
    jq = jq_grid.ravel()
    kxy = center_k[None, :] + iq[:, None] * dk1[None, :] + jq[:, None] * dk2[None, :]
    return iq, jq, kxy[:, 0], kxy[:, 1], center_k


def default_patch_output_path(config: dict[str, Any], patch_cfg: dict[str, Any]) -> Path:
    data_dir = Path(str(config.get("io", {}).get("data_dir", "data")))
    mesh_cfg = dict(patch_cfg.get("mesh", {}))
    model = str(patch_cfg.get("model", config.get("hamiltonian", {}).get("model", "tb"))).lower()
    gap = compact_float_tag(read_onsite_layer_gap_eV(patch_cfg))
    temp = compact_float_tag(float(patch_cfg.get("temperature_K", 0.0)))
    center = str(patch_cfg.get("center", mesh_cfg.get("center", "K")))
    nx = int(patch_cfg.get("Nx", mesh_cfg.get("Nx", 41)))
    ny = int(patch_cfg.get("Ny", mesh_cfg.get("Ny", nx)))
    dx_val = float(patch_cfg.get("dx_Ainv", mesh_cfg.get("dx_Ainv", mesh_cfg.get("dk_Ainv", 0.002))))
    dy_val = float(patch_cfg.get("dy_Ainv", mesh_cfg.get("dy_Ainv", mesh_cfg.get("dk_Ainv", dx_val))))
    dx = compact_float_tag(dx_val)
    dy = compact_float_tag(dy_val)
    return data_dir / f"patch_{model}_gap{gap}_{center}_cart_Nx{nx}_Ny{ny}_dx{dx}_dy{dy}_T{temp}.npz"


def build_patch_data(
    config: dict[str, Any],
    patch_cfg: dict[str, Any],
    st: BilayerGrapheneStructure,
    model: Any,
) -> PatchData:
    iq, jq, kx, ky, center_k = build_patch_grid(st, patch_cfg)
    nk = len(kx)
    h0 = model.build_hk(np.array([kx[0], ky[0]], dtype=float))
    nbasis = h0.shape[0]
    evals = np.empty((nk, nbasis), dtype=float)
    evecs = np.empty((nk, nbasis, nbasis), dtype=np.complex128)

    for i, (kx_i, ky_i) in enumerate(zip(kx, ky, strict=True)):
        h = model.build_hk(np.array([kx_i, ky_i], dtype=float))
        eval_i, evec_i = np.linalg.eigh(h)
        evals[i, :] = eval_i
        evecs[i, :, :] = evec_i

    energy_shift, energy_zero_meta = energy_zero_shift_eV(patch_cfg, iq, jq, evals)
    evals = evals - energy_shift

    mu_values = read_patch_mu_values(patch_cfg)
    reference_mu_eV = float(patch_cfg.get("reference_mu_eV", patch_cfg.get("mu_eV", 0.0)))
    temperature_K = float(patch_cfg.get("temperature_K", 0.0))
    occ = fermi_occupation(evals, reference_mu_eV, temperature_K)
    occ_mu = np.stack([fermi_occupation(evals, mu, temperature_K) for mu in mu_values], axis=0)

    mesh_cfg = dict(patch_cfg.get("mesh", {}))
    nx = int(patch_cfg.get("Nx", mesh_cfg.get("Nx", 41)))
    ny = int(patch_cfg.get("Ny", mesh_cfg.get("Ny", nx)))
    dx = float(patch_cfg.get("dx_Ainv", mesh_cfg.get("dx_Ainv", mesh_cfg.get("dk_Ainv", 0.002))))
    dy = float(patch_cfg.get("dy_Ainv", mesh_cfg.get("dy_Ainv", mesh_cfg.get("dk_Ainv", dx))))
    dk1 = np.array([dx, 0.0], dtype=float)
    dk2 = np.array([0.0, dy], dtype=float)

    kcell_area = _cell_area_2d(dk1, dk2)
    integration_weight = kcell_area / (2.0 * pi) ** 2
    meta = {
        "format": "blg_patch_npz_v1",
        "model": str(patch_cfg.get("model", config.get("hamiltonian", {}).get("model", "tb"))).lower(),
        "onsite_layer_gap_eV": read_onsite_layer_gap_eV(patch_cfg),
        "valley": int(patch_cfg.get("valley", config.get("hamiltonian", {}).get("valley", +1))),
        "reference_mu_eV": reference_mu_eV,
        "mu_values_eV": mu_values.tolist(),
        "temperature_K": temperature_K,
        "mesh": {
            "type": "cartesian",
            "center": str(patch_cfg.get("center", mesh_cfg.get("center", "K"))),
            "Nx": nx,
            "Ny": ny,
            "dx_Ainv": dx,
            "dy_Ainv": dy,
        },
        "nk_total": int(nk),
        "nbands": int(nbasis),
        "basis": "A1,B1,A2,B2",
        "energy_zero": energy_zero_meta,
        "center_k_Ainv": center_k.tolist(),
        "dk1_Ainv": dk1.tolist(),
        "dk2_Ainv": dk2.tolist(),
        "kcell_area_Ainv2": kcell_area,
        "integration_weight": integration_weight,
        "integration_weight_note": "Default chi sum uses d^2k/(2*pi)^2 for each k point.",
    }
    return PatchData(meta, iq, jq, kx, ky, evals, occ, evecs, dk1, dk2, center_k, mu_values, occ_mu)


def save_patch(path: str | Path, patch: PatchData) -> Path:
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        out_path,
        meta_json=np.array(json.dumps(patch.meta, indent=2), dtype=np.str_),
        iq=patch.iq,
        jq=patch.jq,
        kx=patch.kx,
        ky=patch.ky,
        evals=patch.evals,
        occ=patch.occ,
        evecs=patch.evecs,
        dk1=patch.dk1,
        dk2=patch.dk2,
        center_k=patch.center_k,
        mu_values=patch.mu_values if patch.mu_values is not None else np.array([], dtype=float),
        occ_mu=patch.occ_mu if patch.occ_mu is not None else np.array([], dtype=float),
    )
    return out_path


def load_patch(path: str | Path) -> PatchData:
    with np.load(path, allow_pickle=False) as data:
        meta = json.loads(str(data["meta_json"]))
        return PatchData(
            meta=meta,
            iq=data["iq"].astype(np.int64),
            jq=data["jq"].astype(np.int64),
            kx=data["kx"].astype(float),
            ky=data["ky"].astype(float),
            evals=data["evals"].astype(float),
            occ=data["occ"].astype(float),
            evecs=data["evecs"].astype(np.complex128),
            dk1=data["dk1"].astype(float),
            dk2=data["dk2"].astype(float),
            center_k=data["center_k"].astype(float),
            mu_values=data["mu_values"].astype(float) if "mu_values" in data else None,
            occ_mu=data["occ_mu"].astype(float) if "occ_mu" in data else None,
        )


def build_q_shifts(sus_cfg: dict[str, Any]) -> list[tuple[int, int]]:
    if "q_list" in sus_cfg:
        return [(int(q["iq"]), int(q["jq"])) for q in sus_cfg["q_list"]]

    q_mesh = dict(sus_cfg.get("q_mesh", {"Nq_half": 4}))
    nq_half = int(q_mesh.get("Nq_half", q_mesh.get("Nq", 4)))
    if nq_half < 0:
        raise ValueError("susceptibility.q_mesh.Nq_half must be >= 0")
    vals = range(-nq_half, nq_half + 1)
    return [(iq, jq) for iq in vals for jq in vals]


def calculate_susceptibility(patch: PatchData, sus_cfg: dict[str, Any]) -> tuple[FloatArray, dict[str, Any]]:
    eta_eV = float(sus_cfg.get("eta_eV", 1e-4))
    if eta_eV <= 0.0:
        raise ValueError("susceptibility.eta_eV must be positive")

    use_form_factor = bool(sus_cfg.get("use_form_factor", True))
    diagonal_band_only = bool(sus_cfg.get("diagonal_band_only", False))
    band_indices = sus_cfg.get("band_indices", None)
    if band_indices is None or len(band_indices) == 0:
        bands = list(range(patch.evals.shape[1]))
    else:
        bands = [int(i) for i in band_indices]

    weight = float(sus_cfg.get("integration_weight", patch.meta.get("integration_weight", 1.0)))
    index_of = {(int(iq), int(jq)): idx for idx, (iq, jq) in enumerate(zip(patch.iq, patch.jq, strict=True))}
    q_shifts = build_q_shifts(sus_cfg)

    rows: list[list[float]] = []
    for q_i, q_j in q_shifts:
        qvec = q_i * patch.dk1 + q_j * patch.dk2
        chi = 0.0 + 0.0j
        n_pair = 0
        for idx1, (i1, j1) in enumerate(zip(patch.iq, patch.jq, strict=True)):
            idx2 = index_of.get((int(i1) + q_i, int(j1) + q_j))
            if idx2 is None:
                continue
            n_pair += 1
            for b in bands:
                eb = patch.evals[idx1, b]
                fb = patch.occ[idx1, b]
                ub = patch.evecs[idx1, :, b]
                for m in bands:
                    if diagonal_band_only and m != b:
                        continue
                    em = patch.evals[idx2, m]
                    fm = patch.occ[idx2, m]
                    form_factor = 1.0
                    if use_form_factor:
                        um = patch.evecs[idx2, :, m]
                        form_factor = float(abs(np.vdot(ub, um)) ** 2)
                    chi += weight * form_factor * (fm - fb) / complex(eb - em, eta_eV)

        rows.append([q_i, q_j, float(qvec[0]), float(qvec[1]), float(chi.real), float(chi.imag), n_pair, len(patch.iq)])

    meta = {
        "format": "blg_susceptibility_csv_v1",
        "patch_format": patch.meta.get("format", "unknown"),
        "eta_eV": eta_eV,
        "use_form_factor": use_form_factor,
        "diagonal_band_only": diagonal_band_only,
        "band_indices": bands,
        "integration_weight": weight,
        "columns": ["iq", "jq", "qx_Ainv", "qy_Ainv", "chi_re", "chi_im", "nKpair", "nK"],
    }
    return np.asarray(rows, dtype=float), meta


def default_chi_output_path(config: dict[str, Any], patch_path: str | Path, sus_cfg: dict[str, Any]) -> Path:
    data_dir = Path(str(config.get("io", {}).get("data_dir", "data")))
    eta = compact_float_tag(float(sus_cfg.get("eta_eV", 1e-4)))
    q_shifts = build_q_shifts(sus_cfg)
    if "q_list" in sus_cfg and len(q_shifts) <= 3:
        q_tag = "q" + "_".join(f"{iq}_{jq}" for iq, jq in q_shifts)
    else:
        q_mesh = dict(sus_cfg.get("q_mesh", {"Nq_half": 4}))
        q_tag = f"qmesh_Nh{int(q_mesh.get('Nq_half', q_mesh.get('Nq', 4)))}"
    patch_stem = Path(patch_path).stem
    return data_dir / f"chi_{patch_stem}_eta{eta}_{q_tag}.csv"


def save_chi_csv(path: str | Path, rows: FloatArray, meta: dict[str, Any], patch_meta: dict[str, Any]) -> Path:
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    header = "\n".join(
        [
            "metadata = " + json.dumps(meta, sort_keys=True),
            "patch_metadata = " + json.dumps(patch_meta, sort_keys=True),
            "iq,jq,qx_Ainv,qy_Ainv,chi_re,chi_im,nKpair,nK",
        ]
    )
    np.savetxt(
        out_path,
        rows,
        delimiter=",",
        header=header,
        comments="# ",
        fmt=["%d", "%d", "%.12e", "%.12e", "%.12e", "%.12e", "%d", "%d"],
    )
    return out_path
