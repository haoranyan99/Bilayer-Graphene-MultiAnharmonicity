from __future__ import annotations

from dataclasses import dataclass
from math import exp

import numpy as np
from numpy.typing import NDArray

from structure import BilayerGrapheneStructure


FloatArray = NDArray[np.float64]
ComplexArray = NDArray[np.complex128]


@dataclass
class SKParams:
    Vpp_pi0: float = -2.81
    Vpp_sigma0: float = 0.48
    q_pi: float = 3.1451
    a_pi: float = 1.418
    q_sigma: float = 7.428
    a_sigma: float = 3.349
    r_c: float = 6.14
    l_c: float = 0.265
    r_cut: float = 8.5
    nR_max: int = 6
    include_intralayer: bool = True
    include_interlayer: bool = True
    r_min: float = 1e-6
    t_min: float = 1e-12

    @classmethod
    def from_dict(cls, data: dict[str, object] | None) -> "SKParams":
        params = cls()
        if not data:
            return params
        for key, value in data.items():
            if hasattr(params, key):
                setattr(params, key, value)
        params.nR_max = int(params.nR_max)
        return params


@dataclass(frozen=True)
class Hopping:
    m: int
    n: int
    x: int
    y: int
    z: int
    t: float


class BilayerSKModel:
    """pz-only Slater-Koster tight-binding model for AB bilayer graphene."""

    def __init__(
        self,
        structure: BilayerGrapheneStructure,
        params: SKParams | None = None,
        onsite_layer_gap_eV: float = 0.0,
    ) -> None:
        self.structure = structure
        self.params = params or SKParams()
        self.onsite_layer_gap_eV = float(onsite_layer_gap_eV)
        self._hoppings: list[Hopping] | None = None

    @property
    def hoppings(self) -> list[Hopping]:
        if self._hoppings is None:
            self._hoppings = self.generate_hoppings()
        return self._hoppings

    def generate_hoppings(self, first_layer_shear: FloatArray | None = None) -> list[Hopping]:
        base = self.structure.positions.copy()
        if first_layer_shear is not None:
            shear = np.asarray(first_layer_shear, dtype=float)
            if shear.shape != (2,):
                raise ValueError("first_layer_shear must have shape (2,)")
            for i, site in enumerate(self.structure.atoms):
                if site.layer == 0:
                    base[i, :2] += shear

        enlarged, info = self._build_enlarged_atoms(base)
        hops: list[Hopping] = []
        meta = self.structure.atoms
        p = self.params

        for m, r_m in enumerate(base):
            for e, r_e in enumerate(enlarged):
                dr = r_e - r_m
                r = float(np.linalg.norm(dr))
                if r < p.r_min or r > p.r_cut:
                    continue
                n, tx, ty, tz = info[e]
                same_layer = meta[m].layer == meta[n].layer
                if same_layer and not p.include_intralayer:
                    continue
                if (not same_layer) and not p.include_interlayer:
                    continue
                t = self._sk_integral_pz_pz(dr, r)
                if abs(t) >= p.t_min:
                    hops.append(Hopping(m, n, tx, ty, tz, t))
        return hops

    def _build_enlarged_atoms(self, base: FloatArray) -> tuple[FloatArray, list[tuple[int, int, int, int]]]:
        p = self.params
        cells: list[FloatArray] = []
        info: list[tuple[int, int, int, int]] = []
        for x in range(-p.nR_max, p.nR_max + 1):
            for y in range(-p.nR_max, p.nR_max + 1):
                shift2 = x * self.structure.a1 + y * self.structure.a2
                shift3 = np.array([shift2[0], shift2[1], 0.0], dtype=float)
                for n, r in enumerate(base):
                    cells.append(r + shift3)
                    info.append((n, x, y, 0))
        return np.vstack(cells), info

    def soft_cut(self, r: float) -> float:
        x = (r - self.params.r_c) / self.params.l_c
        if x > 50.0:
            return 0.0
        if x < -50.0:
            return 1.0
        return 1.0 / (1.0 + exp(x))

    def vpp_pi(self, r: float) -> float:
        p = self.params
        return p.Vpp_pi0 * exp(p.q_pi * (1.0 - r / p.a_pi)) * self.soft_cut(r)

    def vpp_sigma(self, r: float) -> float:
        p = self.params
        return p.Vpp_sigma0 * exp(p.q_sigma * (1.0 - r / p.a_sigma)) * self.soft_cut(r)

    def _sk_integral_pz_pz(self, dr: FloatArray, r: float) -> float:
        z2_over_r2 = float(dr[2] * dr[2] / (r * r))
        return self.vpp_pi(r) * (1.0 - z2_over_r2) + self.vpp_sigma(r) * z2_over_r2

    def build_hk(self, k: FloatArray, enforce_hermitian: bool = True) -> ComplexArray:
        norb = len(self.structure.atoms)
        h = np.zeros((norb, norb), dtype=np.complex128)
        k = np.asarray(k, dtype=float)

        for hop in self.hoppings:
            r2 = hop.x * self.structure.a1 + hop.y * self.structure.a2
            phase = np.exp(1j * float(np.dot(k, r2)))
            h[hop.m, hop.n] += hop.t * phase

        if self.onsite_layer_gap_eV != 0.0:
            h += self._onsite_layer_gap_matrix()

        if enforce_hermitian:
            h = 0.5 * (h + h.conj().T)
        return h

    def _onsite_layer_gap_matrix(self) -> ComplexArray:
        # onsite_layer_gap_eV is the top-bottom layer potential difference.
        diag = np.array([-0.5, -0.5, 0.5, 0.5], dtype=float) * self.onsite_layer_gap_eV
        return np.diag(diag.astype(np.complex128))

    def bands_at_k(self, k: FloatArray) -> FloatArray:
        return np.linalg.eigvalsh(self.build_hk(k, enforce_hermitian=True))

    def bands_along_path(self, k_list: FloatArray) -> FloatArray:
        return np.vstack([self.bands_at_k(k) for k in k_list])
