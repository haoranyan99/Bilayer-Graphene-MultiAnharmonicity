from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from structure import BilayerGrapheneStructure


FloatArray = NDArray[np.float64]
ComplexArray = NDArray[np.complex128]


@dataclass
class TBParams:
    # McCann-Koshino four-band bilayer graphene parameters in eV.
    gamma0: float = 3.16
    gamma1: float = 0.381
    gamma3: float = 0.38
    gamma4: float = 0.14
    delta_prime: float = 0.022

    @classmethod
    def from_dict(cls, data: dict[str, object] | None) -> "TBParams":
        params = cls()
        if not data:
            return params
        for key, value in data.items():
            if hasattr(params, key):
                setattr(params, key, float(value))
        return params


class BilayerTBModel:
    """Four-band continuum tight-binding Hamiltonian near K/K' for AB bilayer graphene.

    Basis: (A1, B1, A2, B2), matching the structure.py atom order.
    The dimer pair is B1-A2. Momentum k is absolute in A^-1; the model expands
    around K for valley=+1 and K' for valley=-1.
    """

    def __init__(
        self,
        structure: BilayerGrapheneStructure,
        params: TBParams | None = None,
        onsite_layer_gap_eV: float = 0.0,
        valley: int = +1,
    ) -> None:
        self.structure = structure
        self.params = params or TBParams()
        self.onsite_layer_gap_eV = float(onsite_layer_gap_eV)
        self.valley = +1 if valley >= 0 else -1

    def _hbar_v(self, gamma: float) -> float:
        # hbar v_i = sqrt(3) a gamma_i / 2, with a the lattice constant.
        return 0.5 * np.sqrt(3.0) * self.structure.a_A * gamma

    def _q_from_k(self, k: FloatArray) -> FloatArray:
        center = self.structure.K if self.valley > 0 else self.structure.Kp
        return np.asarray(k, dtype=float) - center

    def build_hk(self, k: FloatArray, enforce_hermitian: bool = True) -> ComplexArray:
        p = self.params
        q = self._q_from_k(k)

        xi = float(self.valley)
        pi_op = xi * q[0] + 1j * q[1]
        pi_dag = np.conjugate(pi_op)

        hv0 = self._hbar_v(p.gamma0)
        hv3 = self._hbar_v(p.gamma3)
        hv4 = self._hbar_v(p.gamma4)
        u1 = -0.5 * self.onsite_layer_gap_eV
        u2 = +0.5 * self.onsite_layer_gap_eV

        h = np.array(
            [
                [u1, hv0 * pi_dag, -hv4 * pi_dag, hv3 * pi_op],
                [hv0 * pi_op, u1 + p.delta_prime, p.gamma1, -hv4 * pi_dag],
                [-hv4 * pi_op, p.gamma1, u2 + p.delta_prime, hv0 * pi_dag],
                [hv3 * pi_dag, -hv4 * pi_op, hv0 * pi_op, u2],
            ],
            dtype=np.complex128,
        )

        if enforce_hermitian:
            h = 0.5 * (h + h.conj().T)
        return h

    def bands_at_k(self, k: FloatArray) -> FloatArray:
        return np.linalg.eigvalsh(self.build_hk(k, enforce_hermitian=True))

    def bands_along_path(self, k_list: FloatArray) -> FloatArray:
        return np.vstack([self.bands_at_k(k) for k in k_list])
