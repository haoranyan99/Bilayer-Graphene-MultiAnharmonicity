from __future__ import annotations

from dataclasses import dataclass
from math import pi, sqrt

import numpy as np
from numpy.typing import NDArray


FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class AtomSite:
    layer: int
    sublattice: str
    r: FloatArray


@dataclass
class KPath:
    k_list: FloatArray
    k_dist: FloatArray
    tick_pos: list[float]
    tick_labels: list[str]


class BilayerGrapheneStructure:
    """AB-stacked bilayer graphene geometry using the rhombohedral code convention."""

    def __init__(
        self,
        a0_A: float = 2.46,
        d0_A: float = 3.35,
        pressure_GPa: float = 0.0,
        vacuum_A: float = 10.0,
    ) -> None:
        if a0_A <= 0.0:
            raise ValueError("a0_A must be positive")
        if d0_A <= 0.0:
            raise ValueError("d0_A must be positive")
        if vacuum_A < 0.0:
            raise ValueError("vacuum_A must be non-negative")

        self.n_layer = 2
        self.a0_A = float(a0_A)
        self.d0_A = float(d0_A)
        self.pressure_GPa = float(pressure_GPa)
        self.vacuum_A = float(vacuum_A)

        # Same linear pressure correction as the rhombohedral C++ project.
        self.a_A = self.a0_A - 0.00197 * self.pressure_GPa
        self.d_A = self.d0_A - 0.05710 * self.pressure_GPa
        if self.a_A <= 0.0:
            raise ValueError("pressure makes in-plane lattice constant non-positive")
        if self.d_A <= 0.0:
            raise ValueError("pressure makes interlayer spacing non-positive")

        self.a1 = np.array([0.0, -self.a_A], dtype=float)
        self.a2 = np.array([0.5 * sqrt(3.0) * self.a_A, 0.5 * self.a_A], dtype=float)
        self.b1, self.b2 = self._reciprocal_vectors(self.a1, self.a2)

        self.Gamma = np.array([0.0, 0.0], dtype=float)
        self.M = 0.5 * self.b1
        self.K = (self.b1 + self.b2) / 3.0
        self.Kp = 2.0 * (self.b1 + self.b2) / 3.0

        self._frac_A = np.array([0.0, 0.0], dtype=float)
        self._frac_B = np.array([1.0 / 3.0, 2.0 / 3.0], dtype=float)
        self._frac_C = np.array([2.0 / 3.0, 1.0 / 3.0], dtype=float)
        self.atoms = self._build_ab_atoms()

    @staticmethod
    def _reciprocal_vectors(a1: FloatArray, a2: FloatArray) -> tuple[FloatArray, FloatArray]:
        area = a1[0] * a2[1] - a1[1] * a2[0]
        if abs(area) < 1e-15:
            raise ValueError("primitive vectors are singular")
        b1 = 2.0 * pi * np.array([a2[1], -a2[0]], dtype=float) / area
        b2 = 2.0 * pi * np.array([-a1[1], a1[0]], dtype=float) / area
        return b1, b2

    def frac_to_xy(self, frac: FloatArray) -> FloatArray:
        return frac[0] * self.a1 + frac[1] * self.a2

    def _build_ab_atoms(self) -> list[AtomSite]:
        # Layer 0 uses A,B; layer 1 uses B,C. This is the N=2 slice of ABC stacking.
        frac_pairs = [(self._frac_A, self._frac_B), (self._frac_B, self._frac_C)]
        atoms: list[AtomSite] = []
        for layer, pair in enumerate(frac_pairs):
            z = layer * self.d_A
            for sub, frac in zip(("A", "B"), pair, strict=True):
                xy = self.frac_to_xy(frac)
                atoms.append(AtomSite(layer, sub, np.array([xy[0], xy[1], z], dtype=float)))
        return atoms

    @property
    def positions(self) -> FloatArray:
        return np.vstack([site.r for site in self.atoms])

    def generate_gmkg_path(self, nk_seg: int) -> KPath:
        return self._path_from_points(
            [self.Gamma, self.M, self.K, self.Gamma],
            nk_seg,
            ["Gamma", "M", "K", "Gamma"],
        )

    def generate_gkm_path(self, nk_seg: int, frac_local: float = 1.0) -> KPath:
        points = [
            self.K + frac_local * (self.Gamma - self.K),
            self.K,
            self.K + frac_local * (self.M - self.K),
        ]
        return self._path_from_points(points, nk_seg, ["Gamma", "K", "M"])

    def generate_local_k_mkkp_path(self, nk_seg: int, frac_local: float = 0.1) -> KPath:
        points = [
            self.K + frac_local * (self.M - self.K),
            self.K,
            self.K + frac_local * (self.Kp - self.K),
        ]
        return self._path_from_points(points, nk_seg, ["M", "K", "Kp"])

    def generate_kx_path(
        self,
        nk_seg: int,
        kx_min_Ainv: float,
        kx_max_Ainv: float,
        ky_Ainv: float = 0.0,
    ) -> KPath:
        if nk_seg < 2:
            raise ValueError("nk_seg must be >= 2")
        if kx_min_Ainv >= kx_max_Ainv:
            raise ValueError("kx_min_Ainv must be less than kx_max_Ainv")
        kx = np.linspace(kx_min_Ainv, kx_max_Ainv, nk_seg + 1)
        k_list = self.K + np.column_stack([kx, np.full_like(kx, ky_Ainv)])
        tick_pos = [float(kx[0]), float(kx[len(kx) // 2]), float(kx[-1])]
        return KPath(k_list, kx, tick_pos, [f"{x:g}" for x in tick_pos])

    def generate_local_k_line_path(
        self,
        num_points: int,
        length_Ainv: float,
        direction: str = "kx",
    ) -> KPath:
        if num_points < 2:
            raise ValueError("num_points must be >= 2")
        if length_Ainv <= 0.0:
            raise ValueError("length_Ainv must be positive")

        key = direction.strip().lower()
        if key == "kx":
            unit = np.array([1.0, 0.0], dtype=float)
            labels = [f"-{0.5 * length_Ainv:g}", "K", f"{0.5 * length_Ainv:g}"]
        elif key == "ky":
            unit = np.array([0.0, 1.0], dtype=float)
            labels = [f"-{0.5 * length_Ainv:g}", "K", f"{0.5 * length_Ainv:g}"]
        elif key == "km":
            v = self.M - self.K
            unit = v / np.linalg.norm(v)
            labels = ["M side", "K", "opposite"]
        elif key in {"kgamma", "gamma"}:
            v = self.Gamma - self.K
            unit = v / np.linalg.norm(v)
            labels = ["Gamma side", "K", "opposite"]
        elif key in {"kkp", "kprime", "kp"}:
            v = self.Kp - self.K
            unit = v / np.linalg.norm(v)
            labels = ["K side", "K", "Kp side"]
        else:
            raise ValueError("direction must be kx, ky, km, kgamma, or kkp")

        s = np.linspace(-0.5 * length_Ainv, 0.5 * length_Ainv, num_points)
        k_list = self.K + s[:, None] * unit[None, :]
        k_dist = s + 0.5 * length_Ainv
        tick_pos = [0.0, 0.5 * length_Ainv, length_Ainv]
        return KPath(k_list, k_dist, tick_pos, labels)

    @staticmethod
    def _path_from_points(points: list[FloatArray], nk_seg: int, labels: list[str]) -> KPath:
        if nk_seg < 2:
            raise ValueError("nk_seg must be >= 2")
        if len(points) != len(labels):
            raise ValueError("points and labels must have the same length")

        k_chunks: list[FloatArray] = []
        dist_chunks: list[FloatArray] = []
        tick_pos = [0.0]
        s_acc = 0.0

        for start, end in zip(points[:-1], points[1:], strict=True):
            dk = end - start
            seg_len = float(np.linalg.norm(dk))
            t = np.arange(nk_seg, dtype=float) / nk_seg
            k_chunks.append(start[None, :] + t[:, None] * dk[None, :])
            dist_chunks.append(s_acc + t * seg_len)
            s_acc += seg_len
            tick_pos.append(s_acc)

        k_chunks.append(points[-1][None, :])
        dist_chunks.append(np.array([s_acc], dtype=float))
        return KPath(np.vstack(k_chunks), np.concatenate(dist_chunks), tick_pos, labels)
