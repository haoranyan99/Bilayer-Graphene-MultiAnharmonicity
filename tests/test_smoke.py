from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np

from sk_model import BilayerSKModel
from tb_model import BilayerTBModel
from structure import BilayerGrapheneStructure


def main() -> int:
    st = BilayerGrapheneStructure()
    model = BilayerSKModel(st)
    assert len(st.atoms) == 4
    assert len(model.hoppings) > 0

    hk = model.build_hk(st.K)
    assert hk.shape == (4, 4)
    assert np.allclose(hk, hk.conj().T)

    bands = model.bands_at_k(st.K)
    assert bands.shape == (4,)
    assert np.all(np.isfinite(bands))

    gapped_model = BilayerSKModel(st, onsite_layer_gap_eV=0.02)
    assert np.allclose(
        np.diag(gapped_model.build_hk(st.Gamma)).real,
        np.diag(model.build_hk(st.Gamma)).real + np.array([-0.01, -0.01, 0.01, 0.01]),
    )

    path = st.generate_gmkg_path(5)
    path_bands = model.bands_along_path(path.k_list)
    assert path_bands.shape == (16, 4)

    tb_model = BilayerTBModel(st)
    tb_hk = tb_model.build_hk(st.K)
    assert tb_hk.shape == (4, 4)
    assert np.allclose(tb_hk, tb_hk.conj().T)
    assert np.all(np.isfinite(tb_model.bands_at_k(st.K)))

    print("smoke test passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
