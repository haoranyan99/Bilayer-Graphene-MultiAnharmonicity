from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np

from config import load_config
from model_factory import build_model, build_structure, choose_kpoint, read_onsite_layer_gap_eV
from structure import BilayerGrapheneStructure, KPath


def nice_polyline(points: list[tuple[float, float]]) -> str:
    return " ".join(f"{x:.3f},{y:.3f}" for x, y in points)


def write_svg_band_plot(
    out_path: Path,
    k_dist: np.ndarray,
    bands: np.ndarray,
    tick_pos: list[float],
    tick_labels: list[str],
    ylim: tuple[float, float] | None,
) -> None:
    width = 760
    height = 520
    left = 82
    right = 26
    top = 42
    bottom = 72
    plot_w = width - left - right
    plot_h = height - top - bottom

    xmin = float(k_dist[0])
    xmax = float(k_dist[-1])
    if ylim is None:
        ymin = float(np.min(bands))
        ymax = float(np.max(bands))
        pad = 0.05 * max(ymax - ymin, 1.0)
        ymin -= pad
        ymax += pad
    else:
        ymin, ymax = ylim

    def xmap(x: float) -> float:
        return left + (x - xmin) / (xmax - xmin) * plot_w

    def ymap(y: float) -> float:
        return top + (ymax - y) / (ymax - ymin) * plot_h

    lines: list[str] = []
    lines.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">')
    lines.append('<rect width="100%" height="100%" fill="white"/>')
    lines.append(f'<text x="{width / 2:.1f}" y="24" text-anchor="middle" font-family="Arial" font-size="18">AB bilayer graphene SK bands</text>')
    lines.append(f'<rect x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" fill="none" stroke="black" stroke-width="1"/>')

    for xpos, label in zip(tick_pos, tick_labels, strict=True):
        x = xmap(float(xpos))
        lines.append(f'<line x1="{x:.3f}" y1="{top}" x2="{x:.3f}" y2="{top + plot_h}" stroke="#d0d0d0" stroke-width="1"/>')
        lines.append(f'<text x="{x:.3f}" y="{height - 38}" text-anchor="middle" font-family="Arial" font-size="14">{label}</text>')

    if ymin < 0.0 < ymax:
        y0 = ymap(0.0)
        lines.append(f'<line x1="{left}" y1="{y0:.3f}" x2="{left + plot_w}" y2="{y0:.3f}" stroke="#999" stroke-width="1" stroke-dasharray="5,5"/>')

    for frac in np.linspace(0.0, 1.0, 6):
        energy = ymin + frac * (ymax - ymin)
        y = ymap(float(energy))
        lines.append(f'<line x1="{left - 5}" y1="{y:.3f}" x2="{left}" y2="{y:.3f}" stroke="black" stroke-width="1"/>')
        lines.append(f'<text x="{left - 10}" y="{y + 4:.3f}" text-anchor="end" font-family="Arial" font-size="12">{energy:.2f}</text>')

    for ib in range(bands.shape[1]):
        pts = [(xmap(float(k_dist[i])), ymap(float(bands[i, ib]))) for i in range(len(k_dist))]
        lines.append(f'<polyline points="{nice_polyline(pts)}" fill="none" stroke="black" stroke-width="1.4"/>')

    lines.append(f'<text x="{left + plot_w / 2:.1f}" y="{height - 14}" text-anchor="middle" font-family="Arial" font-size="14">k path</text>')
    lines.append(f'<text x="18" y="{top + plot_h / 2:.1f}" transform="rotate(-90 18 {top + plot_h / 2:.1f})" text-anchor="middle" font-family="Arial" font-size="14">Energy (eV)</text>')
    lines.append("</svg>")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines), encoding="utf-8")


def build_full_kpath(st: BilayerGrapheneStructure, plot_cfg: dict, cli_path: str | None) -> KPath:
    path_name = cli_path or str(plot_cfg.get("path", "GMKG"))
    if path_name == "GMKG":
        return st.generate_gmkg_path(int(plot_cfg.get("Nk_seg", 200)))
    if path_name == "GKM":
        return st.generate_gkm_path(
            int(plot_cfg.get("Nk_seg", 200)),
            float(plot_cfg.get("frac_local", 1.0)),
        )
    if path_name == "localK_MKKp":
        return st.generate_local_k_mkkp_path(
            int(plot_cfg.get("Nk_seg", 200)),
            float(plot_cfg.get("frac_local", 0.1)),
        )
    raise ValueError("full path must be GMKG, GKM, or localK_MKKp")


def build_local_kpath(st: BilayerGrapheneStructure, plot_cfg: dict) -> KPath:
    return st.generate_local_k_line_path(
        int(plot_cfg.get("num_points", 401)),
        float(plot_cfg.get("length_Ainv", 0.2)),
        str(plot_cfg.get("direction", "kx")),
    )


def compact_float_tag(x: float) -> str:
    if abs(x) < 5e-13:
        x = 0.0
    return f"{x:.10g}".replace("-", "m").replace(".", "p")


def should_shift_to_midgap(config: dict, plot_cfg: dict) -> bool:
    if "energy_shift_to_midgap" in plot_cfg:
        return bool(plot_cfg["energy_shift_to_midgap"])
    return bool(config.get("patch", {}).get("energy_shift_to_midgap", False))


def band_energy_shift_eV(st: BilayerGrapheneStructure, model, plot_cfg: dict) -> float:
    ref_name = str(plot_cfg.get("energy_zero_point", plot_cfg.get("center", "K")))
    ref_k = choose_kpoint(st, ref_name)
    evals = np.linalg.eigvalsh(model.build_hk(ref_k))
    if len(evals) < 2:
        return 0.0
    b0 = len(evals) // 2 - 1
    b1 = len(evals) // 2
    return 0.5 * float(evals[b0] + evals[b1])


def output_path_from_config(
    plot_kind: str,
    plot_cfg: dict,
    cli_path: str | None,
    output_override: str | None,
) -> Path:
    if output_override:
        return Path(output_override)

    output_dir = Path(str(plot_cfg.get("output_dir", "figures")))
    model_name = str(plot_cfg.get("model", "sk")).lower()
    gap_tag = compact_float_tag(read_onsite_layer_gap_eV(plot_cfg))

    if plot_kind == "full":
        path_tag = cli_path or str(plot_cfg.get("path", "GMKG"))
        if "Nk_seg" in plot_cfg:
            path_tag += f"_Nk{int(plot_cfg['Nk_seg'])}"
    else:
        direction = str(plot_cfg.get("direction", "kx"))
        length_tag = compact_float_tag(float(plot_cfg.get("length_Ainv", 0.2)))
        npt = int(plot_cfg.get("num_points", 401))
        path_tag = f"localK_{direction}_L{length_tag}_N{npt}"

    return output_dir / f"blg_band_{model_name}_gap{gap_tag}_{path_tag}.svg"


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot AB bilayer graphene SK bands as SVG.")
    parser.add_argument("-c", "--config", default="configs/bilayer_sk.json", help="JSON config path")
    parser.add_argument("--plot", default="localK", choices=["full", "localK"], help="which config section to use")
    parser.add_argument("-o", "--output", default=None, help="full output SVG path; overrides automatic naming")
    parser.add_argument("--path", default=None, choices=["GMKG", "GKM", "localK_MKKp"], help="full-path name; overrides band_plot_full.path")
    parser.add_argument("--nk-seg", type=int, default=None, help="k points per segment; overrides config")
    parser.add_argument("--num-points", type=int, default=None, help="total points for localK_line; overrides config")
    parser.add_argument("--length-Ainv", type=float, default=None, help="total localK_line length in 1/A; overrides config")
    parser.add_argument("--direction", default=None, choices=["kx", "ky", "km", "kgamma", "kkp"], help="localK_line direction; overrides config")
    parser.add_argument("--frac-local", type=float, default=None, help="local path fraction; overrides config")
    parser.add_argument("--ylim", nargs=2, type=float, default=None, metavar=("EMIN", "EMAX"))
    args = parser.parse_args()

    config = load_config(args.config)
    section_name = "band_plot_full" if args.plot == "full" else "band_plot_localK"
    band_cfg = dict(config.get(section_name, {}))
    if args.nk_seg is not None:
        band_cfg["Nk_seg"] = args.nk_seg
    if args.num_points is not None:
        band_cfg["num_points"] = args.num_points
    if args.length_Ainv is not None:
        band_cfg["length_Ainv"] = args.length_Ainv
    if args.direction is not None:
        band_cfg["direction"] = args.direction
    if args.frac_local is not None:
        band_cfg["frac_local"] = args.frac_local

    st = build_structure(config)
    fallback_h_cfg = config.get("hamiltonian", {})
    model = build_model(config, band_cfg or fallback_h_cfg, st)

    if args.plot == "full":
        kpath = build_full_kpath(st, band_cfg, args.path)
    else:
        kpath = build_local_kpath(st, band_cfg)
    bands = model.bands_along_path(kpath.k_list)
    energy_shift = 0.0
    if should_shift_to_midgap(config, band_cfg):
        energy_shift = band_energy_shift_eV(st, model, band_cfg)
        bands = bands - energy_shift

    out_path = output_path_from_config(args.plot, band_cfg, args.path, args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    cfg_ylim = band_cfg.get("ylim", None)
    ylim = tuple(args.ylim) if args.ylim is not None else (tuple(cfg_ylim) if cfg_ylim is not None else None)
    write_svg_band_plot(out_path, kpath.k_dist, bands, kpath.tick_pos, kpath.tick_labels, ylim)

    print(f"Wrote: {out_path}")
    print(f"kpoints: {len(kpath.k_dist)}")
    print(f"energy_shift_to_midgap: {should_shift_to_midgap(config, band_cfg)}")
    print(f"energy_shift_eV: {energy_shift:.12g}")
    print(f"energy span: [{bands.min():.8g}, {bands.max():.8g}] eV")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
