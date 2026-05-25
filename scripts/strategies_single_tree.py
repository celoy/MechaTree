"""Reproduce strategies_single_tree.m from the Eloy et al. 2017 Nat Commun lineage.

Two phases:

  * `--build` (one-time): read S3.dat (50-col genome dump from
    ~/Documents/Arbres/FORTRAN/ART_Analysis_LAD/), cluster the population in
    tag-gene space, pick one champion per species, and write the result to
    `data/S3_champions.json`. The JSON is self-contained — it stores
    the schema, the full tag-gene cloud, the species labels, and each
    champion's gene vectors and fitness metrics.
  * default: read the JSON and render the figures. The .dat file is no
    longer touched.

Figures:
  * species scatter — the two final non-coding "tag" genes as x/y, colored by
    detected species, champions starred.
  * 2x2 strategy panels per champion, mirroring the Nat Commun figure
    (b: Safety, c: Photosensitivity, d: Biomass for segments, e: Biomass for
    seeds). Color scales are fixed: Safety 0..4, others 0..1.

Dump schema (this S3.dat prepends 6 position/size cols before the 40-element
genome, so all column indices below are +6 vs the matlab original):
  cols 0..5   physical prefix (positions, sizes)         -- ignored
  cols 6..8   3 angle genes
  cols 9..18  10 NNbranch genes
  cols 19..36 18 NNreserve genes
  cols 37..45 9 reserved genes (last 2 are the species tag)
  col  46     max moment on leaves (fitness)
  col  47     N seeds (fitness, integer)
  col  48     N leaves (integer)
  col  49     extra metric
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Step 21: the 2-cluster k-means + gap-threshold logic now lives in
# ``mechatree.evolution.curate`` so the new tournament runner can share
# this implementation. The script keeps its CLI surface unchanged.
from mechatree.evolution.curate import detect_species, kmeans2  # noqa: F401

DEFAULT_DAT = Path("/Users/Ch/Documents/Arbres/FORTRAN/ART_Analysis_LAD/S3.dat")
DEFAULT_DATA_DIR = Path(__file__).parent.parent / "data"
DEFAULT_JSON = DEFAULT_DATA_DIR / "S3_champions.json"

# 0-based column slices for the dump on disk.
COL_NNBRANCH = slice(9, 19)
COL_NNRESERVE = slice(19, 37)
COL_TAG = slice(44, 46)  # the two final non-coding genes (matlab cols 39-40)
COL_MOMENT = 46
COL_NSEEDS = 47

GRID_N = 128

SAFETY_VRANGE = (0.0, 4.0)
ALLOC_VRANGE = (0.0, 1.0)


# --------- neural-net forward passes (ports of neural_branch.m / neural_reserve.m) ---------


def neural_branch(nb_leaves: np.ndarray, max_stress: np.ndarray, nbranch: np.ndarray) -> np.ndarray:
    """Vectorized port of neural_branch.m. nbranch is a 10-element gene vector."""
    toto = np.tan((nbranch - 0.5) * np.pi * 0.99)
    M1 = toto[0:6].reshape(3, 2, order="F").copy()
    M2 = toto[6:10].reshape(1, 4, order="F").copy()
    M1[0, 0] = 0.0
    M1[2, 1] = 0.0

    X = np.stack([0.01 * nb_leaves, max_stress])
    Z = np.tensordot(M1, X, axes=([1], [0]))
    Zp_hidden = np.tanh(5.0 * Z) / 3.0
    bias = np.full_like(Zp_hidden[:1], 1.0 / 3.0)
    Zp = np.concatenate([Zp_hidden, bias], axis=0)
    F = np.tensordot(M2, Zp, axes=([1], [0]))[0]
    return np.maximum(0.0, F + 1.0)


def neural_reserve(
    nb_leaves: np.ndarray, vol_relative: np.ndarray, nreserve: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Vectorized port of neural_reserve.m. Returns (Pseeds, Pleaves, Phototropism)."""
    toto = np.tan((nreserve - 0.5) * np.pi * 0.99)
    M1 = toto[0:6].reshape(3, 2, order="F").copy()
    M2 = toto[6:18].reshape(3, 4, order="F").copy()
    M1[0, 0] = 0.0
    M1[2, 1] = 0.0

    X = np.stack([0.01 * nb_leaves, vol_relative])
    Z = np.tensordot(M1, X, axes=([1], [0]))
    Zp_hidden = np.tanh(5.0 * Z) / 3.0
    bias = np.full_like(Zp_hidden[:1], 1.0 / 3.0)
    Zp = np.concatenate([Zp_hidden, bias], axis=0)
    F = np.tensordot(M2, Zp, axes=([1], [0]))

    pleaves = np.minimum(np.maximum(0.0, F[0] + 2.0), 4.0) / 4.0
    pseeds = np.minimum(np.maximum(0.0, F[1] + 2.0), 4.0) / 4.0
    phototro = np.minimum(np.maximum(0.0, F[2] + 2.0), 4.0) / 4.0

    s = pseeds + pleaves
    over = s > 1.0
    pleaves = np.where(over, pleaves / np.where(over, s, 1.0), pleaves)
    pseeds = np.where(over, pseeds / np.where(over, s, 1.0), pseeds)
    return pseeds, pleaves, phototro


# --------- build phase: .dat -> JSON ---------


def build_from_dat(dat_path: Path, json_path: Path) -> None:
    data = np.loadtxt(dat_path)
    assert data.shape[1] == 50, f"{dat_path}: expected 50 cols, got {data.shape[1]}"
    n = data.shape[0]
    tag_genes = data[:, COL_TAG]
    moment_leaves = data[:, COL_MOMENT]
    n_seeds = data[:, COL_NSEEDS]

    labels, centroids = detect_species(tag_genes)
    species_ids = sorted(np.unique(labels).tolist())

    print(f"=== {dat_path.name} ({n} individuals) ===")
    print(f"detected {len(species_ids)} species at centroids:")
    for k, c in zip(species_ids, centroids, strict=False):
        print(
            f"  species {k}: tag center=({c[0]:.3f}, {c[1]:.3f}), {(labels == k).sum()} individuals"
        )

    species_payload = []
    for k in species_ids:
        m = labels == k
        masked = np.where(m, moment_leaves, -np.inf)
        rep = int(np.argmax(masked))
        row = data[rep]
        species_payload.append(
            {
                "species_id": int(k),
                "n_members": int(m.sum()),
                "centroid_tag": centroids[k].tolist() if k < len(centroids) else None,
                "champion_index": rep,
                "champion_tag_genes": row[COL_TAG].tolist(),
                "champion_moment_leaves": float(row[COL_MOMENT]),
                "champion_n_seeds": int(row[COL_NSEEDS]),
                "nn_branch": row[COL_NNBRANCH].tolist(),
                "nn_reserve": row[COL_NNRESERVE].tolist(),
                "full_row": row.tolist(),
            }
        )
        print(
            f"  -> species {k} champion: idx {rep}, "
            f"moment={moment_leaves[rep]:.4g}, n_seeds={int(n_seeds[rep])}"
        )

    payload = {
        "dataset": dat_path.name,
        "source_path": str(dat_path),
        "n_individuals": n,
        "schema": {
            "col_nnbranch": [9, 19],
            "col_nnreserve": [19, 37],
            "col_tag": [44, 46],
            "col_moment_leaves": 46,
            "col_n_seeds": 47,
            "note": "0-based python col indices; +6 offset relative to matlab "
            "strategies_single_tree.m because this dump prepends 6 "
            "position/size cols before the 40-element genome.",
        },
        "tag_genes": tag_genes.tolist(),
        "species_labels": labels.tolist(),
        "species": species_payload,
    }
    json_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(payload))
    print(f"wrote {json_path} ({json_path.stat().st_size // 1024} KB)")


# --------- plot phase: JSON -> figures ---------


def species_scatter(
    tag_genes: np.ndarray,
    labels: np.ndarray,
    reps: list[int],
    title: str,
) -> go.Figure:
    palette = ["rgb(51,115,217)", "rgb(51,191,102)"]  # blue, green
    fig = go.Figure()
    for k in np.unique(labels):
        m = labels == k
        fig.add_trace(
            go.Scatter(
                x=tag_genes[m, 0],
                y=tag_genes[m, 1],
                mode="markers",
                marker=dict(color=palette[int(k) % 2], size=5, opacity=0.6),  # noqa: C408
                name=f"species {k} ({m.sum()} indiv)",
            )
        )
    for k, idx in enumerate(reps):
        fig.add_trace(
            go.Scatter(
                x=[tag_genes[idx, 0]],
                y=[tag_genes[idx, 1]],
                mode="markers",
                marker=dict(  # noqa: C408
                    symbol="star",
                    color="red",
                    size=18,
                    line=dict(color="black", width=1),  # noqa: C408
                ),
                name=f"species {k} champion (idx {idx})",
            )
        )
    fig.update_layout(
        title=title,
        xaxis=dict(  # noqa: C408
            title="tag gene 1 (matlab col 39)", range=[0, 1], scaleanchor="y", scaleratio=1
        ),
        yaxis=dict(title="tag gene 2 (matlab col 40)", range=[0, 1]),  # noqa: C408
        width=700,
        height=700,
        paper_bgcolor="white",
        plot_bgcolor="white",
    )
    return fig


def strategy_figure(nbranch: np.ndarray, nreserve: np.ndarray, title: str) -> go.Figure:
    """2x2 strategy panels in the Nat Commun figure's visual style.

    Layout: b Safety, c Photosensitivity, d Biomass for segments, e Biomass for seeds.
    Axis ranges and colour scales match ``strategies_Forest.m``: nb_leaves
    0..100, max_stress 0..0.5, vol_relative 0..2. Safety contours 0..4,
    allocation contours 0..1.
    """
    i = np.arange(GRID_N)[:, None] * np.ones((1, GRID_N))
    j = np.arange(GRID_N)[None, :] * np.ones((GRID_N, 1))
    nb_leaves = 100.0 * i / (GRID_N - 1)
    max_stress = 0.5 * j / (GRID_N - 1)
    vol_relat = 2.0 * j / (GRID_N - 1)

    safety = neural_branch(nb_leaves, max_stress, nbranch)
    pseeds, pleaves, phototro = neural_reserve(nb_leaves, vol_relat, nreserve)

    fig = make_subplots(
        rows=2,
        cols=2,
        subplot_titles=(
            "<b>b</b> Safety",
            "<b>c</b> Photosensitivity",
            "<b>d</b> Biomass for segments",
            "<b>e</b> Biomass for seeds",
        ),
        horizontal_spacing=0.12,
        vertical_spacing=0.12,
    )

    def add_panel(row, col, X, Y, Z, xlabel, ylabel, vrange):
        vmin, vmax = vrange
        # Distinct colorbars per panel.
        cbar_x = 0.48 if col == 1 else 1.02
        cbar_y = 0.79 if row == 1 else 0.21
        fig.add_trace(
            go.Contour(
                x=X[0, :],
                y=Y[:, 0],
                z=Z,
                contours=dict(start=vmin, end=vmax, size=(vmax - vmin) / 15),  # noqa: C408
                line=dict(color="black", width=0.4),  # noqa: C408
                colorscale="Viridis",
                zmin=vmin,
                zmax=vmax,
                colorbar=dict(x=cbar_x, y=cbar_y, len=0.42, thickness=12),  # noqa: C408
                showscale=True,
            ),
            row=row,
            col=col,
        )
        fig.update_xaxes(title_text=xlabel, row=row, col=col)
        fig.update_yaxes(title_text=ylabel, row=row, col=col)

    add_panel(
        1, 1, nb_leaves, max_stress, safety, "number of foliages", "relative stress", SAFETY_VRANGE
    )
    add_panel(
        1,
        2,
        nb_leaves,
        vol_relat,
        phototro,
        "total number of foliages",
        "relative biomass volume",
        ALLOC_VRANGE,
    )
    add_panel(
        2,
        1,
        nb_leaves,
        vol_relat,
        pleaves,
        "total number of foliages",
        "relative biomass volume",
        ALLOC_VRANGE,
    )
    add_panel(
        2,
        2,
        nb_leaves,
        vol_relat,
        pseeds,
        "total number of foliages",
        "relative biomass volume",
        ALLOC_VRANGE,
    )

    fig.update_layout(
        title_text=title,
        width=950,
        height=850,
        paper_bgcolor="white",
        plot_bgcolor="white",
        showlegend=False,
    )
    return fig


def plot_from_json(json_path: Path, save_dir: Path | None) -> list[go.Figure]:
    payload = json.loads(json_path.read_text())
    dataset = payload["dataset"]
    tag_genes = np.array(payload["tag_genes"])
    labels = np.array(payload["species_labels"])
    species = payload["species"]
    reps = [s["champion_index"] for s in species]

    print(
        f"=== {dataset} ({payload['n_individuals']} individuals, "
        f"{len(species)} species) — from {json_path.name} ==="
    )
    for s in species:
        print(
            f"  species {s['species_id']}: idx {s['champion_index']}, "
            f"moment={s['champion_moment_leaves']:.4g}, n_seeds={s['champion_n_seeds']}"
        )

    figs: list[go.Figure] = []
    stem = Path(dataset).stem

    fig_s = species_scatter(
        tag_genes,
        labels,
        reps,
        f"{dataset}: species scatter "
        f"({payload['n_individuals']} individuals, {len(species)} species)",
    )
    figs.append(fig_s)
    if save_dir is not None:
        save_dir.mkdir(parents=True, exist_ok=True)
        fig_s.write_image(str(save_dir / f"{stem}_species.png"))

    for s in species:
        nbranch = np.array(s["nn_branch"])
        nreserve = np.array(s["nn_reserve"])
        fig = strategy_figure(
            nbranch,
            nreserve,
            f"{dataset} — species {s['species_id']} champion "
            f"(idx {s['champion_index']}, moment={s['champion_moment_leaves']:.3g}, "
            f"n_seeds={s['champion_n_seeds']})",
        )
        figs.append(fig)
        if save_dir is not None:
            fig.write_image(str(save_dir / f"{stem}_species{s['species_id']}.png"))

    if save_dir is not None:
        print(f"saved {len(figs)} PNGs under {save_dir}")
    return figs


# --------- CLI ---------


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument(
        "--json",
        type=Path,
        default=DEFAULT_JSON,
        help=f"champion JSON to plot from (default: {DEFAULT_JSON}).",
    )
    p.add_argument(
        "--build", action="store_true", help="(re)build the JSON from --dat before plotting."
    )
    p.add_argument(
        "--dat",
        type=Path,
        default=DEFAULT_DAT,
        help=f"source .dat file when --build is set (default: {DEFAULT_DAT}).",
    )
    p.add_argument(
        "--save-dir",
        type=Path,
        default=None,
        help="if given, save PNGs here in addition to showing.",
    )
    p.add_argument(
        "--no-show",
        action="store_true",
        help="skip opening the plotly figures (useful with --save-dir for headless runs).",
    )
    args = p.parse_args()

    if args.build or not args.json.exists():
        if not args.json.exists() and not args.build:
            print(f"{args.json} not found — building from {args.dat}")
        build_from_dat(args.dat, args.json)

    figs = plot_from_json(args.json, args.save_dir)

    if not args.no_show:
        for fig in figs:
            fig.show()


if __name__ == "__main__":
    main()
