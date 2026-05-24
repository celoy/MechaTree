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

import matplotlib.pyplot as plt
import numpy as np

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


# --------- species detection ---------


def kmeans2(points: np.ndarray, n_iter: int = 50) -> tuple[np.ndarray, np.ndarray]:
    """Tiny 2-cluster k-means with farthest-point init. Returns (labels, centroids)."""
    c0 = points[0]
    c1 = points[np.argmax(np.linalg.norm(points - c0, axis=1))]
    centroids = np.stack([c0, c1])
    labels = np.zeros(len(points), dtype=int)
    for _ in range(n_iter):
        d = np.linalg.norm(points[:, None, :] - centroids[None, :, :], axis=2)
        new_labels = d.argmin(axis=1)
        if np.array_equal(new_labels, labels):
            break
        labels = new_labels
        for k in (0, 1):
            members = points[labels == k]
            if len(members):
                centroids[k] = members.mean(axis=0)
    return labels, centroids


def detect_species(
    tag_genes: np.ndarray, gap_threshold: float = 0.15
) -> tuple[np.ndarray, np.ndarray]:
    """Cluster into 1 or 2 species. Collapses to 1 if the centroid gap is below threshold."""
    labels, centroids = kmeans2(tag_genes)
    if np.linalg.norm(centroids[0] - centroids[1]) < gap_threshold:
        return np.zeros(len(tag_genes), dtype=int), tag_genes.mean(axis=0, keepdims=True)
    return labels, centroids


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
) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(7, 7))
    palette = np.array([[0.20, 0.45, 0.85], [0.20, 0.75, 0.40]])  # blue, green
    for k in np.unique(labels):
        m = labels == k
        ax.scatter(
            tag_genes[m, 0],
            tag_genes[m, 1],
            c=[palette[k % 2]],
            s=8,
            alpha=0.6,
            edgecolors="none",
            label=f"species {k} ({m.sum()} indiv)",
        )
    for k, idx in enumerate(reps):
        ax.scatter(
            tag_genes[idx, 0],
            tag_genes[idx, 1],
            marker="*",
            c="red",
            s=200,
            edgecolors="black",
            linewidths=0.5,
            label=f"species {k} champion (idx {idx})",
        )
    ax.set_xlabel("tag gene 1 (matlab col 39)")
    ax.set_ylabel("tag gene 2 (matlab col 40)")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.set_title(title)
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()
    return fig


def strategy_figure(nbranch: np.ndarray, nreserve: np.ndarray, title: str) -> plt.Figure:
    """2x2 strategy panels in the Nat Commun figure's visual style.

    Layout: b Safety, c Photosensitivity, d Biomass for segments, e Biomass for seeds.
    Axis ranges and color scales match `strategies_Forest.m` from
    ART_Analysis_LAD/ART_Strategies/ (the variant whose plots match the paper
    figure): nb_leaves 0..100, max_stress 0..0.5, vol_relative 0..2.
    Safety contours 0..4, allocation contours 0..1.
    """
    i = np.arange(GRID_N)[:, None] * np.ones((1, GRID_N))
    j = np.arange(GRID_N)[None, :] * np.ones((GRID_N, 1))
    nb_leaves = 100.0 * i / (GRID_N - 1)
    max_stress = 0.5 * j / (GRID_N - 1)
    vol_relat = 2.0 * j / (GRID_N - 1)

    safety = neural_branch(nb_leaves, max_stress, nbranch)
    pseeds, pleaves, phototro = neural_reserve(nb_leaves, vol_relat, nreserve)

    fig, axes = plt.subplots(2, 2, figsize=(9, 7.5))
    fig.suptitle(title)

    def panel(ax, X, Y, Z, xlabel, ylabel, sub_title, letter, vrange):
        vmin, vmax = vrange
        levels = np.linspace(vmin, vmax, 16)
        cs = ax.contourf(
            X, Y, Z, levels=levels, cmap="viridis", vmin=vmin, vmax=vmax, extend="both"
        )
        ax.contour(X, Y, Z, levels=levels, colors="k", linewidths=0.4)
        fig.colorbar(cs, ax=ax, ticks=np.linspace(vmin, vmax, 5))
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(sub_title)
        ax.text(
            -0.18, 1.05, letter, transform=ax.transAxes, fontsize=14, fontweight="bold", va="top"
        )

    panel(
        axes[0, 0],
        nb_leaves,
        max_stress,
        safety,
        "number of foliages",
        "relative stress",
        "Safety",
        "b",
        SAFETY_VRANGE,
    )
    panel(
        axes[0, 1],
        nb_leaves,
        vol_relat,
        phototro,
        "total number of foliages",
        "relative biomass volume",
        "Photosensitivity",
        "c",
        ALLOC_VRANGE,
    )
    panel(
        axes[1, 0],
        nb_leaves,
        vol_relat,
        pleaves,
        "total number of foliages",
        "relative biomass volume",
        "Biomass for segments",
        "d",
        ALLOC_VRANGE,
    )
    panel(
        axes[1, 1],
        nb_leaves,
        vol_relat,
        pseeds,
        "total number of foliages",
        "relative biomass volume",
        "Biomass for seeds",
        "e",
        ALLOC_VRANGE,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    return fig


def plot_from_json(json_path: Path, save_dir: Path | None) -> list[plt.Figure]:
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

    figs: list[plt.Figure] = []
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
        fig_s.savefig(save_dir / f"{stem}_species.png", dpi=120)

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
            fig.savefig(save_dir / f"{stem}_species{s['species_id']}.png", dpi=120)

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
        help="skip plt.show() (useful with --save-dir for headless runs).",
    )
    args = p.parse_args()

    if args.build or not args.json.exists():
        if not args.json.exists() and not args.build:
            print(f"{args.json} not found — building from {args.dat}")
        build_from_dat(args.dat, args.json)

    plot_from_json(args.json, args.save_dir)

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
