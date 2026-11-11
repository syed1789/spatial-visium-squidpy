#!/usr/bin/env python3
"""Visium spatial transcriptomics with Scanpy + Squidpy.

Dataset: 10x Genomics Visium, adult mouse brain (H&E), as shipped by
`squidpy.datasets.visium_hne_adata()` (Squidpy tutorial dataset, ~700 spots).

Steps: QC -> normalize/HVG -> PCA/neighbors/Leiden -> spatial graph ->
neighborhood enrichment -> co-occurrence -> Moran's I -> ligand-receptor.
All figures and tables are written to --outdir.
"""
import argparse
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc
import squidpy as sq

warnings.filterwarnings("ignore")


def main(outdir: Path, seed: int = 0) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    figdir = outdir / "figures"
    figdir.mkdir(exist_ok=True)
    sc.settings.figdir = figdir
    sc.settings.set_figure_params(dpi=110, facecolor="white")

    def save(name):
        plt.savefig(figdir / f"{name}.png", dpi=130, bbox_inches="tight")
        plt.close("all")

    adata = sq.datasets.visium_hne_adata()
    adata.var_names_make_unique()
    print(f"Loaded: {adata.n_obs} spots x {adata.n_vars} genes")

    # --- QC ---
    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    sc.pl.spatial(adata, color=["total_counts", "n_genes_by_counts"], show=False)
    save("01_qc_spatial")
    sc.pp.filter_genes(adata, min_cells=10)

    # --- Normalize, HVG, embedding, clustering ---
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    sc.pp.pca(adata, random_state=seed)
    sc.pp.neighbors(adata, random_state=seed)
    sc.tl.umap(adata, random_state=seed)
    adata.uns.pop("leiden_colors", None)  # dataset ships stale colors from its own clustering
    sc.tl.leiden(adata, resolution=0.6, random_state=seed,
                 flavor="igraph", n_iterations=2, directed=False)
    sc.pl.umap(adata, color="leiden", show=False)
    save("02_umap_leiden")
    sc.pl.spatial(adata, color="leiden", show=False)
    save("03_spatial_leiden")

    # --- Spatial statistics ---
    sq.gr.spatial_neighbors(adata, coord_type="grid", n_neighs=6)

    sq.gr.nhood_enrichment(adata, cluster_key="leiden", seed=seed)
    sq.pl.nhood_enrichment(adata, cluster_key="leiden", method="average")
    save("04_nhood_enrichment")

    sq.gr.co_occurrence(adata, cluster_key="leiden")
    sq.pl.co_occurrence(adata, cluster_key="leiden", clusters="0")
    save("05_co_occurrence")

    sq.gr.spatial_autocorr(adata, mode="moran", n_perms=100, n_jobs=1, seed=seed)
    moran = adata.uns["moranI"]
    moran.head(50).to_csv(outdir / "moran_top50.csv")
    top = moran.head(6).index.tolist()
    sc.pl.spatial(adata, color=top, ncols=3, show=False)
    save("06_top_spatially_variable_genes")

    # --- Ligand-receptor (CellPhoneDB-style, via Omnipath) ---
    # Needs network access to fetch the Omnipath resource; skipped if unavailable.
    try:
        sq.gr.ligrec(adata, n_perms=100, cluster_key="leiden",
                     use_raw=False, transmitter_params={"categories": "ligand"},
                     receiver_params={"categories": "receptor"}, seed=seed)
        res = adata.uns["leiden_ligrec"]
        res["means"].to_csv(outdir / "ligrec_means.csv")
        res["pvalues"].to_csv(outdir / "ligrec_pvalues.csv")
        # Plot only strong, significant pairs between a few clusters; plotting
        # every pair produces an image ~170k px wide.
        sq.pl.ligrec(adata, cluster_key="leiden", source_groups="0",
                     target_groups=["1", "2", "8"], pvalue_threshold=0.01,
                     means_range=(1.3, float("inf")), alpha=0.01, swap_axes=True)
        save("07_ligand_receptor")
    except Exception as e:  # noqa: BLE001
        print(f"ligrec skipped: {type(e).__name__}: {e}")

    adata.obs[["leiden"]].to_csv(outdir / "cluster_assignments.csv")
    print("Cluster sizes:\n", adata.obs["leiden"].value_counts().sort_index())
    print("Top Moran's I genes:\n", moran.head(10)[["I", "pval_norm_fdr_bh"]])


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", type=Path, default=Path("results"))
    p.add_argument("--seed", type=int, default=0)
    a = p.parse_args()
    main(a.outdir, a.seed)
