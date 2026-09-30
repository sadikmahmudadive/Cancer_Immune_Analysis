#!/usr/bin/env python3
"""Create continuous Hot/Cold panel scores and figures for one processed TCGA project."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


HOT = ["CD8A", "CD8B", "IFNG", "GZMB", "PRF1", "NKG7", "CXCL9", "CXCL10", "CXCL11", "STAT1", "IRF1", "HLA-A", "HLA-B", "HLA-C", "B2M"]
COLD = ["TGFB1", "CXCL12", "POSTN", "HIF1A", "VEGFA", "KDR", "FAP", "COL1A1", "CD163", "CSF1R", "FOXP3", "IL10"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", choices=["SKCM", "BRCA", "PAAD"], required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    processed = root / "data" / "processed"
    figures = root / "results" / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    expression = pd.read_csv(processed / f"{args.project}_27gene_log2_counts.csv", index_col="gene_symbol")
    metadata = pd.read_csv(processed / f"{args.project}_sample_metadata.csv").set_index("sample_id")
    available_hot = [gene for gene in HOT if gene in expression.index]
    available_cold = [gene for gene in COLD if gene in expression.index]
    standard_deviation = expression.std(axis=1, ddof=0).replace(0, float("nan"))
    z_scores = expression.sub(expression.mean(axis=1), axis=0).div(standard_deviation, axis=0)
    scores = pd.DataFrame({
        "sample_id": expression.columns,
        "Hot_score": z_scores.loc[available_hot].mean(axis=0).values,
        "Cold_score": z_scores.loc[available_cold].mean(axis=0).values,
        "hot_genes_used": len(available_hot),
        "cold_genes_used": len(available_cold),
    }).set_index("sample_id").join(metadata[["case_id", "sample_type", "project", "file_uuid"]])
    scores.to_csv(processed / f"{args.project}_Hot_Cold_scores.csv")

    sns.set_theme(style="whitegrid", context="notebook")
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.4))
    sns.histplot(scores["Hot_score"], bins=30, ax=axes[0], color="#c2410c")
    axes[0].set(title=f"{args.project} Hot score distribution", xlabel="Mean gene Z-score")
    sns.histplot(scores["Cold_score"], bins=30, ax=axes[1], color="#2563eb")
    axes[1].set(title=f"{args.project} Cold score distribution", xlabel="Mean gene Z-score")
    sns.scatterplot(data=scores, x="Hot_score", y="Cold_score", ax=axes[2], color="#475569", s=24, edgecolor="white", linewidth=0.3)
    axes[2].axhline(0, color="black", linewidth=0.7)
    axes[2].axvline(0, color="black", linewidth=0.7)
    axes[2].set(title=f"{args.project} continuous Hot vs Cold scores")
    fig.tight_layout()
    fig.savefig(figures / f"{args.project}_Hot_Cold_score_distributions.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    ordered_samples = scores.sort_values(["Hot_score", "Cold_score"]).index
    row_colors = ["#c2410c"] * len(available_hot) + ["#2563eb"] * len(available_cold)
    grid = sns.clustermap(z_scores.loc[available_hot + available_cold, ordered_samples], row_cluster=False, col_cluster=False, row_colors=row_colors, cmap="vlag", center=0, xticklabels=False, yticklabels=True, figsize=(16, 11), cbar_kws={"label": "Gene Z-score"})
    grid.fig.suptitle(f"TCGA-{args.project} 27-gene panel Z-score heatmap", y=1.02)
    grid.savefig(figures / f"{args.project}_27gene_heatmap.png", dpi=200, bbox_inches="tight")
    plt.close(grid.fig)
    print(f"{args.project}: {expression.shape[0]} genes x {expression.shape[1]} samples; saved continuous Hot/Cold scores and figures.")


if __name__ == "__main__":
    main()
