import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import scipy.stats as stats
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Plot specificity vs hits non-hits."
    )

    parser.add_argument(
        "processed_screen_data"
    )

    parser.add_argument(
        "mageck_summary"
    )

    return parser.parse_args()

if __name__ == "__main__":
   args = parse_arguments()

   screen_df = pd.read_csv(args.processed_screen_data)
   screen_df = screen_df[~screen_df.specificity.isna()]
   mageck_df = pd.read_csv(args.mageck_summary)

   fdrs = [0.05, 0.10]
   pvalues = [0.001, 0.005, 0.01, 0.05, 0.10]
   for fdr_threshold in fdrs:
       for pvalue_threshold in pvalues:
            mageck_groups = mageck_df.groupby(["publication", "crispr_system", "cell_line"])
            for identifier, df in mageck_groups:
                publication, crispr_system, cell_line = identifier
                print(publication)
                print(crispr_system)

                pos_hit_genes = set(df[(df["pos|p-value"] < fdr_threshold) & (df["pos|fdr"] < pvalue_threshold)]["gene"])
                neg_hit_genes = set(df[(df["neg|p-value"] < fdr_threshold) & (df["neg|fdr"] < pvalue_threshold)]["gene"])

                cond1 = (screen_df["publication"] == publication) & (screen_df["cell_line"] == cell_line)
                cond2 = (screen_df["crispr_system"] == crispr_system)

                pos_hit_df = screen_df[cond1 & cond2 & (screen_df["gene"].isin(pos_hit_genes))]
                neg_hit_df = screen_df[cond1 & cond2 & (screen_df["gene"].isin(neg_hit_genes))]
                pos_non_hit_df = screen_df[cond1 & cond2 & ~(screen_df["gene"].isin(pos_hit_genes))]
                neg_non_hit_df = screen_df[cond1 & cond2 & ~(screen_df["gene"].isin(neg_hit_genes))]

                mean_gene_spec_pos_hit_df = pos_hit_df.groupby('gene').specificity.mean()
                mean_gene_spec_neg_hit_df = neg_hit_df.groupby('gene').specificity.mean()
                mean_gene_spec_pos_non_hit_df = pos_non_hit_df.groupby('gene').specificity.mean()
                mean_gene_spec_neg_non_hit_df = neg_non_hit_df.groupby('gene').specificity.mean()

                fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(8, 4.5))

                if len(neg_hit_df) > 2:
                    res = stats.mannwhitneyu(mean_gene_spec_neg_hit_df, mean_gene_spec_neg_non_hit_df, method='asymptotic')

                    ax = axes[0]
                    sns.ecdfplot(data=mean_gene_spec_neg_hit_df, ax=ax, label=f"Negative hits (n={len(mean_gene_spec_neg_hit_df)})")
                    sns.ecdfplot(data=mean_gene_spec_neg_non_hit_df, ax=ax, label=f"Non-hits (n={len(mean_gene_spec_neg_non_hit_df)})")
                    ax.set_xlim(0, 1)

                    if res.pvalue < 2.22e-16:
                        ax.text(0.55, 0.25,f"p-value < 2.22e-16", transform=ax.transAxes)
                    else:
                        ax.text(0.55, 0.25,f"p-value = {res.pvalue:.6f}", transform=ax.transAxes)

                    ax.legend(loc='upper left')
                    ax.set_xlabel('Mean Gene Specificity')
                    ax.set_ylabel('F(x)')

                if len(pos_hit_df) > 2:
                    res = stats.mannwhitneyu(mean_gene_spec_pos_hit_df, mean_gene_spec_pos_non_hit_df, method='asymptotic')

                    ax = axes[1]
                    sns.ecdfplot(data=mean_gene_spec_pos_hit_df, ax=ax, label=f"Positive hits (n={len(mean_gene_spec_pos_hit_df)})")
                    sns.ecdfplot(data=mean_gene_spec_pos_non_hit_df, ax=ax, label=f"Non-hits (n={len(mean_gene_spec_pos_non_hit_df)})")

                    if res.pvalue < 10e-6:
                        ax.text(0.55, 0.25,f"p-value < 2.22e-16", transform=ax.transAxes)
                    else:
                        ax.text(0.55, 0.25,f"p-value = {res.pvalue:.6f}", transform=ax.transAxes)

                    ax.legend(loc='upper left')
                    ax.set_xlabel('Mean Gene Specificity')
                    ax.set_ylabel(None)
                    ax.set_yticks([])

                fig.tight_layout()

                os.makedirs(f"figures/pval_{pvalue_threshold}_fdr_{fdr_threshold}", exist_ok=True)
                fig.savefig(f"figures/pval_{pvalue_threshold}_fdr_{fdr_threshold}/{publication}_{crispr_system}_{cell_line}_specificity.nocontrols.pdf")
