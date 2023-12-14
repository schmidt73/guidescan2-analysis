import os.path
import sys

import pandas as pd
import numpy as np
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Summarize MAGeCK results into CSV."
    )

    parser.add_argument(
        "--controls", default=False, action="store_true",
    )

    parser.add_argument(
        "mageck_gene_summaries", nargs="+"
    )

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()

    dfs = []
    for gene_summary_file in args.mageck_gene_summaries:
        name = os.path.basename(gene_summary_file)

        if args.controls:
            publication, crispr_system, cell_line, replicate = name[:-len(".gene_summary.txt")].split("_")
        else:
            publication, crispr_system, cell_line, replicate = name[:-len(".nocontrols.gene_summary.txt")].split("_")

        gene_summary = pd.read_csv(gene_summary_file, sep='\t')

        # neg_hits = (gene_summary['neg|fdr'] < 0.05) & (gene_summary['neg|p-value'] < 0.01)
        # pos_hits = (gene_summary['pos|fdr'] < 0.05) & (gene_summary['pos|p-value'] < 0.01)

        hits_df = gene_summary[['id', 'neg|fdr', 'neg|p-value', 'neg|lfc', 'pos|fdr', 'pos|p-value', 'pos|lfc']]
        hits_df['publication'] = publication
        hits_df['crispr_system'] = crispr_system
        hits_df['cell_line'] = cell_line
        hits_df['replicate'] = replicate

        dfs.append(hits_df)

    pd.concat(dfs).rename(columns={"id": "gene"}).to_csv(sys.stdout)
