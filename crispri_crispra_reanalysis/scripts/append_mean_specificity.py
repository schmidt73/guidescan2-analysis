import sys
import pandas as pd
import numpy as np
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Appends mean specificity column"
    )

    parser.add_argument(
        "mageck_results"
    )

    parser.add_argument(
        "grna_specificities"
    )

    return parser.parse_args()

if __name__ == "__main__":
   args = parse_arguments()

   mageck_results = pd.read_csv(args.mageck_results, index_col=0)
   grna_specificities = pd.read_csv(args.grna_specificities)
   grna_specificities = grna_specificities.groupby(["publication", "crispr_system", "cell_line", "replicate", "gene"]).specificity.mean()

   mageck_results["mean_specificity"] = mageck_results.apply(
       lambda r: grna_specificities[(r["publication"], r["crispr_system"], r["cell_line"], r["replicate"], r["gene"])],
       axis=1
   )

   mageck_results = mageck_results[["publication", "crispr_system", "cell_line", "replicate", "gene",
                                    "mean_specificity", "neg|fdr", "neg|p-value", "neg|lfc", "pos|fdr",
                                    "pos|p-value", "pos|lfc"]]

   mageck_results.sort_values(by=["publication", "crispr_system", "cell_line", "replicate", "gene"]).to_csv(sys.stdout)
