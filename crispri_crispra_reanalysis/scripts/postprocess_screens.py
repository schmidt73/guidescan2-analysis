import sys

import pandas as pd
import numpy as np
import argparse

from tqdm import tqdm

def parse_arguments():
    parser = argparse.ArgumentParser(
        description=""
    )

    parser.add_argument(
        "preprocessed_grnas_csv", help="CSV file containing grnas."
    )

    parser.add_argument(
        "preprocessed_marson_grnas", help="Marson2022 CSV file containing grnas."
    )

    parser.add_argument(
        "guidescan2_processed", help="Guidescan2 processed gRNAs."
    )

    return parser.parse_args()

if __name__ == "__main__":
   args = parse_arguments()

   grnas_df = pd.read_csv(args.preprocessed_grnas_csv, index_col=0)[[
       "sgRNA_ID", "gene", "sgRNA", "cell_line", "replicate", "publication", "crispr_system"
   ]]
   grnas_df["sgRNA_ID"] = grnas_df["sgRNA_ID"].str.split(",").map(lambda x: x[0])

   marson_df = pd.read_csv(args.preprocessed_marson_grnas)[[
       "sgRNA_ID", "gene", "sgRNA", "cell_line", "replicate", "publication", "crispr_system"
   ]]
   grnas_df = pd.concat([marson_df, grnas_df])

   specificity_df = pd.read_csv(args.guidescan2_processed)[["id", "specificity"]].rename(columns={"id": "sgRNA_ID"})
   out_df = grnas_df.merge(specificity_df, how='left', on='sgRNA_ID')
   out_df = out_df.drop_duplicates(subset=["sgRNA_ID","gene","sgRNA","cell_line","replicate","publication","crispr_system"])
   out_df.to_csv(sys.stdout, index=False)
