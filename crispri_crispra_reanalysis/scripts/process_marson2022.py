import sys

import pandas as pd
import numpy as np
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Processes Marson2022 data to kmers and CSV file."
    )

    return parser.parse_args()

if __name__ == "__main__":
   args = parse_arguments()

   crispra_ifng = pd.read_csv("data/inputs/marson2022/Marson2022_CRISPRa_IFNG_1.sgrna_summary.txt", sep="\t")
   crispra_il2 = pd.read_csv("data/inputs/marson2022/Marson2022_CRISPRa_IL2_1.sgrna_summary.txt", sep="\t")
   crispri_ifng = pd.read_csv("data/inputs/marson2022/Marson2022_CRISPRi_IFNG_1.sgrna_summary.txt", sep="\t")
   crispri_il2 = pd.read_csv("data/inputs/marson2022/Marson2022_CRISPRi_IL2_1.sgrna_summary.txt", sep="\t")

   crispra_ifng.rename(columns={"sgrna": "sgRNA_ID"}, inplace=True)
   crispra_il2.rename(columns={"sgrna": "sgRNA_ID"}, inplace=True)
   crispri_ifng.rename(columns={"sgrna": "sgRNA_ID"}, inplace=True)
   crispri_il2.rename(columns={"sgrna": "sgRNA_ID"}, inplace=True)

   crispra_ifng["sgRNA"] = crispra_ifng.sgRNA_ID.str.split("_").map(lambda x: x[1])
   crispra_il2["sgRNA"] = crispra_il2.sgRNA_ID.str.split("_").map(lambda x: x[1])
   crispri_ifng["sgRNA"] = crispri_ifng.sgRNA_ID.str.split("_").map(lambda x: x[1])
   crispri_il2["sgRNA"] = crispri_il2.sgRNA_ID.str.split("_").map(lambda x: x[1])

   crispra_ifng["cell_line"] = "IFNG"
   crispra_il2["cell_line"] = "IL2"
   crispri_ifng["cell_line"] = "IFNG"
   crispri_il2["cell_line"] = "IL2"

   crispra_ifng["replicate"] = 1
   crispra_il2["replicate"] = 1
   crispri_ifng["replicate"] = 1
   crispri_il2["replicate"] = 1

   crispra_ifng["crispr_system"] = "CRISPRa"
   crispra_il2["crispr_system"] = "CRISPRa"
   crispri_ifng["crispr_system"] = "CRISPRi"
   crispri_il2["crispr_system"] = "CRISPRi"

   df = pd.concat([crispra_ifng, crispra_il2, crispri_il2, crispri_ifng]).rename(columns={"Gene": "gene"})
   df["publication"] = "Marson2022"
   df[["sgRNA_ID", "sgRNA", "gene", "cell_line", "replicate", "crispr_system", "publication"]].to_csv(
       "data/preprocessed_marson2022_grnas.csv",
       index=False
   )

   guidescan_data = df.rename(columns={"sgRNA_ID": "id", "sgRNA": "sequence"})[["id", "sequence"]]
   guidescan_data["PAM"] = "NGG"
   guidescan_data["chromosome"] = "1"
   guidescan_data["position"] = "1"
   guidescan_data["sense"] = "+"

   guidescan_data.to_csv(
       "data/preprocessed_marson2022_kmers.csv",
       index=False
   )
