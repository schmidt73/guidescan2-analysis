import pandas as pd
import numpy as np
import argparse
import sys

def parse_arguments():
    parser = argparse.ArgumentParser(
        description=("Pre-processes screen data into a format suitable" + 
                     "for MAGeCK and Guidescan2")
    )

    return parser.parse_args()

def preprocess_doench2018():
    sheets = pd.read_excel("data/inputs/doench2018_crispri_counts.xlsx", sheet_name=None)
    def process_doench2018_sheet(counts, grnas):
        counts.columns = pd.MultiIndex.from_arrays([counts.columns, counts.iloc[0].values])
        counts.drop([0], inplace=True)
        counts.rename(columns={"Unnamed: 0": "sgRNA", "sgRNA Sequence": ""}, inplace=True)
        counts = pd.melt(counts, id_vars=["sgRNA", ("pDNA", "pDNA")])
        counts['cell_line'] = counts.variable_0.str.split(".").map(lambda x: x[0])
        counts['replicate'] = counts.variable_1.str[3:].map({"A": 1, "B": 2, "C": 3})
        counts['timepoint'] = 1
        counts.drop(columns=["variable_0", "variable_1"], inplace=True)

        counts_final = counts.drop(columns=[("pDNA", "pDNA")])
        counts_pdna = counts.drop(columns=["value"])
        counts_pdna['timepoint'] = 0
        counts_pdna.rename(columns={("pDNA", "pDNA"): "read_count"}, inplace=True)
        counts_final.rename(columns={"value": "read_count"}, inplace=True)
        
        counts = pd.concat([counts_pdna, counts_final])
        counts = counts.merge(grnas, left_on="sgRNA", right_on="sgRNA Sequence", how="left")
        counts.drop(columns=["sgRNA Sequence", "Annotated Gene ID"], inplace=True)
        counts.rename(columns={"Annotated Gene Symbol": "gene"}, inplace=True)

        counts["crispr_system"] = "CRISPRi"
        counts["sgRNA"] = counts['sgRNA'].str.upper()
        return counts

    seta_counts = process_doench2018_sheet(sheets['SetA raw reads'], sheets['SetA sgRNA annotations'])
    seta_counts["publication"] = "Doench2018-SetA"
    setb_counts = process_doench2018_sheet(sheets['SetB raw reads'], sheets['SetB sgRNA annotations'])
    setb_counts["publication"] = "Doench2018-SetB"

    counts = pd.concat([seta_counts, setb_counts])
    counts["sgRNA_ID"] = counts.sgRNA + "_" + counts.gene
    return counts

def preprocess_weissman2016():
    crispra_counts = pd.read_csv("data/inputs/weissman2016_crispra_counts.csv")
    crispra_grnas = pd.read_csv("data/inputs/weissman2016_crispra_grnas.csv")
    crispri_counts = pd.read_csv("data/inputs/weissman2016_crispri_counts.csv")
    crispri_grnas = pd.read_csv("data/inputs/weissman2016_crispri_grnas.csv")

    def preprocess(counts, grnas):
        counts = pd.melt(counts, id_vars="sgId")
        counts['timepoint'] = counts.variable.str.split(",").map(lambda x: x[0])
        counts['timepoint'] = counts.timepoint.map({"Endpoint": 1, "T0": 0})
        counts['replicate'] = counts.variable.str.split(",").map(lambda x: x[1].strip()[3:]).astype(int)
        counts = counts.merge(grnas, how='left', left_on='sgId', right_on='sgID')

        counts.drop(columns=['variable', 'sgID'], inplace=True)

        counts.rename(columns={
            "value": "read_count",
            "sgId": "sgRNA_ID",
            "protospacer sequence": "sgRNA"
        }, inplace=True)

        counts["publication"] = "Weissman2016"
        counts["sgRNA"] = counts["sgRNA"].str.upper()
        return counts

    crispra_counts = preprocess(crispra_counts, crispra_grnas)
    crispra_counts['crispr_system'] = 'CRISPRa'

    crispri_counts = preprocess(crispri_counts, crispri_grnas)
    crispri_counts['crispr_system'] = 'CRISPRi'

    counts = pd.concat([crispri_counts, crispra_counts])
    counts['cell_line'] = 'K562'
    return counts

def preprocess_lim2017():
    crispri_counts = pd.read_csv("data/inputs/lim2017_crispri_counts.csv", header=[0,1,2])\
                       .rename(columns={"Unnamed: 0_level_0": "sgRNA_ID",
                                        "Unnamed: 0_level_1": "",
                                        "Unnamed: 0_level_2": ""})
    crispri_counts = pd.melt(crispri_counts, id_vars="sgRNA_ID")
    crispri_counts.rename(columns={
        "variable_0": "cell_line",
        "variable_1": "timepoint",
        "variable_2": "replicate",
        "value": "read_count"
    }, inplace=True)

    crispri_counts["timepoint"] = crispri_counts.timepoint.map({"T0": 0, "T12": 1})
    crispri_counts["replicate"] = crispri_counts.replicate.map({"Rep1": 1, "Rep2": 2})
    crispri_counts["publication"] = "Lim2017"
    crispri_counts["crispr_system"] = "CRISPRi"
    crispri_counts.dropna(inplace=True)

    crispri_grnas_1 = pd.read_csv("data/inputs/lim2017_crispri_grnas_1.csv")
    crispri_grnas_2 = pd.read_csv("data/inputs/lim2017_crispri_grnas_2.csv")

    crispri_grnas = pd.concat([crispri_grnas_1, crispri_grnas_2])

    crispri_counts = crispri_counts.merge(
        crispri_grnas[["sgRNA ID", "sgRNA sequence", "Gene ID"]],
        left_on="sgRNA_ID", right_on="sgRNA ID", how="left"
    )

    crispri_counts.rename(columns={
        "Gene ID": "gene",
        "sgRNA sequence": "sgRNA"
    }, inplace=True)

    crispri_counts.drop(columns=["sgRNA ID"], inplace=True)
    crispri_counts['sgRNA'] = crispri_counts['sgRNA'].str.upper()
    crispri_counts["sgRNA_ID"] = crispri_counts["sgRNA_ID"] + "_" + crispri_counts["sgRNA"]
    return crispri_counts.drop_duplicates()

def preprocess_weissman2021():
    crispri_hek293t = pd.read_csv("data/inputs/weissman2021_crispri_hek293t_counts.csv")
    crisproff_hek293t = pd.read_csv("data/inputs/weissman2021_crisproff_hek293t_counts.csv")
    crisproff_hek293t_mutant = pd.read_csv("data/inputs/weissman2021_crisproff_hek293t_mutant_counts.csv")
    crisproff_grnas = pd.read_csv("data/inputs/weissman2021_crisproff_grnas.csv")
    crisproff_grnas.rename(columns={"id": "sgRNA_ID", "sequence": "sgRNA"}, inplace=True)
    crisproff_grnas["sgRNA"] = crisproff_grnas["sgRNA"].str.upper()

    count_dfs = [("HEK293TMUT", "CRISPRoff", crisproff_hek293t_mutant),
           ("HEK293T", "CRISPRoff", crisproff_hek293t),
           ("HEK293T", "CRISPRi", crispri_hek293t)]
    dfs = []
    for cell_line, crispr_system, count_df in count_dfs:
        count_df.drop(columns=["Unnamed: 5", "Rep1", "Rep2", "ave_Rep1_Rep2"], inplace=True)
        count_df.rename(columns={"id": "sgRNA_ID"}, inplace=True)
        count_df = count_df.melt(id_vars="sgRNA_ID")
        count_df['timepoint'] = count_df.variable.str.split(",").map(lambda x: x[0])
        count_df['timepoint'] = count_df.timepoint.map({"T0": 0, "Tfinal": 1})
        count_df['replicate'] = count_df.variable.str.split(",").map(lambda x: x[1].strip())
        count_df['replicate'] = count_df.replicate.map({"Rep1": 1, "Rep2": 2})
        count_df = count_df.drop(columns="variable").rename(columns={"value": "read_count"})
        count_df['cell_line'] = cell_line
        count_df['crispr_system'] = crispr_system
        count_df['publication'] = 'Weissman2021'
        count_df = count_df.merge(crisproff_grnas, on="sgRNA_ID", how="left")
        dfs.append(count_df)
    return pd.concat(dfs)

if __name__ == "__main__":
    args = parse_arguments()

    lim2017 = preprocess_lim2017()
    weissman2016 = preprocess_weissman2016()
    doench2018 = preprocess_doench2018()
    weissman2021 = preprocess_weissman2021()

    all_data = pd.concat([lim2017, weissman2016, doench2018, weissman2021])
    all_data.to_csv("data/preprocessed_screen_counts.csv")

    control_names = ['CONTROL', 'CTRL', 'negative_ctrl', 'negative_control']
    conditions = all_data.groupby(["publication", "crispr_system", "cell_line", "replicate"])
    for ((publication, crispr_system, cell_line, replicate), count_df) in conditions: 
        initial_timepoint = count_df[count_df['timepoint'] == 0]
        final_timepoint = count_df[count_df['timepoint'] == 1]

        mageck_df = initial_timepoint.merge(final_timepoint, on='sgRNA_ID', suffixes=('_initial', '_final'))
        mageck_df = mageck_df[['sgRNA_ID', 'gene_initial', 'read_count_initial', 'read_count_final']]
        mageck_df.rename(columns={
            "sgRNA_ID": "sgRNA",
            "gene_initial": "gene",
            "read_count_initial": "initial",
            "read_count_final": "final"
        }, inplace=True)

        mageck_ctrls_df = mageck_df[mageck_df['gene'].isin(control_names)][["sgRNA"]]

        mageck_df.to_csv(f"data/mageck_counts/{publication}_{crispr_system}_{cell_line}_{replicate}_counts.txt", index=False, sep='\t')
        mageck_ctrls_df.to_csv(f"data/mageck_counts/{publication}_{crispr_system}_{cell_line}_{replicate}_controls.txt", index=False, sep='\t')

    # all_data = all_data[all_data["publication"] == "Weissman2021"]
    guidescan_data = all_data.rename(columns={
        "sgRNA_ID": "id",
        "sgRNA": "sequence"
    })[["id", "sequence"]]

    guidescan_data["id"] = guidescan_data["id"].str.split(",").map(lambda x: x[0])
    guidescan_data["pam"] = "NGG"
    guidescan_data["chromosome"] = "1"
    guidescan_data["position"] = "1"
    guidescan_data["sense"] = "+"

    guidescan_data.to_csv("data/preprocessed_screen_kmers.csv", index=False)
