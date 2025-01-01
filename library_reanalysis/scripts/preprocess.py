import pandas as pd
import numpy as np

"""
| Filename | Species | Name | Lab | DOI |
|---------|---------|------|-----|-----|
| `teichman_mouse_v2_grnas.csv` | Mouse | Teichmann CRISPRko | Teichmann Lab | 10.1016/j.cell.2018.11.044 |
| `m1_grna_sequence.txt` | Mouse | Liu CRISPRko M1 | Liu Lab | 10.1158/2159-8290.CD-20-0812 | 
| `m2_grna_sequence.txt` | Mouse | Liu CRISPRko M2 | Liu Lab | 10.1158/2159-8290.CD-20-0812 |
| `guidescan2_mouse.csv` | Mouse | GuideScan Mouse | Pritykin Lab | 10.1101/2022.05.02.490368 |
| `tkov3_guide_sequence.csv` | Human | TKOv3 | Moffat Lab | 10.1016/j.celrep.2019.02.041 |
| `minlibcas9_library.csv` | Human | MinLibCas9 | Garnett Lab | 10.1186/s13059-021-02268-4 | 
| `dgrna_metadata.csv` | Human | minimized double gRNA CRISPRko | Parts Lab | 10.1101/859652 | 
| `guidescan2_human.csv` | Mouse | GuideScan Mouse | Pritykin Lab | 10.1101/2022.05.02.490368 |
| `gattinara_human.csv` | Human | Gattinara | Doench Lab | 10.1038/s41467-020-14620-6 |
| `gouda_human.csv` | Mouse | Gouda | Doench Lab | 10.1038/s41467-020-14620-6 |
"""

""" Preprocess the data from the different  input 
files into a single dataframe with columns:
    - id: unique identifier for each gRNA
    - library: name of the library 
    - species: species of the library
    - gene: gene targeted by the gRNA
    - sequence: gRNA sequence
"""
def preprocess_guidescan2(species):
    guidescan2_df = pd.read_csv(f"inputs/{species}/guidescan2_{species}.csv") 
    guidescan2_df = guidescan2_df[guidescan2_df['type'] == 'gene_targeting']
    guidescan2_df = guidescan2_df.rename(columns={"gRNA": "sequence"})
    guidescan2_df['library'] = 'guidescan2'
    guidescan2_df['species'] = species
    guidescan2_df = guidescan2_df[['id', 'library', 'species', 'gene', 'sequence']]
    return guidescan2_df

def preprocess_gouda():
    gouda_df = pd.read_csv("inputs/mouse/gouda_sequences.txt", sep='\t')
    gouda_df = gouda_df.rename(columns={"Barcode Sequence": "sequence", "Annotated Gene Symbol": "gene"})
    gouda_df['library'] = 'gouda'
    gouda_df['species'] = 'mouse'
    gouda_df['id'] = gouda_df['gene'] + '_' + (gouda_df.groupby('gene').cumcount() + 1).astype(str)
    gouda_df = gouda_df[['id', 'library', 'species', 'gene', 'sequence']]
    return gouda_df

def preprocess_gattinara():
    gattinara_df = pd.read_csv("inputs/human/gattinara_sequences.txt", sep='\t')
    gattinara_df = gattinara_df.rename(columns={"Barcode Sequence": "sequence", "Annotated Gene Symbol": "gene"})
    gattinara_df['library'] = 'gattinara'
    gattinara_df['species'] = 'human'
    gattinara_df['id'] = gattinara_df['gene'] + '_' + (gattinara_df.groupby('gene').cumcount() + 1).astype(str)
    gattinara_df = gattinara_df[['id', 'library', 'species', 'gene', 'sequence']]
    return gattinara_df

def preprocess_minlabcas9():
    minlabcas9_df = pd.read_csv("inputs/human/minlibcas9_library.csv")
    minlabcas9_df = minlabcas9_df[minlabcas9_df['Library'] == 'KosukeYusa']
    minlabcas9_df = minlabcas9_df[~minlabcas9_df['sgRNA_ID'].str.contains('CTRL')]
    minlabcas9_df = minlabcas9_df.rename(columns={"WGE_Sequence": "sequence", "sgRNA_ID": "id", "Approved_Symbol": "gene"})
    minlabcas9_df['library'] = 'minlabcas9'
    minlabcas9_df['species'] = 'human'
    minlabcas9_df = minlabcas9_df[['id', 'library', 'species', 'gene', 'sequence']]
    minlabcas9_df['sequence'] = minlabcas9_df['sequence'].str[0:20]
    return minlabcas9_df

def preprocess_minimized_dgrna():
    dgrna_df = pd.read_csv("inputs/human/dgrna_metadata.csv")
    dgrna_df = dgrna_df.rename(columns={"gRNA_sequence": "sequence", "Gene_name": "gene"})
    dgrna_df = dgrna_df[~dgrna_df['gene'].isna()]
    dgrna_df = dgrna_df[~dgrna_df['gene'].str.contains('control')]
    dgrna_df['library'] = 'minimized_dgrna'
    dgrna_df['species'] = 'human'
    dgrna_df['id'] = dgrna_df['gene'] + '_' + dgrna_df['sequence']
    dgrna_df = dgrna_df[['id', 'library', 'species', 'gene', 'sequence']]
    return dgrna_df

def preprocess_tkov3():
    tkov3_df = pd.read_csv("inputs/human/tkov3_guide_sequence.csv")
    tkov3_df = tkov3_df[~tkov3_df['TARGET EXON'].str.contains('CTRL')]
    tkov3_df = tkov3_df.rename(columns={"SEQUENCE": "sequence", "GUIDE_ID": "id", "GENE": "gene"})
    tkov3_df['library'] = 'tkov3'
    tkov3_df['species'] = 'human'
    tkov3_df = tkov3_df[['id', 'library', 'species', 'gene', 'sequence']]
    return tkov3_df

def preprocess_liu():
    liu_df_m1 = pd.read_csv("inputs/mouse/m1_grna_sequence.txt", sep='\t')
    liu_df_m2 = pd.read_csv("inputs/mouse/m2_grna_sequence.txt", sep='\t')
    liu_df = pd.concat([liu_df_m1, liu_df_m2])
    liu_df = liu_df.rename(columns={"sgrna": "id", "seq": "sequence", "gene": "gene"})
    liu_df['library'] = 'liu'
    liu_df['species'] = 'mouse'
    liu_df['sequence'] = liu_df['sequence'] + 'G'
    liu_df = liu_df[~liu_df['gene'].isna()]
    liu_df = liu_df[~liu_df['gene'].str.contains('CTRL')]
    liu_df = liu_df[['id', 'library', 'species', 'gene', 'sequence']]
    return liu_df

def preprocess_teichmann():
    teichmann_df = pd.read_csv("inputs/mouse/teichman_mouse_v2_grnas.csv")
    teichmann_df = teichmann_df.rename(columns={"guide_sequence": "sequence", "gRNA_ID": "id", "gene": "gene"})
    teichmann_df['library'] = 'teichmann'
    teichmann_df['species'] = 'mouse'
    teichmann_df['sequence'] = teichmann_df['sequence'] + 'G'
    teichmann_df = teichmann_df[['id', 'library', 'species', 'gene', 'sequence']]
    return teichmann_df

if __name__ == "__main__":
    dfs = []
    for species in ['human', 'mouse']:
        df = preprocess_guidescan2(species)
        dfs.append(df)

    df = preprocess_minlabcas9()
    dfs.append(df)

    df = preprocess_minimized_dgrna()
    dfs.append(df)

    df = preprocess_tkov3()
    dfs.append(df)

    df = preprocess_liu()
    dfs.append(df)

    df = preprocess_teichmann()
    dfs.append(df)

    df = preprocess_gouda()
    dfs.append(df)

    df = preprocess_gattinara()
    dfs.append(df)

    all_libs_df = pd.concat(dfs)

    # compute the summary statistics of the number of gRNAs per gene in each library:
    print(all_libs_df.groupby(['library', 'species', 'gene']).size().reset_index(name='counts').groupby(['library', 'species'])['counts'].describe())

    all_libs_df.to_csv("outputs/all_libraries.csv", index=False)

    all_libs_df['pam'] = 'NGG'
    all_libs_df['chromosome'] = '1'
    all_libs_df['position'] = 1
    all_libs_df['sense'] = '+'
    all_libs_df = all_libs_df[['id','sequence','pam','chromosome','position','sense','species']]
    all_libs_df[all_libs_df['species'] == 'mouse'].drop(columns=['species']).to_csv("outputs/all_libraries_mouse_kmers.csv", index=False)
    all_libs_df[all_libs_df['species'] == 'human'].drop(columns=['species']).to_csv("outputs/all_libraries_human_kmers.csv", index=False)
