import pandas as pd
import argparse
import random

from funcy import ilen

A549_COLUMNS = [
    'A549_Vehicle_Lung_1_A01',
    'A549_Vehicle_Lung_2_A05'
]

THRESHOLD         = 5000

NUM_NON_ESSENTIAL = 135
NUM_ESSENTIAL     = 100
NUM_CONTROLS      = 100

def build_control_guides(all_libraries):
    seeds = [random.randint(0, 1000), random.randint(0, 1000)]
    
    non_targeting_controls = all_libraries[all_libraries['Type'] == 'non_targeting_control'].sample(NUM_CONTROLS, random_state=seeds[0])
    safe_targeting_controls = all_libraries[all_libraries['Type'] == 'safe_targeting_control'].sample(NUM_CONTROLS, random_state=seeds[1])
    return pd.concat([non_targeting_controls, safe_targeting_controls])

def build_non_essential_guides(all_libraries, non_essential_genes):
    non_essential_guides = all_libraries[all_libraries['Gene'].isin(non_essential_genes)].copy()

    bassik_ne_guides = non_essential_guides[non_essential_guides['Library'] == 'Bassik2017'].copy()
    bassik_ne_guides['Num Off-Targets'] = bassik_ne_guides[['0 Off-Targets', '1 Off-Targets', '2 Off-Targets', '3 Off-Targets']].sum(axis=1)
    bassik_ne_genes = bassik_ne_guides.groupby('Gene')['Num Off-Targets'].apply(lambda df: df.nlargest(2))
    bassik_ne_genes = set(bassik_ne_genes.groupby('Gene').mean().nlargest(NUM_NON_ESSENTIAL).index)
    bassik_ne_guides = bassik_ne_guides[bassik_ne_guides['Gene'].isin(bassik_ne_genes)]

    non_essential_guides = non_essential_guides.groupby(['Library', 'Gene']) \
                                               .Specificity \
                                               .apply(lambda df: df.nsmallest(2)) \
                                               .groupby(['Library', 'Gene']) \
                                               .mean() \
                                               .groupby('Library') \
                                               .nsmallest(NUM_NON_ESSENTIAL).index.tolist()

    non_essential_guides = [(x, z) for (x, _, z) in non_essential_guides]
    non_essential_guides = pd.DataFrame(non_essential_guides, columns=['Library', 'Gene'])
    non_essential_guides = non_essential_guides[non_essential_guides['Library'] != 'Guidescan']
    non_essential_guides = non_essential_guides[non_essential_guides['Library'] != 'Bassik2017']
    non_essential_guides = non_essential_guides.merge(all_libraries, how='inner')

    non_essential_guides = pd.concat([non_essential_guides, bassik_ne_guides])

    gs_ne_genes = set(non_essential_guides['Gene'])
    gs_ne_guides = all_libraries[(all_libraries['Library'] == 'Guidescan') & (all_libraries['Gene'].isin(gs_ne_genes))]

    non_essential_guides = pd.concat([non_essential_guides, gs_ne_guides]).drop(columns=['Num Off-Targets'])
    non_essential_guides['Type'] = 'non_essential_gene_targeting'

    return non_essential_guides

def build_essential_guides(all_libraries, essential_genes):
    essential_guides = all_libraries[all_libraries['Gene'].isin(essential_genes)].copy()
    essential_guides['Type'] = 'essential_gene_targeting'
    return essential_guides

def load_all_df(args):
    library_renaming_map = {
        'Doench': 'Root2016',
        'Hart': 'Moffat2015',
        'Bassik': 'Bassik2017',
        'Zuber2020': 'Elling2020',
        'Sabatini': 'Sabatini2015',
    }

    guidescan_library = pd.read_csv(args.guidescan_library_csv)
    guidescan_library['Library'] = 'Guidescan'
    guidescan_library['0 Off-Targets'] = 1

    soa_df = pd.read_csv(args.state_of_the_art_csv).drop_duplicates(
        subset=['sgRNA', 'library']
    )

    soa_df.rename(columns={
        'pos': 'Pos',
        'gene': 'Gene',
        'specificity': 'Specificity',
        'cutting_efficiency': 'Cutting Efficiency',
        'chr': 'Chr',
        'strand': 'Strand',
        'library': 'Library',
        'identifier': 'Identifier',
    }, inplace=True)

    soa_df.drop(columns=[
        'score', 'distance', 'absolute_pos'
    ], inplace=True)

    soa_df['Library'] = soa_df['Library'].apply(lambda lib: library_renaming_map[lib])

    return pd.concat([guidescan_library, soa_df])

def get_non_essential_genes(args, expressed_genes, common_genes):
    non_essential_df    = pd.read_csv(args.non_essential_genes_csv)
    non_essential_genes = (set(non_essential_df['Nonessential Genes (NE)']) - expressed_genes) & common_genes
    return set(non_essential_genes)

def get_essential_genes(args, expressed_genes):
    essential_genes = set(map(
        lambda x: x.split()[0],
        open(args.essential_genes_txt).readlines()
    ))

    essential_genes = essential_genes & expressed_genes
    essential_genes = sorted(essential_genes) # to ensure deterministic operation
    essential_genes = random.sample(essential_genes, NUM_ESSENTIAL)
    return set(essential_genes)

def get_expressed_genes(common_genes, args):
    gene_expression_df = pd.read_csv(args.gene_expression_csv)

    gene_expression_df[A549_COLUMNS] = gene_expression_df[A549_COLUMNS] / gene_expression_df[A549_COLUMNS].max()
    gene_expression_df['A549_average'] = gene_expression_df[A549_COLUMNS].mean(axis=1)

    gene_expression_df = gene_expression_df[gene_expression_df['Gene_name'].isin(common_genes)]
    expressed_genes = gene_expression_df.sort_values(by=['A549_average'], ascending=False)\
                                        .iloc[:THRESHOLD]['Gene_name']
    return set(expressed_genes)

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Designs the essentialome screen.'
    )

    parser.add_argument(
        'guidescan_library_csv',
        help='Guidescan library CSV file.'
    )

    parser.add_argument(
        'state_of_the_art_csv',
        help=('CSV file containing Guidescan processed state-of-the-art '
              'sgRNAs.')
    )
    
    parser.add_argument(
        'gene_expression_csv',
        help='Gene expression CSV'
    )

    parser.add_argument(
        'essential_genes_txt',
        help=('TXT file containing essential genes from TableS2 '
              'of https://doi.org/10.1534/g3.117.041277')
    )

    parser.add_argument(
        'non_essential_genes_csv',
        help=('CSV file containing non-essentials from Supplementary'
              'Dataset S1 of https://doi.org/10.15252/msb.20145216')
    )

    parser.add_argument(
        '-o', '--outfile',
        default='out.csv',
        help='Name of output CSV.'
    )

    return parser.parse_args()

def main():
    args = parse_arguments()

    random.seed(74) # easter egg: email me the closest prime number to this seed

    all_libraries       = load_all_df(args)
    common_genes        = set.intersection(*all_libraries.groupby('Library').Gene.agg(lambda x: set(x)))

    expressed_genes     = get_expressed_genes(common_genes, args)
    essential_genes     = get_essential_genes(args, expressed_genes)
    non_essential_genes = get_non_essential_genes(args, expressed_genes, common_genes)

    essential_guides = build_essential_guides(
        all_libraries, essential_genes
    )

    non_essential_guides = build_non_essential_guides(
        all_libraries, non_essential_genes
    )

    controls = build_control_guides(all_libraries)

    screen_df = pd.concat([essential_guides, non_essential_guides, controls])
    unnamed_columns = list(filter(lambda s: s.startswith('Unnamed'), screen_df.columns))
    screen_df.drop(columns=unnamed_columns, inplace=True)
    screen_df['G/C Content'] = screen_df['sgRNA'].apply(
        lambda sgrna: ilen(filter(lambda n: n in ['G', 'C'], sgrna)) / len(sgrna)
    )

    screen_df = screen_df.rename(columns={
        '0 Off-Targets': 'Distance 0 Matches',
        '1 Off-Targets': 'Distance 1 Matches',
        '2 Off-Targets': 'Distance 2 Matches',
        '3 Off-Targets': 'Distance 3 Matches',
    })

    screen_df.to_csv(args.outfile)

if __name__ == '__main__':
    main()
