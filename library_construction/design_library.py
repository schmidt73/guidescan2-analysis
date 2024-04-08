import argparse
import numpy as np
import pandas as pd
import sys

from funcy import ilen

TIERS = [(1, 0.7), (0.7, 0.45), (0.45, 0.30), (0.30, 0.15), (0.15, 0)]
VERBOSE = False

def log(*args, **kwargs):
    if VERBOSE:
        print(*args, file=sys.stderr, **kwargs)

def rank(xs, x):
    for i in range(len(xs)):
        h, l = xs[i]
        if x <= h and x > l:
            return len(xs) - i
    return 0

def contains_monopolymer(seq, k=4):
    count, last = 1, seq[0]
    for i in range(1, len(seq)):
        if count >= k:
            return True

        if seq[i] == last:
            count += 1
            continue
        
        count = 0
        last = seq[i]
        
    if count >= k:
        return True 
    
    return False

def rank_guides(guide):
    gc_cont = guide['G/C Content']
    
    if guide['type'] == 'exon':
        return -1
    
    if gc_cont < 0.20 or gc_cont > 0.80:
        return -2
    
    sgrna = guide['sgRNA']
    if contains_monopolymer(sgrna, 4):
        return -3
    
    return min(guide['min_specificity'], 1.25 * guide['cutting_efficiency'])

def min_revision(guidescan_df):
    ranker = lambda g: rank_guides(g)

    guidescan_df = guidescan_df.copy()
    guidescan_df['rank'] = guidescan_df.apply(ranker, axis=1)

    log(
        f'Number of G/C content filtered sgRNAs: '
        f'{len(guidescan_df[~guidescan_df["rank"].isin([-1, -2])])}'
    )


    log(
        f'Number of monopolymer filtered sgRNAs: '
        f'{len(guidescan_df[~guidescan_df["rank"].isin([-1, -2, -3])])}'
    )

    g = guidescan_df.groupby('gene').apply(lambda df: df.nlargest(6, ['rank']))
    
    return g.reset_index(drop=True)

def build_safe_targeting_controls(args):
    safe_targeting_controls = pd.read_csv(args.safe_targeting_csv)
    safe_targeting_controls.loc[safe_targeting_controls['antisense'], 'Strand'] = '-'
    safe_targeting_controls.loc[~safe_targeting_controls['antisense'], 'Strand'] = '-'
    safe_targeting_controls['sgRNA'] = safe_targeting_controls['sgRNA'].str[:-3]
    safe_targeting_controls['PAM'] = 'NGG'

    safe_targeting_controls = safe_targeting_controls.drop(
        columns=['end', 'antisense']
    ).rename(
        columns={
            'specificity': 'Specificity',
            'cutting_efficiency': 'Cutting Efficiency',
            'safe_targeting_region': 'Safe Targeting Region',
            'start': 'Pos',
            'chr': 'Chr'
        }
    )

    safe_targeting_controls['Type'] = 'safe_targeting_control'
    safe_targeting_controls = safe_targeting_controls.sort_values(by='Specificity', ascending=False)
    safe_targeting_controls = safe_targeting_controls.iloc[:5000]

    return safe_targeting_controls

def build_non_targeting_controls(args):
    non_targeting_controls = pd.read_table(args.non_targeting_txt, header=None)
    non_targeting_controls = non_targeting_controls.rename(columns={0: 'sgRNA'})

    non_targeting_controls['Type'] = 'non_targeting_control'
    non_targeting_controls['0 Off-Targets'] = 0
    non_targeting_controls['1 Off-Targets'] = 0
    non_targeting_controls['2 Off-Targets'] = 0
    non_targeting_controls['3 Off-Targets'] = 0

    return non_targeting_controls

def build_gene_targeting_guides(args):
    guidescan_df = pd.read_csv(args.guidescan_csv)

    log(f'Number of exon filtered sgRNAs: {len(guidescan_df[guidescan_df["type"] == "exon"])}')
    log(f'Number of CDS filtered sgRNAs: {len(guidescan_df[guidescan_df["type"] == "CDS"])}')

    guidescan_df['gene'] = guidescan_df['gene'].str[5:]
    guidescan_df['identifier'] = guidescan_df['identifier'].str[5:]
    guidescan_df['sgRNA'] = guidescan_df['sgRNA'].str[:-3]
    guidescan_df['PAM'] = 'NGG'
    guidescan_df = guidescan_df[~guidescan_df['5pG Specificity'].isna()]
    guidescan_df['min_specificity'] = guidescan_df[['specificity', '5pG Specificity']].min(axis=1)

    guidescan_df = guidescan_df[guidescan_df['cutting_efficiency'] > 0.25]

    log(
        f'Number of cutting efficiency (> 0.25) filtered sgRNAs: '
        f'{len(guidescan_df[guidescan_df["type"] == "CDS"])}'
    )

    guidescan_df = guidescan_df[guidescan_df['specificity'] > 0.20]

    log(
        f'Number of specificity (> 0.20) filtered sgRNAs: '
        f'{len(guidescan_df[guidescan_df["type"] == "CDS"])}'
    )

    guidescan_df['G/C Content'] = guidescan_df['sgRNA'].apply(
        lambda sgrna: ilen(filter(lambda n: n in ['G', 'C'], sgrna)) / len(sgrna)
    )

    gene_targeting = min_revision(guidescan_df)

    gene_targeting.rename(columns={
        'gene': 'Gene',
        'identifier': 'Identifier',
        'cutting_efficiency': 'Cutting Efficiency',
        'specificity': 'Specificity',
        'type': 'Cuts In',
        'strand': 'Strand',
        'id': 'Cutting Region ID',
        'chr': 'Chr',
        'pos': 'Pos',
    }, inplace=True)

    gene_targeting['Type'] = 'gene_targeting'

    gene_targeting.drop(columns=[
        'min_specificity', 'start'
    ], inplace=True)

    return gene_targeting

def parse_arguments():
    parser = argparse.ArgumentParser(
        description=
        'Constructs the final Guidescan library from the Guidescan '
        'processed CDS/exon cutting sgRNAs.'
    )

    parser.add_argument(
        'guidescan_csv',
        help='Guidescan processed CSV'
    )

    parser.add_argument(
        'safe_targeting_csv',
        help='Guidescan processed safe-targeting controls CSV'
    )

    parser.add_argument(
        'non_targeting_txt',
        help='Guidescan processed non-targeting controls TXT'
    )

    parser.add_argument(
        '-o', '--outfile',
        default='out.csv',
        help='Name of output CSV.'
    )

    parser.add_argument(
        '-v', '--verbose', action='store_true',
        help='Verbose mode.'
    )
    
    return parser.parse_args()

def main():
    args = parse_arguments()

    global VERBOSE
    VERBOSE = args.verbose

    non_targeting_controls = build_non_targeting_controls(args)
    safe_targeting_controls = build_safe_targeting_controls(args)
    gene_targeting_guides = build_gene_targeting_guides(args)

    full_library = pd.concat([
        gene_targeting_guides,
        non_targeting_controls,
        safe_targeting_controls,
    ])

    unnamed_columns = list(filter(lambda s: s.startswith('Unnamed'), full_library.columns))
    full_library.drop(columns=unnamed_columns, inplace=True)
    full_library.drop(columns=['rank'], inplace=True)
    full_library['G/C Content'] = full_library['sgRNA'].apply(
        lambda sgrna: ilen(filter(lambda n: n in ['G', 'C'], sgrna)) / len(sgrna)
    )

    full_library.to_csv(args.outfile)

if __name__ == '__main__':
    main()
