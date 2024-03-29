import argparse
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('new_scored_libraries')
    parser.add_argument('old_scored_libraries')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    scored_library_df = pd.read_csv(args.new_scored_libraries)
    old_scored_library_df = pd.read_csv(args.old_scored_libraries)
    old_scored_library_df = old_scored_library_df.rename(columns={'gRNA': 'sequence'})

    scored_library_df = scored_library_df[scored_library_df['species'] == 'mouse']
    scored_library_df = scored_library_df[['library', 'id', 'gene', 'sequence', 'specificity']]
    scored_library_df['old'] = False
    print(scored_library_df.sample(100))
    old_scored_library_df = old_scored_library_df[['library', 'id', 'gene', 'sequence', 'specificity']]
    old_scored_library_df['old'] = True
    old_scored_library_df = old_scored_library_df[old_scored_library_df['library'] != 'Elling2020-Screen']
    old_scored_library_df = old_scored_library_df[old_scored_library_df['library'] != 'Weissman2016-CRISPRi']
    old_scored_library_df = old_scored_library_df[old_scored_library_df['library'] != 'Weissman2016-CRISPRa']

    concat_df = pd.concat([scored_library_df, old_scored_library_df])
    concat_df = concat_df[(concat_df['sequence'].str.len() == 20) | (concat_df['library'] != 'Bassik2017')]
    rename_dict = {
        'guidescan2': 'GuideScan2', 
        'minlabcas9': 'Garnett2021', 
        'minimized_dgrna': 'Parts2019', 
        'tkov3': 'Moffat2019', 
        'liu': 'Liu2021', 
        'teichmann': 'Teichman2019',
        'Elling2020-Library': 'Elling2020',
    }
    concat_df['library'] = concat_df['library'].map(lambda x: rename_dict[x] if x in rename_dict.keys() else x)

    fig, ax = plt.subplots(figsize=(5, 4))
    sns.violinplot(x='library', y='specificity', data=concat_df, ax=ax)
    ax.set_xlabel('Library')
    ax.set_ylabel('Specificity')

    # rotate x-axis labels
    for tick in ax.get_xticklabels():
        tick.set_rotation(60)
    fig.tight_layout()

    fig.savefig('figures/mouse_all_libraries_specificity.pdf')

    fig, ax = plt.subplots(figsize=(5, 4), ncols=1)
    sns.violinplot(x='library', y='specificity', data=concat_df[concat_df['old'] == False], ax=ax)
    ax.set_xlabel('Library')
    ax.set_ylabel('Specificity')

    for tick in ax.get_xticklabels():
        tick.set_rotation(60)

    fig.tight_layout()
    fig.savefig('figures/mouse_new_libraries_specificity.pdf')

    plt.show()