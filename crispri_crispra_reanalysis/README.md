# Guidescan2 CRISPRi/a Analysis

This repository contains the scripts, data, and instructions 
required to perform a analysis of existing CRISPRi/a screens 
using Guidescan2. 

## Input Files

The input files are located in the `inputs` directory. They are pulled from 
several publications which perform genome-wide perturbation screens using CRISPRi/a/off. 
The files are located under `data/inputs` and are described in the table below.

| File Name | Source (Publication) | Description |
| --------- | -------------------- | ----------- |
| `doench2018_crispri_counts.xlsx` | *Optimized libraries for CRISPR-Cas9 genetic screens with multiple modalities* | Measured gRNA abundances and gRNA sequences |
| `lim2017_crispri_counts.csv` | *CRISPRi-based genome-scale identification of functional long noncoding RNA loci in human cells* | Measured gRNA abundances |
| `lim2017_crispri_grnas_1.csv` | *CRISPRi-based genome-scale identification of functional long noncoding RNA loci in human cells* | gRNA sequences (Set 1) |
| `lim2017_crispri_grnas_2.csv` | *CRISPRi-based genome-scale identification of functional long noncoding RNA loci in human cells* | gRNA sequences (Set 2) |
| `weissman2016_crispra_counts.csv` | *Compact and highly active next-generation libraries for CRISPR-mediated gene repression and activation* | Measured gRNA abundances (CRISPRi) |
| `weissman2016_crispra_grnas.csv` | *Compact and highly active next-generation libraries for CRISPR-mediated gene repression and activation* | gRNA sequences (CRISPRi) |
| `weissman2016_crispri_counts.csv` | *Compact and highly active next-generation libraries for CRISPR-mediated gene repression and activation* | Measured gRNA abundances (CRISPRa) |
| `weissman2016_crispri_grnas.csv` | *Compact and highly active next-generation libraries for CRISPR-mediated gene repression and activation* | gRNA sequences (CRISPRa) |
| `weissman2021_crispri_hek293t_counts.csv` | *Genome-wide programmable transcriptional memory by CRISPR-based epigenome editing* | Measured gRNA abundances (CRISPRi) |
| `weissman2021_crisproff_grnas.csv` | *Genome-wide programmable transcriptional memory by CRISPR-based epigenome editing* | gRNA sequences |
| `weissman2021_crisproff_hek293t_counts.csv` | *Genome-wide programmable transcriptional memory by CRISPR-based epigenome editing* | Measured gRNA abundances |
| `weissman2021_crisproff_hek293t_mutant_counts.csv` | *Genome-wide programmable transcriptional memory by CRISPR-based epigenome editing* | Measured gRNA abundances (mutant) |

In short, for each dataset, there are two files: one containing the gRNA sequences 
and one containing the gRNA abundances. The gRNA sequence files contain the true gRNA 
sequences used in the screen. The gRNA abundance files contain the measured gRNA 
abundances for each gRNA in the screen in both the initial and final timepoints
across several cell lines.

For the Marson et al. 2022 publication titled *CRISPR activation and interference screens decode stimulation responses in primary human T cells*
the raw abundance data was not available. Instead, the authors provided the 
MAGeCK output files for the CRISPRi and CRISPRa screens. These files are located
in the `data/inputs/marson2022` directory.

## Instructions

First, we need to preprocess the input files into a format that Guidescan2 and MAGeCK can 
understand. This is done using the `preprocess_inputs.py` script. From the root
directory of this repository, run the following command:

```
$ python preprocess_inputs.py
```

This will create a new directory called `data/mageck_counts` which contains the
gRNA abundnace files in a format that MAGeCK can understand, split by 
publication, cell line, perturbation type, timepoint, and replicate. It will
also create a master file called `data/preprocessed_screen_counts.csv` which contains
all of the gRNA abundances for all of the screens, publications, cell lines, etc. in a 
single file. Finally, it will create a file called `data/preprocessed_screen_kmers.csv`
which contains the gRNA sequences tagged by ID and the corresponding 20-mer sequences.
This file can be read directly by Guidescan2 to compute specificity scores.

Next, we need to run MAGeCK on the preprocessed gRNA abundance files. This is done
using `snakemake`, which performs the analysis in parallel across the different
experimental conditions. To run the analysis, run the following command from the
root directory of this repository:

```
$ snakemake --cores all
```

This will create a new directory called `data/mageck_results` which contains the
MAGeCK output files for each of the screens. To summarize the MAGeCK results, run
the following commands:

```
$ python scripts/summarize_mageck_results.py results/mageck/*.gene_summary.txt --controls > results/mageck_results_table.with_controls.csv
$ python scripts/summarize_mageck_results.py results/mageck/*.nocontrols.gene_summary.txt > results/mageck_results_table.nocontrols.csv
```

These summarize the MAGeCK results for each of the screens, into two files.

Next, we need to run Guidescan2 on the preprocessed gRNA sequences. This is done
by directly invoking the Guidescan2 command line tool `guidescan` on the 
preprocessed gRNA kmers file using the appropriate genome index.

```
$ guidescan enumerate [GENOME_INDEX] -f data/preprocessed_screen_kmers.csv --format sam -o results/guidescan2_processed_grnas.sam -a NAG`
```

This will create a file called `results/guidescan2_processed_grnas.sam` which contains
the Guidescan2 specificity scores for each of the gRNAs in the preprocessed gRNA kmers
in a compact SAM format. To decode the SAM file into a CSV file, run the following
command:

```
$ python scripts/decode_sam.py results/guidescan2_processed_grnas.sam [GENOME_FASTA] --mode succinct > results/guidescan2_processed_grnas.csv
```

Finally, append the mean gene specificity scores to the MAGeCK results table using
the following command:

```
$ python scripts/append_mean_specificity.py results/mageck_results_table.with_controls.csv results/guidescan2_processed_grnas.csv > results/mageck_results_table.with_controls.guidescan2.csv
$ python scripts/append_mean_specificity.py results/mageck_results_table.nocontrols.csv results/guidescan2_processed_grnas.csv > results/mageck_results_table.nocontrols.guidescan2.csv
```
