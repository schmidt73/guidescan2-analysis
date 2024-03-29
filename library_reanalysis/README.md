# CRISPRko gRNA library reanalysis

We first pulled the following gRNA libraries from (addgene)[www.addgene.org]. They
are located in `inputs/mouse/` and `inputs/human` and are summarized in the 
table below.

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

Then, we pre-process the gRNA libraries into a unified format
to be used by the `guidescan` tool. In particular, we tag each gRNA
with the library name, species, target gene, unique identifier, and 
sequence. To do this, we run the following script:

```bash
$ python scripts/preprocess.py
                           count      mean       std  min  25%  50%  75%    max
library         species
guidescan2      human    18232.0  5.911200  0.604103  1.0  6.0  6.0  6.0    6.0
                mouse    20908.0  5.864167  0.714533  1.0  6.0  6.0  6.0    6.0
liu             mouse    18705.0  4.965410  1.887652  1.0  5.0  5.0  5.0  243.0
minimized_dgrna human    19256.0  3.054217  0.265119  1.0  3.0  3.0  3.0    8.0
minlabcas9      human    17389.0  1.954339  0.208755  1.0  2.0  2.0  2.0    2.0
teichmann       mouse    18424.0  4.897416  0.744960  1.0  5.0  5.0  5.0   17.0
tkov3           human    18053.0  3.929984  0.404214  1.0  4.0  4.0  4.0    4.0
```

This computes the number of gRNAs in each library and the mean, standard deviation,
and quartiles of the number of gRNAs per gene. It will also create three output files in 
the `outputs/` directory:
- `all_libraries.csv`: all gRNAs from all libraries
- `all_libraries_mouse_kmers.csv`: mouse gRNAs from all libraries in Guidescan2 kmer format
- `all_libraries_human_kmers.csv`: human gRNAs from all libraries in Guidescan2 kmer format

Now, we are ready to run `guidescan`:
```bash
$ guidescan enumerate /n/fs/ragr-research/projects/guidescan2/genomes/GCF_000001405.40_GRCh38.p14_genomic.fna.index \
                      --kmers-file outputs/all_libraries_human_kmers.csv --alt-pam NAG --mismatches 3 \
                      --output outputs/all_libraries_human_kmers_guidescan.sam  --format sam
```

