# Guidescan2 Performance Analysis

This repository contains the scripts, data, and instructions 
required to perform the performance analysis of Guidescan2
to other CRISPR-Cas9 guideRNA off-target enumeration tools.

## Instructions

### Downloading the data

To start, first download the `hg38` reference genome as a FASTA file. 
In particular, we will use the `GRCh38.p14` version of the genome, with 
RefSeq assembly ID `GCF_000001405.40`. Place the FASTA file, named
`GCF_000001405.40_GRCh38.p14_genomic.fna` in the directory `genomes/`
at the root of this repository. Create the `genomes/` directory if it 
does not already exist.

Next, to construct the indices for Crispritz, we will need to 
split the FASTA file into separate files for the individual scaffolds. 
To do this, we will use the `faSplit` tool from the UCSC Genome Browser. 
Once the tool is installed, run the following command to split the genome
into individual chromosomes:

```
$ mkdir genomes/chromosomes
$ faSplit byname genomes/GCF_000001405.40_GRCh38.p14_genomic.fna genomes/chromosomes/
```

### Construct and timing the indices

Create a directory `indexing_times/` at the root of this repository.
This directory will contain the timing information for the index
Next, construct the genome index for CRISPRitz. To do this,
first install [CRISPRitz](https://github.com/pinellolab/CRISPRitz). 
Then, run the following sequence of commands, which will create 
the index and time the process:

```
$ mkdir -p indices/crispritz
$ echo "NNNNNNNNNNNNNNNNNNNNNRG 3" > indices/crispritz/pam_NRG.txt
$ /usr/bin/time -v crispritz.py index-genome indices/crispritz/hg38 genomes/chromosomes/ \
                indices/crispritz/pam_NRG.txt -bMax 0 -th 4 > indexing_times/crispritz_timing.txt
```

Similarly, install [FlashFry](https://github.com/mckennalab/FlashFry) and run the following 
commands to construct the index and time the process. We assume that the `flashfry.jar` 
file is found in the current directory.

```
$ mkdir -p indices/flash_fry
$ /usr/bin/time -v java -Xmx4g -jar flashfry.jar index 
                              --tmpLocation ./tmp --database indices/flash_fry/hg38\
                              --reference genomes/GCF_000001405.40_GRCh38.p14_genomic.fna.gz --enzyme spcas9\
                              > flash_fry_timing.txt
```

Finally, install [GuideScan2](https://github.com/pritykinlab/guidescan-cli) and run 
the following commands to construct the index and time the process. We assume that 
the `guidescan` binary is available on the `PATH`.

```
$ /usr/bin/time -v guidescan index genomes/GCF_000001405.40_GRCh38.p14_genomic.fna > guidescan2_timing.txt
```

### Timing the off-target enumeration

The timing script will run the off-target enumeration for each tool
across a range of different parameters. To run the timing script,
first appropriately set the location of GuideScan2, FlashFry, and 
CRISPRitz at the top of the `Snakefile`:

```
crispritz_binary  = ''
guidescan2_binary = ''
flashfry_jar      = ''
```

Then run `snakemake` from the root of the
directory. Ensure that enough cores are allocated to ensure no
resource contention occurs.

```
$ snakemake -j 4
```

Once snakemaker has completed, the timing information will be
located in the `timing_results/` directory for each tool and parameter
setting. To summarize the results, run the following command:

```
$ python scripts/summarize.py --experiment-type normal timing_results/crispritz/ \
                              timing_results/flashfry/ timing_results/guidescan2/ \
                              > timing_results/summary.txt
``
`
### Checking that off-targets match across tools


