from Bio import SeqIO

import random
import gffutils
import os
import pandas as pd
import numpy as np
import argparse

NUCS = list("ACTG")
NUC_MAP = {"A": "T", "T": "A", "C": "G", "G": "C"}

def uniformly_random_guides(n, k = 20):
    for _ in range(n):
        yield "".join(random.sample(NUCS, 1)[0] for _ in range(k))

def revcom(dna):
    return "".join(list(map(lambda n: NUC_MAP[n], list(dna)))[::-1])

def generate_pam_set(pam):
    pam_stack = [pam]

    while any(["N" in pam for pam in pam_stack]):
        pam = pam_stack.pop(0)

        if "N" not in pam:
            pam_stack.append(pam)
            continue 

        for nuc in NUCS:
            pam_stack.append(pam.replace("N", nuc, 1))

    return pam_stack

def get_guides_across_region(chrm, start, end, pam, k=20, forward=True):
    index = start

    while True:
        index = chrm.find(pam, index)

        if index == -1 or index > end:
            break 

        if forward:
            kmer = chrm[index - k:index]
            position = index - k
        else:
            kmer = chrm[index + len(pam):index + k + len(pam)]
            position = index

        index += 1

        if position < 0:
            continue

        yield (kmer.upper(), position)

def get_random_guides(regions):
    for region in regions:
        chrm, start, end = region.seqid, region.start, region.end

        if args.chr2acc:
            if chrm not in acc2chrm:
                continue

        forward_guides = []
        for pam in generate_pam_set("NGG"):
            forward_guides += list(get_guides_across_region(record_dict[chrm].seq, start, end, pam))

        reverse_guides = []
        for pam in generate_pam_set("NGG"):
            reverse_guides += list(get_guides_across_region(record_dict[chrm].seq, start, end, revcom(pam), forward=False))

        forward_guides = [(x, y, '+') for (x, y) in forward_guides]
        reverse_guides = [(x, y, '-') for (x, y) in reverse_guides]

        guides = reverse_guides + forward_guides
        if not guides:
            continue

        if args.chr2acc:
            chrm = acc2chrm[chrm]

        guide = random.sample(guides, 1)[0]
        yield guide, chrm, region

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Selects gRNAs for performance analysis."
    )

    parser.add_argument(
        "reference_genome", help="FASTA reference genome."
    )

    parser.add_argument(
        "annotation_db", help="Genomic annotations in GFF format."
    )

    parser.add_argument(
        "analyzed_library", help="Guidescan2 analyzed gRNA library CSV."
    )

    parser.add_argument(
        "--chr2acc", help="Re-map accessions to chromosome names."
    )

    parser.add_argument(
        "-o", help="Output prefix"
    )

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()

    if args.chr2acc:
        acc2chrm = {}
        with open(args.chr2acc, "r") as f:
            next(f)
            for line in f:
                chrm, acc = line.split()
                acc2chrm[acc] = chrm

    random.seed(73)
    record_dict = SeqIO.to_dict(SeqIO.parse(args.reference_genome, "fasta"))
    db = gffutils.FeatureDB(args.annotation_db)

    lnc_regions = list(db.features_of_type("lnc_RNA"))
    random_lnc_regions = random.sample(lnc_regions, 10_000)

    random_lnc_guides = []
    for guide, chrm, lnc_region in get_random_guides(random_lnc_regions):
        random_lnc_guides.append({
            "Identifier": f"{lnc_region.attributes['gene'][0]}:{chrm}:{guide[1]}",
            "gRNA": str(guide[0]),
            "Position": str(guide[1]),
            "Chromosome": chrm,
            "Gene": lnc_region.attributes["gene"][0],
            "PAM": "NGG",
            "Strand": guide[2]
        })

    random_lnc_df = pd.DataFrame(random_lnc_guides)

    mirna_regions = list(db.features_of_type("miRNA"))
    random_mirna_regions = random.sample(mirna_regions, 1_000)

    random_mirna_guides = []
    for guide, chrm, mirna_region in get_random_guides(random_mirna_regions):
        random_mirna_guides.append({
            "Identifier": f"{mirna_region.attributes['gene'][0]}:{chrm}:{guide[1]}",
            "gRNA": str(guide[0]),
            "Position": str(guide[1]),
            "Chromosome": chrm,
            "Gene": mirna_region.attributes["gene"][0],
            "PAM": "NGG",
            "Strand": guide[2]
        })

    random_mirna_df = pd.DataFrame(random_mirna_guides)

    cds_regions = list(db.features_of_type("CDS"))
    random_cds_regions = random.sample(cds_regions, 10_000)

    random_cds_guides = []
    for guide, chrm, cds_region in get_random_guides(random_cds_regions):
        random_cds_guides.append({
            "Identifier": f"{cds_region.attributes['gene'][0]}:{chrm}:{guide[1]}",
            "gRNA": str(guide[0]),
            "Position": str(guide[1]),
            "Chromosome": chrm,
            "Gene": cds_region.attributes["gene"][0],
            "PAM": "NGG",
            "Strand": guide[2]
        })

    random_cds_df = pd.DataFrame(random_cds_guides)

    unif_random_guides = [
        {"Identifier": f"{g}_{i}", "gRNA": g, "PAM": "NGG"} for (i, g) in enumerate(uniformly_random_guides(10_000))
    ]

    analyzed_library_df = pd.read_csv(args.analyzed_library)
    analyzed_library_df = analyzed_library_df[analyzed_library_df['chr'].isin(acc2chrm)]
    analyzed_library_df = analyzed_library_df.drop(
        columns=["Unnamed: 0", "score", "absolute_pos", "distance"]
    ).rename(
        columns={"identifier": "Identifier",
                 "gene": "Gene",
                 "sgRNA": "gRNA",
                 "strand": "Strand",
                 "chr": "Chromosome",
                 "pos": "Position",
                 "specificity": "Specificity",
                 "cutting_efficieny": "Cutting Efficieny"}
    )

    analyzed_library_df["Chromosome"] = analyzed_library_df["Chromosome"].apply(
        lambda acc: acc2chrm[acc]
    )

    analyzed_library_df.sample(10000).to_csv(f"{args.o}_published_guides.csv")

    for n in [20, 19, 18, 17]:
        analyzed_library_df[analyzed_library_df['gRNA'].str.len() == n]\
            .sample(10000)\
            .to_csv(f"{args.o}_published_n{n}_guides.csv")

    unif_random_guides_df = pd.DataFrame(unif_random_guides)
    random_mirna_df.to_csv(f"{args.o}_mirna_guides.csv")
    random_lnc_df.to_csv(f"{args.o}_lnc_guides.csv")
    random_cds_df.to_csv(f"{args.o}_cds_guides.csv")
    unif_random_guides_df.to_csv(f"{args.o}_unif_guides.csv")
