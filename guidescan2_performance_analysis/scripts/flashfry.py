import argparse
import pandas as pd
import numpy as np

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Converts gRNAs to a format suitable for the FlashFry tool."
    )

    parser.add_argument(
        "grnas_csv", help="gRNAs to sub-sample from."
    )

    parser.add_argument(
        "-n", type=int, help="Number of gRNAs to select.", default=1000
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()
    
    np.random.seed(73)

    df = pd.read_csv(args.grnas_csv)
    grnas = df.sample(args.n)

    # need to generate PAM set manually
    for (_, grna) in grnas.iterrows():
        g, identifier, pam = grna.loc["gRNA"], grna.loc["Identifier"], grna.loc["PAM"]
        print(f">{identifier}_TGG")
        print(f"{g}TGG")
