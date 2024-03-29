import sys
import argparse
import pandas

def parse_args():
    parser = argparse.ArgumentParser(description="Attach scores to library file")
    parser.add_argument("scores1", help="Scores file")
    parser.add_argument("scores2", help="Scores file")
    parser.add_argument("input", help="Input file")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    scores1_df = pandas.read_csv(args.scores1, header=None, names=["id", "specificity"])
    scores2_df = pandas.read_csv(args.scores2, header=None, names=["id", "specificity"])
    scores_df = pandas.concat([scores1_df, scores2_df])
    library_df = pandas.read_csv(args.input)

    merged_df = pandas.merge(library_df, scores_df, on="id", how="left")
    merged_df.to_csv(sys.stdout, index=False)
