#!/usr/bin/env python
import argparse
import hashlib
import logging
from pathlib import Path

import numpy as np
import pandas as pd

RANDOM_SEED = 42


def compute_checksum(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def preprocess(input_path: Path, output_path: Path) -> None:
    np.random.seed(RANDOM_SEED)

    counts = pd.read_csv(input_path, sep="\t", index_col=0)
    counts = counts.apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)
    counts = counts.loc[counts.sum(axis=1) > 0]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    counts.to_csv(output_path, sep="\t")


def main() -> None:
    parser = argparse.ArgumentParser(description="Preprocess a raw count matrix.")
    parser.add_argument("input", type=Path, help="Path to raw count matrix (TSV)")
    parser.add_argument("output", type=Path, help="Path to write processed matrix (TSV)")
    parser.add_argument("--log", default="logs/preprocess.log", help="Log file path")
    args = parser.parse_args()

    Path(args.log).parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        filename=args.log,
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )

    preprocess(args.input, args.output)
    checksum = compute_checksum(args.output)
    logging.info("Preprocessed %s -> %s sha256=%s", args.input, args.output, checksum)
    print(f"{args.output}\tsha256={checksum}")


if __name__ == "__main__":
    main()
