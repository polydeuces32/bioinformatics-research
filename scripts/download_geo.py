#!/usr/bin/env python
import argparse
import hashlib
import logging
from pathlib import Path

import GEOparse


def compute_checksum(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def download(accession: str, out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    gse = GEOparse.get_GEO(geo=accession, destdir=str(out_dir), silent=True)

    samples = {name: gsm.table for name, gsm in gse.gsms.items() if not gsm.table.empty}
    if not samples:
        raise RuntimeError(f"No sample tables found for {accession}")

    first = next(iter(samples.values()))
    id_col = first.columns[0]
    value_col = "VALUE" if "VALUE" in first.columns else first.columns[1]

    matrix = None
    for name, table in samples.items():
        series = table.set_index(id_col)[value_col].rename(name)
        matrix = series.to_frame() if matrix is None else matrix.join(series, how="outer")

    out_path = out_dir / f"{accession}_counts.tsv"
    matrix.to_csv(out_path, sep="\t")
    return out_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Download a count matrix from GEO.")
    parser.add_argument("--accession", required=True, help="GEO accession (e.g. GSE12345)")
    parser.add_argument("--out-dir", default="data/raw", help="Destination directory")
    parser.add_argument("--log", default="logs/download_geo.log", help="Log file path")
    args = parser.parse_args()

    Path(args.log).parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        filename=args.log,
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )

    out_path = download(args.accession, Path(args.out_dir))
    checksum = compute_checksum(out_path)
    logging.info("Downloaded %s -> %s sha256=%s", args.accession, out_path, checksum)
    print(f"{out_path}\tsha256={checksum}")


if __name__ == "__main__":
    main()
