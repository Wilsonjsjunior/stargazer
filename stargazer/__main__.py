import argparse
import datetime
import sys
import timeit
import logging

from .version import __version__
from .genotype import genotype

def get_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "dt",
        help="data type (wgs, ts, chip)"
    )
    parser.add_argument(
        "gb",
        help="genome build (hg19, hg38)"
    )
    parser.add_argument(
        "tg",
        help="target gene"
    )
    parser.add_argument(
        "vcf",
        help="VCF file"
    )
    parser.add_argument(
        "out",
        help="output project directory"
    )
    parser.add_argument(
        "--cg",
        help="control gene or region",
        metavar="STR"
    )
    parser.add_argument(
        "--gdf",
        help="GDF file",
        metavar="FILE"
    )
    parser.add_argument(
        "--ref",
        help="reference VCF file",
        metavar="FILE"
    )
    parser.add_argument(
        "--sl",
        nargs="*",
        help="sample list",
        metavar="STR"
    )
    parser.add_argument(
        "--dp",
        action="store_true",
        help="output more detailed plots"
    )
    parser.add_argument(
        "--imp",
        action="store_true",
        help="impute ungenotyped markers"
    )

    return parser

def main():
    start_time = timeit.default_timer()

    parser = get_parser()
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    logger.info(f"Stargazer v{__version__}")
    logger.info(f"""Command: '{" ".join(sys.argv)}'""")
    logger.info(f"Time: {datetime.datetime.now()}")

    genotype(
        args.dt,
        args.gb,
        args.tg,
        args.vcf,
        args.out,
        args.cg,
        args.gdf,
        args.ref,
        args.sl,
        args.dp,
        args.imp,
    )

    stop_time = timeit.default_timer()
    elapsed_time = str(
        datetime.timedelta(seconds=(stop_time - start_time))
    ).split(".")[0]

    logger.info(f"Elapsed time: {elapsed_time}")
    logger.info("Stargazer finished")

if __name__ == "__main__":
    main()