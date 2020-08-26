import argparse
import tarfile


def extract(file: str, out: str):
    with tarfile.open(file, "r") as t:
        t.extractall(out)


def _parse_arg():
    parser = argparse.ArgumentParser(description=".tar file extractor")
    parser.add_argument(
        "input",
        metavar="I",
        type=str,
        help="The path to the input .tar file to extract",
    )
    parser.add_argument(
        "output",
        metavar="O",
        type=str,
        help="""The path of the output folder to extract items into. 
        If folder does not exist, a new one will automatically be created
        E.g. ./out folder out will automatically be created if it doesn't exist""",
    )
    return parser.parse_args()


def main():
    args = _parse_arg()
    extract(args.input, args.output)


if __name__ == "__main__":
    main()

