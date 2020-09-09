import re
import argparse
import datetime
from typing import List, Union
from io import TextIOWrapper
from ete3 import Tree, TreeNode


def convert_to_nexus(file: str, out: str, setting_path: str):
    out += f'{datetime.date.today().strftime("%Y-%m-%d")}_Full'
    taxon_out = get_taxons_out(file)
    write_to_file(taxon_out, out, write=True)
    tree_txt = get_figtree(file)
    write_to_file(tree_txt, out)
    settings = get_settings(setting_path)
    write_to_file(settings, out)


def get_taxons_out(file: str):
    def build_taxon_format(taxons: List[str]):
        taxons.sort()
        out = f"#NEXUS\nbegin taxa;\n\tdimensions ntax={len(taxons)};\n\ttaxlabels\n"
        for taxon in taxons:
            out += f"\t{taxon[1:-1]}\n"
        out += f";\nend;\n\n"
        return out

    def get_taxons(tree: TreeNode):
        taxons = []
        for leaf in tree.iter_leaf_names():
            taxons.append(leaf)
        return taxons

    tree = Tree(newick=file,)
    taxons = get_taxons(tree)
    output = build_taxon_format(taxons)
    return output


def get_figtree(file: str):
    with open(file, "r") as f:
        tree_txt = f.read()
    result = ""
    last_loc = 0
    for i, match in enumerate(re.finditer(r"'|\)\d+", tree_txt)):
        result += f"{tree_txt[last_loc:match.start()]}"
        if tree_txt[match.start() : match.end()] != "'":
            result += f"{tree_txt[match.start():match.start()+ 1]}"
            result += f"[&label={tree_txt[match.start() +1: match.end()]}]"
        last_loc = match.end()
    result += tree_txt[last_loc:]
    result = f"begin trees;\n\ttree tree_1 = [&R] " + result + f"end;\n\n"
    return result


def get_settings(setting_path: str):
    with open(setting_path, "r") as f:
        settings = f.read()
    return settings


def write_to_file(contents: str, out: str, write: bool = False):
    mode = "w" if write else "a"
    with open(out, mode) as f:
        f.write(contents)


def _parse_arg():
    """Command line arguments"""
    parser = argparse.ArgumentParser(
        description="A convertor from .nwk file to figtree output file while preserving color labels and strips quoted tree label",
        epilog="If you notice any issues, please open one over at https://github.com/ElasticBottle/GISAID-Analysis-Update ",
    )
    parser.add_argument(
        "--input",
        "-i",
        type=str,
        dest="i",
        default="./input.nwk",
        help="The path to the input file. Remember to include filename and extension! Defaults to ./input.nwk",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        dest="o",
        default=f"./",
        help="The path to the output file. Should have no extensions. Deafult to today's date + '_Full'",
    )
    parser.add_argument(
        "--settings",
        "-s",
        type=str,
        dest="setting_path",
        default="HOME_ann/BII/biipsashare/winston/GISAID-Analysis-Update/fig_tree_settings.txt",
        help="The path to the setting file on your machine",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_arg()
    convert_to_nexus(args.i, args.o, args.setting_path)
