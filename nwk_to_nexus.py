import re
from typing import List, Union
from io import TextIOWrapper
from ete3 import Tree


def load_tree(file: str, out: str):
    tree = Tree(newick=file)
    tree_txt = ""
    with open(file, "r") as f:
        tree_txt = f.read()
    taxons = []
    for leaf in tree.iter_leaf_names():
        taxons.append(leaf)
    write_to_file(taxons, out, 1)
    write_to_file(tree_txt, out, 2)
    with open("./fig_tree_settings.txt", "r") as f:
        settings = f.read()
    write_to_file(settings, out, 3)


def write_to_file(contents: Union[str, List], out: str, index: int):
    with open(out, "a") as f:
        if index == 1:
            contents.sort()
            f.write(
                f"#NEXUS\nbegin taxa;\n\tdimensions ntax={len(contents)};\n\ttaxlabels\n"
            )
            for content in contents:
                f.write(f"\t{content[1:-1]}\n")
            f.write(f";\nend;\n\n")
        elif index == 2:
            f.write(f"begin trees;\n\ttree tree_1 = [&R] ")
            write_nexus_with_labels(contents, f)
            f.write(f"\nend;\n\n")
        elif index == 3:
            f.write(f"{contents}")


def write_nexus_with_labels(tree_txt: str, out_f: TextIOWrapper):
    last_loc = 0
    for i, match in enumerate(re.finditer(r"'|\)\d+", tree_txt)):
        out_f.write(f"{tree_txt[last_loc:match.start()]}")
        if tree_txt[match.start() : match.end()] != "'":
            out_f.write(f"{tree_txt[match.start():match.start()+ 1]}")
            out_f.write(f"[&label={tree_txt[match.start() +1: match.end()]}]")
        last_loc = match.end()


load_tree(
    r"D:\Datasets\GISAID_Update_Analysis\nwk_conversion_sample\sample\fnb_all_rapidnj.nwk",
    r"D:\Datasets\GISAID_Update_Analysis\nwk_conversion_sample\sample\test_out\test",
)
