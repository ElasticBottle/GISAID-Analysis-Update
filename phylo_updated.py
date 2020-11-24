import argparse
import os
import json
import time
from functools import partial
from typing import Dict, List, Tuple, Union

import numpy as np
import pandas as pd

from ete3 import NodeStyle, Tree, TreeNode, TreeStyle

ALT_NAME = 0
CLADE_COL= 0
MIN_DENSITY = 1
MAX_DENSITY= 2
MIN_COV = 3
MAX_COV = 4

def tree_style() -> TreeStyle:
    style = TreeStyle()
    style.mode = "c"  # draw tree in circular mode
    # style.scale = 20
    # style.branch_vertical_margin = 15
    style.root_opening_factor = 0.04
    style.show_leaf_name = False
    style.show_branch_length = False
    style.show_branch_support = False
    return style


def node_style(color: str = None) -> NodeStyle:
    style = NodeStyle()
    style["fgcolor"] = "#0f0f0f"
    style["size"] = 0
    style["shape"] = "circle"
    style["vt_line_color"] = color if color is not None else "#343434"
    style["hz_line_color"] = color if color is not None else "#343434"
    style["vt_line_type"] = 0  # 0 solid, 1 dashed, 2 dotted
    style["hz_line_type"] = 0
    style["vt_line_width"] = 2 if color is None else 0
    style["hz_line_width"] = 2 if color is None else 0
    return style


def node_bg_style(style: Dict, bg_color: str = None):
    style["bgcolor"] = bg_color
    return style


def _parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Generates hcov-19 Phylogentic tree based on GISAID bi-weekly analysis update",
        epilog="If you notice any issues, please open one over at https://github.com/ElasticBottle/GISAID-Analysis-Update ",
    )
    parser.add_argument(
        "file",
        type=str,
        default="",
        help="file path of the NEWICK tree to parse",
    )
    parser.add_argument(
        "-o",
        "--out",
        type=str,
        dest="out",
        metavar="",
        default="./out.png",
        help="The output path of the image to be stored. Default ./out.png. Script assumes extension is 3 letters long, .png, .svg etc.",
    )
    parser.add_argument(
        "--dpi",
        default=300,
        type=int,
        dest="dpi",
        help="Dpi of the output image. Default 300",
    )
    parser.add_argument(
        "-w",
        "--width",
        default=5000,
        type=int,
        dest="width",
        help="Specifies the width of the image in pixel, the height is then automatically derived. Default 5000",
    )
    parser.add_argument(
        "-d",
        "--display",
        action="store_false",
        dest="no_display",
        help="Use this flag if you have a display",
    )
    parser.add_argument(
        "-q",
        "--quoted",
        action="store_true",
        dest="is_quoted",
        help="Use this flag if newwick file contains quotes around the names",
    )
    parser.add_argument(
        "-ccov",
        "--clade-color-density-coverage",
        type=str,
        dest="clade_color_density_coverage",
        help="path to json file mapping clade delineators to a list[colour(String), density(float), coverage(float)] values. Clade Delineators should occur at the start of sequence input",
    )
    parser.add_argument(
        "-off",
        "--offset",
        dest="offset",
        default = -1,
        type=int,
        help="Number of characters to skip from the back when extracting color code",
    )
    parser.add_argument(
        "-ro",
        "--root-on",
        dest="root_on",
        type=str,
        help="Specifies the clade to root the tree on. Value should occur at the start of sequence name",
    )
    parser.add_argument(
        "-cm",
        "--color-marker",
        dest="color_marker",
        type=str,
        default="",
        help="""Specify the marker to split the string on so that it can obtain the a list in [CLADE_NAME, COLOR] format. 

        Color should be in hex code (either FFFFF or #FFFFF) and is the last thing in sequence name, right after the colo-marker""",
    )
    parser.add_argument(
        "-rf",
        "--rapid-fasttree-color-marker",
        dest="color_marker",
        const="color=",
        action ="store_const",        
        help="Sets the marker to split for rapidnj and fasttree for extracting color information. DO NOT use -cm or -iq with this option",
    )
    parser.add_argument(
        "-iq",
        "--iqtree-color-marker",
        dest="color_marker",
        const="____color__",
        action="store_const",
        help="Sets the marker to split on for iqtree for extracting color information. DO NOT use -cm or -rf with this option",
    )
    args = parser.parse_args()
    return args


def main():
    args = _parse_args()
    f_slash = args.file.rfind("/")
    b_slash = args.file.rfind("\\")
    print(
        f"generating image for {args.file[f_slash + 1 if f_slash != -1 else b_slash + 1:]}, this will take a while. ~2-5 mins"
    )
    if args.no_display:
        os.environ["QT_QPA_PLATFORM"] = "offscreen"
    
   
    # mapping of clade names to either: Alternate name || [colour, density, min_coverage, max_coverage]
    clade_details = {}
    with open(args.clade_color_density_coverage) as f:
        clade_details = json.load(f)
    alt_names = {clade: desired_name[ALT_NAME] for clade, desired_name in clade_details.items() if len(desired_name) == 1}
    clade_list_details = {clade: col_dens_cov for clade, col_dens_cov in clade_details.items() if len(col_dens_cov) > 1}
    clades = list(clade_list_details.keys())

    start = time.time()
    tree = get_tree(file=args.file, is_quoted=args.is_quoted)
    print(f"{len(tree.get_leaves())} taxons present")

    color_taxons(tree, args.color_marker, args.offset)
    clade_counts = get_clade_count(tree, clades = clades,alternate_names=alt_names)
    tree = root_on(
        tree = tree, 
        clade = args.root_on,
        clade_details = clade_list_details[args.root_on], 
        clade_total = clade_counts[args.root_on],
        clades = clades,
        alternate_names=alt_names,
    )
    node_cache = get_node_details(tree, clades = clades, alternate_names = alt_names)
    color_clades(tree, clades = clades, clade_details = clade_list_details, clade_counts=clade_counts, node_cache = node_cache)
    tree.render(args.out, dpi=args.dpi, w=args.width, tree_style=tree_style())
    end = time.time()

    print(f"Time taken: {end - start:.2f}s")


def get_tree(file: str, is_quoted: bool) -> TreeNode:
    """
    Creates a newick tree where leafs are colored base on the tag in name, and background base on the clade grouping.

    Args:
        - file(str): the path to the newick tree
    """
    tree = Tree(newick=file, quoted_node_names=is_quoted)
    return tree

def color_taxons(tree: TreeNode, color_marker: str, offset:int):
    """
    Sets the node style for all nodes in [tree]

    Args:

        - tree(TreeNode): the input tree to set node_style for
        - color_marker(Str): the delineator to split taxon sequence name to extract 
            color information. Color info should be the last thing right after [color_marker].

            Color should also be in Hex format either FFFFF or #FFFFFF
        - offset(int): how many characters off the back to chop off.

    """
    for node in tree.traverse():
        node.img_style = node_style()
        if node.is_leaf():
            color_taxon(node, color_marker=color_marker, offset= offset)

def color_taxon(node:TreeNode, color_marker: str, offset: int):
    """
    Sets the NodeStyle for [n]

    Args:

        - node(TreeNode): the leaf to set node_style for
        - color_marker(Str): the delineator to split taxon sequence name to extract 
            color information. Color info should be the last thing right after [color_marker].

            Color should also be in Hex format either FFFFF or #FFFFFF
        - offset(int): how many characters off the back to chop off.
    """
    split = node.name.split(color_marker)
    if len(split) > 1:
        color = split[-1][:offset]
        if color.startswith("#"):
            node.img_style = node_style(color=color)
        else:
            node.img_style = node_style(color= f"#{color}")

def get_clade_count(tree:TreeNode, clades: List[str], alternate_names: Dict[str, str]) -> Dict[str, int]:
    """
    Returns the total number of clades for each clade in [clades].

    Sequence not matching any of [clades] will be added under "other" in return dictionary

    Args:

        - tree (TreeNode): the tree in which clades are to be counted for
        - clades (List[str]): List containing the clade names
            Clade names should occur at the start of the sequence name.
        - altername_names (Dict[str, str]): Maps alternate name found in sequence to desired name.

    Returns:

        - Dict[str, int]: Maps the clade to the total number present.
    """
    def get_belonging(leaf: TreeNode, clades: Dict[str, str]):
        """
        """
        result = [clades[clade] for clade in clades.keys() if leaf.name.startswith(clade)]
        if len(result) == 0:
            return "other"
        return result[0]
    

    clade_dict = {clade: clade for clade in clades}
    clade_dict.update(alternate_names)

    leaf_clades = list(map(partial(get_belonging, clades = clade_dict), tree.get_leaves()))
    serialize = pd.Series(data = leaf_clades, dtype = str)
    return serialize.value_counts().to_dict()

        
def root_on(tree: TreeNode, clade: str, clade_details: List, clade_total:int, clades: List[str], alternate_names: Dict[str, str]) -> TreeNode:
    """
    Roots the given tree on clade

    Args:

        - tree (TreeNode): The tree in which will be rooted
        - clade (str): The clade to root tree on. Tree must contain clade.
        - clade_details (List): Contains [color, density, coverage] in that order
        - clade_total (int): The total number leaves in the tree belonging to [clade]
        - clades (List[str]): List containing the clade names
            Clade names should occur at the start of the sequence name.
        - altername_names (Dict[str, str]): Maps alternate name found in sequence to desired name.
        
    Returns:        

        - TreeNode: A tree with clade as the outgroup
    """
    node = get_max_ancestor(tree, clade, clade_details, clade_total, clades, alternate_names)
    tree.set_outgroup(node)
    return tree

def get_node_details(tree: TreeNode, clades:List[str], alternate_names: Dict[str, str]) -> Dict[TreeNode, Dict[str, int]]:
    """
    Computes the clade distribution at each node of the tree.

    Args

        - tree (TreeNode): The tree in which will be rooted
        - clades (List[str]): List containing the clade names
            Clade names should occur at the start of the sequence name.
        - altername_names (Dict[str, str]): Maps alternate name found in sequence to desired name.
    
    Returns:

        - Dict[TreeNode, Dict[str, int]]: DIctionary mapping a node to the clade_distribution under it.
    """
    result = {}
    for node in tree.traverse():
        if node.is_leaf():
            continue
        result[node] = get_clade_count(node, clades= clades, alternate_names=alternate_names)
    return result

def get_max_ancestor(tree: TreeNode, clade: Union[str, List[str]], clade_details: List, clade_total:int, clades: List[str] = None, alternate_names: Dict[str, str] = None, node_cache: Dict = None):
    """
    Finds the best ancestor that fulfills the clade_details
    
    Args:

        - tree (TreeNode): The tree in which will be rooted
        - clade (str): The clade to root tree on. Tree must contain sequence belonging to clade.
        - clade_details (List): Contains [color, density, min_coverage, max_coverage] in that order
        - clade_total (int): The total number leaves in the tree belonging to [clade]
        - clades (List[str]): List containing the clade names
            Clade names should occur at the start of the sequence name.
            Should be specified if [node_cache] is not specified
        - alternate_names (Dict[str, str]): Maps alternate name found in sequence to desired name.
            Should be specified if [node_cache] is not specified
        - node_cache (Dict): Maps TreeNode to the clade distribution of its descnedants.
            Should be specified if you are coloring the tree

    Returns:

        - TreeNode: The node corresponding to the best fit ancestor for the given clade
    """
    def get_details(distribution: Dict[str, int], clade:str, clade_total: int) -> Tuple[float, float]:
        total = sum(distribution.values())
        clade_num = float(distribution.get(clade, 0))
        return (clade_num / total, clade_num / clade_total)
    
    def get_qualifying_nodes (tree: TreeNode,clade: Union[str, List[str]], clade_total: int, clade_details: List, clades: List[str] = None, alternate_names: Dict[str, str] = None, node_cache = None, ):
        qualifying_list = {}
        for node in tree.traverse():
            if node.is_leaf():
                continue
            
            children_distribution = {}
            if node_cache is None:
                children_distribution = get_clade_count(node, clades= clades, alternate_names=alternate_names)
            else:
                children_distribution = node_cache[node]
            density, coverage = get_details(children_distribution, clade, clade_total)
            if density >= clade_details[MIN_DENSITY] and density <= clade_details[MAX_DENSITY] and coverage >= clade_details[MIN_COV] and coverage <= clade_details[MAX_COV]:
                qualifying_list[node] = (density, coverage, density + coverage)
        return qualifying_list
        
    SCORE = 2

    qualifying_list = get_qualifying_nodes(tree, clade= clade, clade_details= clade_details, clade_total= clade_total, clades = clades, alternate_names= alternate_names, node_cache= node_cache)
    
    sorted_nodes = list(reversed(sorted(qualifying_list.keys(), key=lambda x: qualifying_list.get(x)[SCORE])))
    try:
        print(
                f"{clade} coverage: {(qualifying_list[sorted_nodes[0]][1])* 100:.2f}%, density of coverage: {qualifying_list[sorted_nodes[0]][0] * 100:.2f}%"
            )
    except:
        raise Exception(f"Could not find any node to color for clade {clade} that fulfills criteria: {clade_details} ")
    return sorted_nodes[0]

def color_clades(tree: TreeNode, clades:List[str], clade_details: Dict[str, List], clade_counts: Dict[str, int], node_cache: Dict[TreeNode, Dict[str, int]]):
    """
    Colors a single node for each item in [clades]

    Args:

        - tree (TreeNode): The tree in which will be rooted
        - clade (str): The clade to root tree on. Tree must contain sequence belonging to clade.
        - clade_cache (Dict): Maps TreeNode to the clade distribution of its descnedants.
    """
    for clade in clades:
        max_ancestor = get_max_ancestor(tree, clade = clade, clade_details=clade_details[clade], clade_total=clade_counts[clade], node_cache = node_cache)
        max_ancestor.img_style = node_bg_style(max_ancestor.img_style, clade_details[clade][CLADE_COL])
if __name__ == "__main__":
    main()
