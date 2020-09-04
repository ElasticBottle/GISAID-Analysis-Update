import argparse
import time
from collections import defaultdict
from typing import Dict, Iterable, List, Tuple, Union

import numpy as np

from ete3 import NodeStyle, Tree, TreeNode, TreeStyle


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


def node_bg_style(style: Dict, clade: str, bg_color: str = None):
    color_scheme = {
        "G_": "#fdb3b3",
        "GH": "#feab8d",
        "GR": "#fe8d8d",
        "S_": "#b3ffb3",
        "V_": "#f0b3f0",
        "L_": "White",
        "O_": "White",
    }
    style["bgcolor"] = (
        bg_color if bg_color is not None else color_scheme.get(clade, None)
    )
    return style


class PhyloIMG:
    def __init__(self, file: str):
        self.file = file
        self.clade_list = [
            "G_",
            "GH",
            "GR",
            "S_",
            "V_",
            "L_",
            "O_",
        ]
        self.clades = {
            "G_": [],
            "GH": [],
            "GR": [],
            "S_": [],
            "V_": [],
            "L_": [],
            "O_": [],
        }  # clade -> taxons
        self.leaf_clade_group = {}  # taxon_name -> clade CURRENTLY NOT USED
        self.colored_nodes = (
            []
        )  # list of internal nodes whose bgcolor has been changed CURRENTLY NOT USED
        self.nodes_per_clade = {}  # internal_node -> number of leaves in each clade
        self.tree = self.get_tree(self.file)

    def get_tree(self, file: str) -> Tree:
        """
        Creates a newick tree where leafs are colored base on the tag in name, and background base on the clade grouping.

        Args:
            - file(str): the path to the newick tree
        """
        tree = Tree(newick=file)
        self._set_node_style(tree)
        self._root_on_s(tree)
        self._color_clades(tree)
        return tree

    def _set_node_style(self, tree: TreeNode):
        """
        Sets the node style for all nodes in the tree
        """
        for n in tree.traverse():
            n.img_style = node_style()
            if n.is_leaf():
                self._color_node(n)

    def _color_node(self, n: TreeNode):
        """
        Sets the NodeStyle for [n] and updates the clades dictionary
        """
        clade, color = self._get_clade_color(n.name)
        n.img_style = node_style(color)
        self._update_clades(clade, n)

    def _get_clade_color(self, name: str) -> Tuple[str, Union[str, None]]:
        """
        Args:
            - name(str): In one of the following formats ('' is part of the string):
                'GR_USAMNMDH1413_20200623_NorthAmerica_[&!color=ff00cc]'
                'GR_EnglandLONDD4621_20200408_Europe_'
        Returns:
            - str: the clade that the taxonomy belongs too
            - str / None : the color for that clade if any
        """
        name_color = name.split("color=")
        clade = name_color[0][1:3]
        color = None
        if len(name_color) == 2:
            color = name_color[1][:-2]
        return (clade, color)

    def _update_clades(self, clade: str, node: TreeNode):
        """
        Update [self.clades] and [self.leaf_clade_group] with the taxon and clade respectively
        
        Args:
            - clade(str): The clade of the taxonomy
            - node(TreeNode): The node to be updated
        """
        clade_group = self.clades.get(clade, self.clades["O_"])
        clade_group.append(node.name)
        if clade in self.clade_list:
            self.clades[clade] = clade_group
            self.leaf_clade_group[node.name] = clade
        else:
            self.clades["O_"] = clade_group
            self.leaf_clade_group[node.name] = "O_"

    def _root_on_s(self, tree: TreeNode):
        s_ancestor, _ = self.max_ancestor(
            tree, "S_", len(self.clades["S_"]), coverage=0.6, rooting=True
        )
        tree.set_outgroup(s_ancestor)

    def find_num_nodes(self, node: TreeNode) -> Dict[str, int]:
        """
        Finds the total number of taxons in each clade from [node]

        Args:
            - node(TreeNode): The node assumed to be the root
        """
        num_clade = {}
        for leaf in node.iter_leaf_names():
            clade, _ = self._get_clade_color(leaf)
            num_clade[clade] = num_clade.get(clade, 0) + 1
        return num_clade

    def calculate_density(self, node_per_clade: Dict[str, int], clade: str) -> float:
        """
        Calculate the percentage of [clade] for a given [node]
        """
        node_count = np.sum(list(node_per_clade.values()))
        # print(
        #     node_per_clade, node_per_clade[clade] / node_count * 100,
        # )
        return node_per_clade[clade] / node_count * 100

    def max_ancestor(
        self,
        tree: TreeNode,
        clade: str,
        clade_total: int,
        coverage: float = 0.6,
        rooting: bool = False,
    ) -> Tuple[TreeNode, float]:
        """
        Finds the internal node which captures [coverage] of taxons from [clade] and returns the node with the highest density of taxons from [clade]
       
        Example:
        Tree has clade distribution of:
        {'GH': 947, 'O_': 304, 'G_': 887, 'S_': 186, 'GR': 1402, 'L_': 77, 'V_': 108}

        For clade GH, with coverage of 0.5, find the internal nodes where:
            - 473.5 < Number of taxons belonging to GH < 947 (clade total)
        Return the internal node where proportion of GH to other clade is the highest
        """
        min_size = coverage * clade_total
        node_details = {}  # node -> density of [clade]
        for leaf in tree.traverse():
            if leaf.is_leaf() or leaf.is_root():
                continue

            # Saving calculations to prevent extra calculations
            node_per_clade = self.nodes_per_clade.get(leaf, None)
            if node_per_clade is None:
                node_per_clade = self.find_num_nodes(leaf)
                if not rooting:
                    self.nodes_per_clade[leaf] = node_per_clade

            # Finding internal nodes that satisfy coverage criteria
            if node_per_clade.get(clade, 0) >= min_size:
                density = self.calculate_density(node_per_clade, clade)
                coverage = node_per_clade[clade] / clade_total * 100
                node_details[leaf] = (density, coverage, node_per_clade)

        # Retrieving the node with densest collection of [clade]
        sorted_nodes = sorted(node_details.keys(), key=node_details.get)
        max_node = sorted_nodes[-1]
        for node in reversed(sorted_nodes):
            if (
                node_details[node][0] > 90
                and node_details[node][1] > node_details[max_node][1]
            ):
                max_node = node
        return (max_node, node_details[max_node])

    def _base_clade_color(self, tree, clade: str):
        """
        Colors a single node for the specified clade

        Args:
            - tree(TreeNode): the tree whose background is to be coloured
            - clade(str): the clade to color
        """
        coverage = (
            0.9
            if clade == "G_"
            else 0.85
            if clade == "GR"
            else 0.84
            if clade == "GH"
            else 0.6
        )
        clade_total = len(self.clades[clade])
        node, density_coverage_breakdown = self.max_ancestor(
            tree, clade=clade, clade_total=clade_total, coverage=coverage
        )
        node.img_style = node_bg_style(node.img_style, clade)
        print(
            f"{clade} coverage: {(density_coverage_breakdown[1]):.2f}%, density of coverage: {density_coverage_breakdown[0]:.2f}%"
        )
        self.colored_nodes.append(node)

    def _color_clades(self, tree: TreeNode):
        """
        Colors the 5 major clades (G, GH, GR, S, V) in the tree

        Args:
            tree(TreeNode) - The tree whose background colors are to be set
        """
        print(f"{len(self.leaf_clade_group.keys())} taxons present")
        for clade in self.clade_list:
            if clade == "L_" or clade == "O_":
                continue
            self._base_clade_color(tree, clade)

    def __call__(self, output: str, dpi: int = 300, width: int = 15000):
        self.tree.render(output, dpi=dpi, w=width, tree_style=tree_style())


def _parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Generates hcov-19 Phylogentic tree based on GISAID bi-weekly analysis update",
        epilog="If you notice any issues, please open one over at https://github.com/ElasticBottle/GISAID-Analysis-Update ",
    )
    parser.add_argument(
        "file", type=str, default="", help="file path of the NEWICK tree to parse",
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
        "-d",
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
        help="Specifies the width of the image in pixel, the height is then automatically derived. Default 15000",
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
    start = time.time()
    img = PhyloIMG(file=args.file)
    img(
        output=args.out, dpi=args.dpi, width=args.width,
    )
    end = time.time()
    print(f"Time taken: {end - start:.2f}s")


if __name__ == "__main__":
    main()

# "./data/input_sample/fnb_all_rapidnj.nwk"
# "phylo_tree_0width.png"


# t = Tree(newick="./test/samples/phylo/fnb_all_rapidnj.nwk")
# root = t.get_tree_root()
# test = list(root.iter_leaves())[0]
# test.up
