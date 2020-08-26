from ete3 import Tree, TreeNode, TreeStyle, NodeStyle
from typing import Union, Tuple, Dict, List, Iterable
from collections import defaultdict
import argparse


def tree_style() -> TreeStyle:
    style = TreeStyle()
    style.mode = "c"  # draw tree in circular mode
    # style.scale = 20
    # style.branch_vertical_margin = 15
    style.show_leaf_name = False
    style.show_branch_length = False
    style.show_branch_support = False
    return style


def node_style(
    clade: str, color: str, is_leaf: bool, bg_color: str = None,
) -> NodeStyle:

    style = NodeStyle()
    style["fgcolor"] = "#0f0f0f"
    style["size"] = 0
    style["shape"] = "circle"
    style["vt_line_color"] = color if color is not None else "#000000"
    style["hz_line_color"] = color if color is not None else "#000000"
    style["vt_line_type"] = 0  # 0 solid, 1 dashed, 2 dotted
    style["hz_line_type"] = 0
    style["vt_line_width"] = 20
    style["hz_line_width"] = 20
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
        self.clades_list = [
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
        }
        self.leaf_clade_group = {}  # taxon_name -> clade
        self.tree = self.get_tree(self.file, test=False)

    def _update_clades(self, clade: str, node: TreeNode):
        """
        Update [self.clades] and [self.leaf_clade_group] with the taxon and clade respectively
        
        Args:
            - clade(str): The clade of the taxonomy
            - node(TreeNode): The node to be updated
        """
        clade_list = self.clades.get(clade, self.clades["O_"])
        clade_list.append(node.name)
        if clade in self.clades_list:
            self.clades[clade] = clade_list
            self.leaf_clade_group[node.name] = clade
        else:
            self.clades["O_"] = clade_list
            self.leaf_clade_group[node.name] = "O_"

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

    def _color_node(self, n: TreeNode):
        """
        Sets the NodeStyle for [n] and updates the clades dictionary
        """
        clade, color = self._get_clade_color(n.name)
        n.img_style = node_style(clade, color, n.is_leaf())
        self._update_clades(clade, n)

    def _root_on_s(self, tree: TreeNode):
        """
        Roots the tree base on clade S
        """
        s_ancestor = tree.get_common_ancestor(*self.clades["S_"])
        tree.set_outgroup(s_ancestor)

    def _get_majority_clade(self, children: Iterable) -> str:
        """
        Finds the clade group occuring in the largest quantity from a particular node

        Args:
            -children (Iterable): The leaf from the particular node of interest to search

        Returns:
            str: the clade that occurs most frequently
        """
        count = defaultdict(int)
        for child in children:
            if not child.is_leaf():
                continue
            count[self.leaf_clade_group[child.name]] += 1
        return max(count, key=count.get)

    def _color_clade(self, n: TreeNode):
        """
        Sets the ['bgcolor'] of a node's [img_style] to be the appropriate color based on which clade group appears most frequently

        Args:
            n(TreeNode) - The node whose background color is to be set
        """
        if n.is_leaf() or n.is_root():
            return
        clade = self._get_majority_clade(n.iter_leaves())
        if clade == "G_" or clade == "GH" or clade == "GR":
            ancestors = n.get_ancestors()
            if len(ancestors) > 2:
                clade = self._get_majority_clade(n.up.up.iter_leaves())
        n.img_style = node_bg_style(n.img_style, clade)

    def get_tree(self, file: str, test: bool = False) -> Tree:
        """
        
        """
        tree = Tree(newick=file)

        for i, n in enumerate(tree.iter_leaves()):
            self._color_node(n)
        for i, n in enumerate(tree.traverse()):
            self._color_clade(n)
            if test and i == 5:
                break
        self._root_on_s(tree)
        tree.write(format=3, outfile="test.nwk")
        return tree

    def __call__(self, output: str, dpi: int = 300, width: int = 15000):
        ext = output[-4:]
        self.tree.write(format=3, outfile=output.replace(ext, ".nwk"))
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
        nargs=1,
        default=300,
        type=int,
        dest="dpi",
        help="Dpi of the output image. Default 300",
    )
    parser.add_argument(
        "-w",
        "--width",
        nargs=1,
        default=15000,
        dest="width",
        help="Specifies the width of the image in pixel, the height is then automatically derived. Default 15000",
    )
    args = parser.parse_args()
    return args


def main():
    args = _parse_args()
    img = PhyloIMG(file=args.file,)
    print("generating image, this will take a while. ~2-5 mins")
    img(
        output=args.out, dpi=args.dpi, width=args.width,
    )


if __name__ == "__main__":
    main()

# "./data/input_sample/fnb_all_rapidnj.nwk"
# "phylo_tree_0width.png"

