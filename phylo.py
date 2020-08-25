from ete3 import Tree, TreeNode, TreeStyle, NodeStyle
from typing import Union, Tuple
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


def node_style(clade: str, color: str, is_leaf: bool,) -> NodeStyle:
    color_scheme = {
        "G_": "#fdb3b3",
        "GH": "#feab8d",
        "GR": "#fe8d8d",
        "S_": "#f0b3f0",
        "V_": "#b3ffb3",
        "L_": "White",
        "O_": "White",
    }
    style = NodeStyle()
    style["fgcolor"] = "#0f0f0f"
    style["size"] = 0
    style["shape"] = "circle"
    style["vt_line_color"] = color if color is not None else "#000000"
    style["hz_line_color"] = color if color is not None else "#000000"
    style["vt_line_type"] = 0  # 0 solid, 1 dashed, 2 dotted
    style["hz_line_type"] = 0
    style["bgcolor"] = color_scheme.get(clade, "White")
    style["vt_line_width"] = 20
    style["hz_line_width"] = 20
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
        self.tree = self.get_tree(self.file)

    def _update_clades(self, clade: str, node: TreeNode):
        """
        Update [self.clades] with the latest [clade] of the taxonomy
        
        Args:
            - clade(str): The clade of the taxonomy
            - node(TreeNode): The node to be updated
        """
        clade_list = self.clades.get(clade, self.clades["O_"])
        clade_list.append(node.name)
        if clade in self.clades_list:
            self.clades[clade] = clade_list
        else:
            self.clades["O_"] = clade_list

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
        s_ancestor = tree.get_common_ancestor(*self.clades["S_"])

        tree.set_outgroup(s_ancestor)

    def _color_ancestors(self, tree: TreeNode):
        # s_ancestor.img_style = node_style("S_", None, False)
        for clade in self.clades_list:
            ancestor = tree.get_common_ancestor(*self.clades.get(clade, []))
            ancestor.img_style = node_style(clade, None, False)

    def get_tree(self, file: str, test: bool = False) -> Tree:
        """

        """
        tree = Tree(newick=file)

        for i, n in enumerate(tree.iter_leaves()):
            self._color_node(n)
            if test and i == 1000:
                break
        self._root_on_s(tree)
        # self._color_ancestors(tree)

        return tree

    def __call__(self, output: str, dpi: int, width: int):
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
        metavar="",
        default="./out.png",
        help="The output path of the image to be stored. Default ./out.png",
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
    img(
        output=args.out, dpi=args.dpi, width=args.width,
    )


if __name__ == "__main__":
    main()
# "./data/input_sample/fnb_all_rapidnj.nwk"
# "phylo_tree_0width.png"

