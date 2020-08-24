from ete3 import Tree, TreeStyle, NodeStyle
from typing import Union


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
        "G_": "#fa0000",
        "GH": "#ff9933",
        "GR": "#ff3333",
        "S_": "LightGreen",
        "V_": "Violet",
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
    return style
    # style["vt_line_width"] = 2
    # style["hz_line_width"] = 2


def get_clade_color(name: str):
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


def get_tree(file: str, test: bool = False) -> Tree:
    tree = Tree(newick=file)
    for i, n in enumerate(tree.iter_leaves()):
        clade, color = get_clade_color(n.name)
        n.img_style = node_style(clade, color, n.is_leaf())
        if test and i == 100:
            break
    return tree


tree = get_tree("./data/input_sample/fnb_all_rapidnj.nwk")
tree.render("phylo_tree.png", dpi=300, w=10000, tree_style=tree_style())

