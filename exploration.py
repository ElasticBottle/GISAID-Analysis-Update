#%%
from typing import Tuple, Union, List, Dict

import numpy as np
import pandas as pd

from ete3 import Tree, TreeNode

file = "./test/samples/input.nwk"
tree = Tree(newick=file, quoted_node_names=True)
clades_list = [
    "G_",
    "GH",
    "GR",
    "S_",
    "V_",
    "L_",
    "O_",
]
clades = {
    "G_": [],
    "GH": [],
    "GR": [],
    "S_": [],
    "V_": [],
    "L_": [],
    "O_": [],
}
leaf_clade_group = {}


def get_clade_color(name: str) -> Tuple[str, Union[str, None]]:
    """
        Args:
            - name(str): In one of the following formats:
                GR_USAMNMDH1413_20200623_NorthAmerica_[&!color=ff00cc]
                GR_EnglandLONDD4621_20200408_Europe_
        Returns:
            - str: the clade that the taxonomy belongs too
            - str / None : the color for that clade if any
        """
    name_color = name.split("color=")
    clade = name_color[0][0:2]
    color = None
    if len(name_color) == 2:
        color = name_color[1][:-1]
    return (clade, color)


def update_clades(clade: str, node: TreeNode):
    """
    Update clades and leaf_clade_group with the taxon and clade respectively
    
    Args:
        - clade(str): The clade of the taxonomy
        - node(TreeNode): The node to be updated
    """
    clade_list = clades.get(clade, clades["O_"])
    clade_list.append(node)
    if clade in clades_list:
        clades[clade] = clade_list
        leaf_clade_group[node.name] = clade
    else:
        clades["O_"] = clade_list
        leaf_clade_group[node.name] = "O_"


for leaf in tree.iter_leaves():
    if leaf.name == "WuhanWIV04_20191230_China_":
        root = leaf
    clade, color = get_clade_color(leaf.name)
    update_clades(clade, leaf)

#%%
import matplotlib.pyplot as plt


def mean_std(clades: Dict[str, List]) -> Dict[str, List]:
    result = {}  # clade -> mean_std
    for k, v in clades.items():
        dist = []
        for node in v:
            dist.append(node.get_distance(root))

        result[k] = dist
    return result


# details = mean_std(clades)

# plt.boxplot(list(details.values()), labels=list(details.keys()), whis=[0, 100])
# plt.show()

#%%
children = []
for node in tree.traverse():
    children.append(len(list(node.iter_leaves())))
plt.boxplot(children, showfliers=False)
plt.show()
#%%
# how many leaves are there at each branch level?
count = {}
for child in children:
    count[(child)] = count.get((child), 0) + 1
for k in sorted(count.keys()):
    print(k, count[k])
len(count)
# %%
# How many sister does each level have?
for i, node in enumerate(tree.iter_leaves()):
    sisters = node.get_sisters()
    sisters.append(node)
    clades_sis = list(map(lambda x: get_clade_color(x.name)[0], sisters))
    if len(pd.Series(clades_sis)) > 2:
        print(pd.Series(clades_sis))

#%%
# for each clade ancestor, how much noise is there?
for clade, clade_list in clades.items():
    print(clade)
    ancestor = tree.get_common_ancestor(*clade_list)
    noise = {}
    for leafs in ancestor.iter_leaves():
        c, _ = get_clade_color(leafs.name)
        if c in clades_list:
            noise[c] = noise.get(c, 0) + 1
        else:
            noise["O_"] = noise.get("O_", 0) + 1
    print(noise)

#%%
# Is there an ancestor that can maximise for the number of taxons belonging to a clade?

from typing import Dict
import numpy as np


def find_num_nodes(tree: TreeNode) -> Dict[str, int]:
    """
    Finds the total number of node per clade
    """
    num_clade = {}
    for leaf in tree.iter_leaf_names():
        clade, _ = get_clade_color(leaf)
        num_clade[clade] = num_clade.get(clade, 0) + 1
    return num_clade


def calculate_density(node: TreeNode, clade: str) -> float:
    """
    Calculate the percentage of [clade] for a given [node]
    """
    node_per_clade = find_num_nodes(node)
    node_count = np.sum(list(node_per_clade.values()))
    print(
        node_per_clade, node_count, node_per_clade[clade] / node_count * 100,
    )
    return node_per_clade[clade] / node_count * 100


def max_ancestor(tree: TreeNode, clade: str, clade_total: int):
    """
    For each unique clade
        traverse the tree top find the internal node with the highest proportion of said clade, provided that it encompasses more than half of the said clade groups.
    Example:
    {'GH': 947, 'O_': 304, 'G_': 887, 'S_': 186, 'GR': 1402, 'L_': 77, 'V_': 108}

    For GH, find the internal node where:
        - 473.5 < GH < 947, and the proportion of GH to other clade is the highest
    """
    min_size = clade_total // 2
    node_density = {}
    for leaf in tree.traverse():
        if leaf.is_leaf() or leaf.is_root():
            continue
        if find_num_nodes(leaf).get(clade, 0) >= min_size:
            density = calculate_density(leaf, clade)
            node_density[leaf] = density
    sorted_nodes = sorted(node_density.keys(), key=node_density.get)
    max_node = sorted_nodes[-1]
    return (max_node, node_density[max_node])


nodes_per_clade = find_num_nodes(tree)
for clade in clades_list:
    print(clade)
    node, density = max_ancestor(tree, clade, nodes_per_clade[clade])
    print((node, density))
