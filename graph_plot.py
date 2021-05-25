import argparse
from collections import OrderedDict
from typing import Dict, List

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PolyCollection
import numpy as np
import pandas as pd
import seaborn as sns


def remove_borders():
    _, ax = plt.subplots()
    right_side = ax.spines["right"]
    right_side.set_visible(False)
    top_side = ax.spines["top"]
    top_side.set_visible(False)


def plot_3d(index: pd.Index, df: pd.DataFrame, out: str):
    def make_tick_index(index, scale):
        result = []
        for idx in range(len(index) * scale):
            if (idx % scale) == 0:
                result.append(index[int(idx / scale)])
            else:
                result.append("")
        print(result)
        return result

    pallette = OrderedDict(
        [
            ("L", "#d9d9d9"),
            ("O", "#7f7f7f"),
            ("S", "#c6e0b4"),
            ("V", "#ffccff"),
            ("G", "#ffcccc"),
            ("G+RBDx", "#ffcccc"),
            ("GH", "#f5b183"),
            ("GH+RBDx", "#f5b183"),
            ("GR", "#fe7c80"),
            ("GR+RBDx", "#fe7c80"),
            ("GRY", "#fed88d"),
            ("GRY+RBDx", "#fed88d"),
        ]
    )
    to_hatch = ["G+RBDx", "GH+RBDx", "GR+RBDx", "GRY+RBDx"]
    non_g_clades = [
        "L",
        "O",
        "S",
        "V",
    ]
    non_g_clade_len = len(non_g_clades)
    g_clades = df.columns.difference(non_g_clades).tolist()
    g_clade_len = len(g_clades)
    order = {
        "L": 4,
        "O": 5,
        "S": 6,
        "V": 7,
        "G, G+RBDx": 3,
        "GH, GH+RBDx": 2,
        "GRY, GRY+RBDx": 1,
        "GR, GR+RBDx": 0,
    }
    # pallette = OrderedDict(
    #     [
    #         ("S", "#c6e0b4"),
    #         ("L", "#d9d9d9"),
    #         ("V", "#ffccff"),
    #         ("G", "#ffcccc"),
    #         ("GH", "#f5b183"),
    #         ("GR", "#fe7c80"),
    #         ("O", "#7f7f7f"),
    #         # ("G+RBDx", "#ffcccc"),
    #         # ("GH+RBDx", "#f5b183"),
    #         # ("GR+RBDx", "#fe7c80"),
    #     ]
    # )
    # order = {
    #     "S": 0,
    #     "L": 1,
    #     "V": 2,
    #     "G": 3,
    #     "GH": 4,
    #     "GR": 5,
    #     "O": 6,
    # }

    fig = plt.figure()
    ax = plt.gca(projection="3d")
    xs = np.arange(0, len(index), 1)
    max_count = 6200

    # Creating the polygons to plot
    vert = [None] * non_g_clade_len
    for idx, series in df.items():
        if idx in non_g_clades:
            ys = series.tolist()
            vert[order[idx] - 3] = list(zip(list(reversed(xs)), ys))
    non_g_pallette = [x for x in pallette.items() if x[0] in non_g_clades]
    colors = sorted(non_g_pallette, key=lambda x: order.get(x[0]))
    colors = dict((x, y) for x, y in colors)
    non_g_poly = PolyCollection(vert, facecolors=colors.values())
    ax.add_collection3d(non_g_poly, zs=range(3, non_g_clade_len + 3), zdir="y")

    g_poly = []
    # Plotting G polygons
    if g_clade_len == 6:
        for i in range(0, g_clade_len - 1, 2):
            value = [df[g_clades[i]].tolist(), df[g_clades[i + 1]].tolist()]
            stacks = ax.stackplot(
                list(reversed(xs)),
                *value,
                colors=[x[1] for x in pallette.items() if x[0] == g_clades[i]],
            )
            for stack, key in zip(stacks, [g_clades[i], g_clades[i + 1]]):
                stack.set_alpha(0.8)
                ax.add_collection3d(stack, zs=2 - (i / 2), zdir="y")

                if key in to_hatch:
                    stack.set_hatch("++")

            g_poly.append(stacks)
    else:
        raise NotImplementedError("Only G Gh GR not implemented")

    # setting the x axis and labels
    ax.set_xlabel("Date", labelpad=12)
    ax.set_xticks(list(range(len(index))))
    ax.set_xticklabels(reversed(index), rotation=45, ha="right")

    # setting the y axis and labels
    # ax.set_ylabel("Clade", labelpad=2)
    ax.set_yticks([-1] + list(range((len(order)))))
    ax.set_yticklabels(
        [""] + list(sorted(order.keys(), key=lambda x: order.get(x))),
        ha="left",
        rotation=-20,
    )

    # setting the z axis and label
    ax.zaxis.set_rotate_label(False)
    ax.set_zlabel("Count", rotation=0, labelpad=7)
    ax.set_zlabel("Count")
    ax.set_zlim3d(0, max_count)

    # adjusting the tick distance from axis
    ax.tick_params(axis="x", which="major", pad=-5)
    ax.tick_params(axis="y", which="major", pad=-2)

    # make the panes transparent
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

    # make the grid lines transparent
    # ax.xaxis._axinfo["grid"]["color"] = (1, 1, 1, 0)
    ax.yaxis._axinfo["grid"]["color"] = (1, 1, 1, 0)
    # # ax.zaxis._axinfo["grid"]["color"] = (1, 1, 1, 0)

    # setting the view point
    ax.view_init(elev=15.0, azim=135)
    # plt.show()
    plt.savefig(out, dpi=320)
    return fig, ax


def plot_stacked_area(index: pd.Index, labels: List, values: List, out: str):
    pallette = OrderedDict(
        [
            ("O", "#7f7f7f"),
            ("S", "#c6e0b4"),
            ("L", "#d9d9d9"),
            ("V", "#ffccff"),
            ("G", "#ffcccc"),
            ("G+RBDx", "#ffcccc"),
            ("GV", "#f08bb5"),
            ("GV+RBDx", "#f08bb5"),
            ("GH", "#f5b183"),
            ("GH+RBDx", "#f5b183"),
            ("GR", "#fe7c80"),
            ("GR+RBDx", "#fe7c80"),
            ("GRY", "#fed88d"),
            ("GRY+RBDx", "#fed88d"),
        ]
    )
    to_hatch = ["G+RBDx", "GH+RBDx", "GR+RBDx", "GV+RBDx", "GRY+RBDx"]

    fig, ax = plt.subplots()
    stacks = ax.stackplot(
        index, *values, labels=labels, colors=list(pallette.values()),
    )

    for stack, key in zip(stacks, pallette.keys()):
        if key in to_hatch:
            stack.set_hatch("++")

    y_ticks = list(map(lambda x: f"{x:.0f}%", np.arange(110, step=10)))
    plt.yticks(
        ticks=np.arange(110, step=10), labels=y_ticks,
    )
    plt.xticks(rotation=90)
    ax.margins(y=0)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + 0.2, box.width, box.height * 0.85])
    ax.legend(ncol=5, loc="lower center", 
              # Adjusted in April/May 2021.
              bbox_to_anchor=(0.6, -0.45))
    
    for item in [fig, ax]:
        item.patch.set_visible(False)
        
    # Added in April/May 2021.
    ## Add another x-axis!
    sns.despine(ax=ax, top=True, bottom=False, 
                right=False, trim=True,
                offset={"right": 75, "bottom":15})
    #ax.axvline(ax.get_xlim()[1]*0.9999, color='k',lw=2)
    ax.tick_params(labelbottom=True, labeltop=False, 
                    labelleft=True, labelright=True,
                    bottom=True, top=False, 
                    left=True, right=True)
                    
    ## Annotate.
    anno_palette = {
        "G": "#ffcccc",
        "GR": "#ff7c80",
        "GH": "#f4b183",
        "GRY": "#fed88d"
    }
    
    anno_labels = {
        "G": "Includes B.1.617+\nand B.1.525",
        "GR": "Includes P.1",
        "GH": "Includes B.1.351\nand B.1.429",
        "GRY": "GRY = B.1.1.7"
    }
    
    anno_positions = {
        "GRY": 0.6,
        "GR": 0.23,
        "GH": 0.12,
        "G": 0.025,
    }
    
    for clade in anno_labels.keys():
        ax.annotate(
            anno_labels[clade],
            xy=(0.97, anno_positions[clade]), xycoords='axes fraction',
            xytext=(0, 0), textcoords='offset points',
            fontsize=7,
            bbox=dict(boxstyle="round", 
                      fc=anno_palette[clade],
                      ec=anno_palette[clade])
        )
        
    ax.annotate('Small and/or Variable\nSample Sizes',
        fontsize=8, color='#8da9db',
        xy=(0.78, 0.22), xycoords='figure fraction',
        xytext=(0.82, 0.22), textcoords='figure fraction',
        arrowprops=dict(facecolor='#8da9db', 
                        edgecolor='#8da9db',
                        width=0.5),
        horizontalalignment='left', verticalalignment='center')
        
    fig.savefig(out, dpi=300, bbox_inches='tight')
    
    
    return fig, ax


def generate_clade_progression(file: str, out: str, is3d: bool):
    """
    Generates a Stacked area plot for the clade progression tsv file at [file] and stores the result in [out]

    Args:
        - file(str): The input file for the clade progression. Must be in .tsv format
            Sample format:
                Per month:

                Month	G	GH	GR	L	O	S	V
                Dec19	0	0	0	17	1	0	0
                Jan20	14	1	0	230	78	167	5
                Feb20	91	36	62	277	265	224	147
                Mar20	6390	6897	5844	2185	1361	3346	3015
                Apr20	7368	6276	9736	977	1038	1048	1286
                May20	2480	2628	4454	107	496	214	90
                Jun20	887	1210	2784	18	201	96	26
                Jul20	490	360	1027	3	72	1	5
                Aug20	6	17	67	0	1	1	0

                Cumulative:

                Month	G	GH	GR	L	O	S	V
                Dec19	0	0	0	17	1	0	0
                Jan20	14	1	0	247	79	167	5
                Feb20	105	37	62	524	344	391	152
                Mar20	6495	6934	5906	2709	1705	3737	3167
                Apr20	13863	13210	15642	3686	2743	4785	4453
                May20	16343	15838	20096	3793	3239	4999	4543
                Jun20	17230	17048	22880	3811	3440	5095	4569
                Jul20	17720	17408	23907	3814	3512	5096	4574
                Aug20	17726	17425	23974	3814	3513	5097	4574
        -out(str): the output path for the generated graph
    """
    df = pd.read_csv(file, sep="\t", skiprows=1, index_col=0,)
    insert_order = {
        "Other": 0,
        "S": 1,
        "L": 2,
        "V": 3,
        "G": 4,
        "G+RBDx": 5,
        "GV": 6,
        "GV+RBDx": 7,
        "GH": 8,
        "GH+RBDx": 9,
        "GR": 10,
        "GR+RBDx":11,
        "GRY": 12,
        "GRY+RBDx":13,
    }

    labels = OrderedDict(
        [
            ("G", "G"),
            ("Gn", "G+RBDx"),
            ("GH", "GH"),
            ("GHn", "GH+RBDx"),
            ("GR", "GR"),
            ("GRn", "GR+RBDx"),
            ("GRY", "GRY"),
            ("GRYn", "GRY+RBDx"),
            ("GV", "GV"),
            ("GVn", "GV+RBDx"),
            ("Other", "O"),
            ("S", "S"),
            ("L", "L"),
            ("V", "V"),
        ]
    )

    # Retreive the monthly dataframe
    monthly_df = df.loc[:"Cumulative:", :].iloc[:-1]
    monthly_df = monthly_df.astype(float)
    # If ever need to retrieve the cumulative_dataFrame
    # cum_df = df.loc["Month":, :].iloc[1:]
    # cum_df = cum_df.astype(float)

    monthly_df.rename(columns=labels, inplace=True)
    hor_sum = monthly_df.sum(axis=1)

    # Find the percentage of the the clades // mutation
    clade_df_percentage = monthly_df.div(hor_sum, axis=0) * 100

    # Separate out the individual rows for plotting
    values = [None] * len(insert_order)

    # Making the values for the clades
    for label, value in clade_df_percentage.items():
        values[insert_order.get(label, 0)] = value

    remove_borders()
    fig, axes = (
        plot_3d(monthly_df.index, monthly_df, out,)
        if is3d
        else plot_stacked_area(
            clade_df_percentage.index,
            sorted(labels.values(), key=lambda x: insert_order.get(x, 0)),
            values,
            out,
        )
    )



def make_df(file: str, continents: List[str]) -> Dict[str, pd.Series]:
    """
    Builds the dataFrame given a file

    Args:
        - file(str): the input path to the tsv file containing the relevant info.
            Example formatting:
                Europe (+104)
                G	47
                GH	8
                GR	49
                L	0
                O	0
                S	0
                V	0

                Asia (+40)
                G	8
                GH	0
                GR	2
                L	0
                O	26
                S	4
                V	0

                Africa (+315)
                G	48
                GH	0
                GR	267
                L	0
                O	0
                S	0
                V	0

                Oceania (+2)
                G	0
                GH	2
                GR	0
                L	0
                O	0
                S	0
                V	0

                NorthAmerica (+484)
                G	112
                GH	311
                GR	44
                L	7
                O	1
                S	9
                V	0

                SouthAmerica (+26)
                G	12
                GH	3
                GR	10
                L	0
                O	0
                S	1
                V	0

        - continents(str): a list of continents that dataframes needs to be build for (6 in total)
    """

    def get_clade_value(f):
        clades, values = [], []
        for _ in range(10):
            clade, value = f.readline().split("\t")
            assert clade is not None and value is not None
            assert len(clade) > 0 and len(value) > 0
            value = value[:-1]
            clades.append(clade)
            values.append(value)
        return (clades, values)

    dfs = OrderedDict()  # continent -> dataFrame of clades
    with open(file, "r") as f:
        for line in f:
            words = line.split("\t")
            if any(continent in words[0] for continent in continents):
                label = words[0][:-1]
                clades, values = get_clade_value(f)

                dfs[label] = pd.Series(values, index=clades).astype(int)
    assert len(dfs) == 6
    return dfs


def make_pie_chart(index: int, values: pd.Series, title: str, axs):
    """
    Builds a pie chart at [index] of a 2, 3 subplot with [values] and [title]

    Args:
        - index(int): the location of the pie chart. Increases right to left and up to down
        - values(pd.Series): The pandas series contain both the labels of the individual portion of the chart and the actual values for each clade
        - title(str): the title of the chart
        - axs: The maplotlib Axs to plot the pie chart on
    """

    # def get_index(index: int) -> List[int]:
    #     index -= 1
    #     if index >= 6:
    #         index -= 6
    #         return [1, index]
    #     return [0, index]

    labels = OrderedDict(
        [
            ("Other", "O"),
            ("S", "S"),
            ("L", "L"),
            ("V", "V"),
            ("G", "G"),
            ("Gn", "G+RBDx"),
            ("GH", "GH"),
            ("GHn", "GH+RBDx"),
            ("GR", "GR"),
            ("GRn", "GR+RBDx"),
            ("GRY", "GRY"),
            ("GRYn", "GRY+RBDx"),
            ("GV", "GV"),
            ("GVn", "GV+RBDx"),
        ]
    )
    # Must correspond to order in data
    colors = {
        "G": "#ffcccc",
        "G+RBDx": "#ffcccc",
        "GH": "#f4b183",
        "GH+RBDx": "#f4b183",
        "GR": "#ff7c80",
        "GR+RBDx": "#ff7c80",
        "GRY": "#fed88d",
        "GRY+RBDx": "#fed88d",
        "GV": "#f08bb5",
        "GV+RBDx": "#f08bb5",
        "L": "#d9d9d9",
        "O": "#808080",
        "S": "#70ad47",
        "V": "#ff99ff",
    }
    to_hatch = ["G+RBDx", "GH+RBDx", "GR+RBDx", "GRY+RBDx", "GV+RBDx"]
    #indices = get_index(index)
    # print(title, values)
    values.rename(index=labels, inplace=True)
    # print("renamed", values)
    normalize = values.values.sum(axis=0) != 0

    # axs[indices[0], indices[1]].set_title(label=title, fontdict={"fontsize": 10})
    # pie = axs[indices[0], indices[1]].pie(
    #     values.values, normalize=normalize, labeldistance=None, colors=colors.values()
    # )
    
    axs[index-1].set_title(label=title.replace(" (+", "\n(+ "), fontdict={"fontsize": 9})
    pie = axs[index-1].pie(
        values.values, normalize=normalize, labeldistance=None, colors=colors.values()
    )
    
    plotted_val = values[values > 0]
    plotted_clades = sorted(
        plotted_val.index, key=lambda x: plotted_val[x], reverse=True
    )

    wedges = list(map(lambda x: (x, x.theta2 - x.theta1), pie[0]))
    wedges = sorted(wedges, key=lambda x: x[1], reverse=True)
    for wedge, clade in zip(wedges, plotted_clades):
        # print(wedge, clade)
        # print(f"{clade} in {to_hatch} = {clade in to_hatch}")
        if clade in to_hatch:
            wedge[0].set_hatch("+++")

    wedges = {}
    for key, value in colors.items():
        patch = mpatches.Patch(
            facecolor=value, hatch="+++" if key in to_hatch else "", label=key
        )
        wedges[key] = patch
    return wedges


def generate_geoclade_progression(file: str, out: str):
    """
    Generates a pie chart showing the number of infection per continent separated based on clades

    Args:
        - file(str): the input path to the tsv file containing the relevant info
        - out(str): the output path for the chart to be saved in
    """
    continents = [
        "Europe",
        "Asia",
        "Africa",
        "Oceania",
        "NorthAmerica",
        "SouthAmerica",
    ]
    #continents_index = [2, 3, 5, 6, 1, 4]
    continents_index = [3, 2, 1, 4, 5, 6]
    dfs = make_df(file, continents)
    fig, axs = plt.subplots(1, 6, figsize=(10,4),
                            #gridspec_kw={'wspace':1.75}
                        )
    patches = {} 
    for i, item in enumerate(dfs.items()):
        title, chart_labels = item
        graph_index = continents_index[i]
        patches.update(make_pie_chart(graph_index, chart_labels, title, axs))
    legend_ax = axs[1]
    fig.legend(
        handles=list(patches.values()),
        labels=list(patches.keys()),
        ncol=6,
        loc="upper left",
        #bbox_to_anchor=(0.5, 0.05),
        bbox_to_anchor=(0.1, 0.05, 0.8, 0.25),
        mode='expand',
        borderaxespad=0
    )
    
    # for i in np.arange(1,6):
    #     axs[1,i].set_visible(False)
        
    fig.savefig(out, dpi=300, bbox_inches='tight')


def _parse_args():
    """ handles parsing of cli arguments"""
    parser = argparse.ArgumentParser(
        description="Generates the clade progression or clade graph depending specification",
        epilog="If you notice any issues, please open one over at https://github.com/ElasticBottle/GISAID-Analysis-Update ",
    )
    parser.add_argument(
        "c_file",
        metavar="clade_prog_file",
        type=str,
        help="""Input file for the clade progression. File format should be .tsv
        Pass a random string if setting -genc False""",
    )
    parser.add_argument(
        "g_file",
        metavar="geoclade_file",
        type=str,
        help="""Input file for the geo clade progression. File format should be .tsv
        Pass a random string if setting -geng False""",
    ),
    parser.add_argument(
        "--c_out",
        "-co",
        type=str,
        metavar="",
        dest="c_out",
        default="./clade_prog.png",
        help="The output path of the clade progression graph. Must include file format. E.g. './clade_prog.png'",
    )
    parser.add_argument(
        "--g_out",
        "-go",
        type=str,
        metavar="",
        dest="g_out",
        default="./geoclade.png",
        help="The output path of the GeoClade progression graph. Must include file format. E.g. './outg.png'",
    )
    parser.add_argument(
        "--genc",
        "-gc",
        dest="gc",
        action="store_false",
        help="Flag to disbale generation of the clade progression graph. If enabled, specify gibberish string for clade_prog_file input",
    )
    parser.add_argument(
        "--geng",
        "-gg",
        dest="gg",
        action="store_false",
        help="Flag to disable the generation of the geoclade graph. If enabled, specify gibberish string for geo_clade_file input",
    )

    parser.add_argument(
        "-c3d",
        "--clade3d",
        dest="c3d",
        action="store_true",
        help="Uses 3d stack area plot instead of regular 3d plot",
    )
    return parser.parse_args()


def main():
    args = _parse_args()
    if args.gc:
        generate_clade_progression(args.c_file, args.c_out, args.c3d)
    if args.gg:
        generate_geoclade_progression(args.g_file, args.g_out)


if __name__ == "__main__":
    main()
# r"./data/input_sample/clade_progression.tsv"
