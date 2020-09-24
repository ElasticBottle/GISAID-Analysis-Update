import argparse
from collections import OrderedDict
from typing import Dict, List

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns


def remove_borders():
    _, ax = plt.subplots()
    right_side = ax.spines["right"]
    right_side.set_visible(False)
    top_side = ax.spines["top"]
    top_side.set_visible(False)


def plot_stacked_area(index: pd.Index, labels: List, values: List, out: str):
    pallette = OrderedDict(
        [
            ("O", "#7f7f7f"),
            ("S", "#c6e0b4"),
            ("L", "#d9d9d9"),
            ("V", "#ffccff"),
            ("G", "#ffcccc"),
            ("G+S477X", "#ffcccc"),
            ("GH", "#f5b183"),
            ("GH+S477X", "#f5b183"),
            ("GR", "#fe7c80"),
            ("GR+S477X", "#fe7c80"),
        ]
    )
    to_hatch = ["G+S477X", "GH+S477X", "GR+S477X"]

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
    ax.legend(ncol=5, loc="lower center", bbox_to_anchor=(0.5, -0.4))
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    return fig, ax


def plot_overlap_line_graph(fig, axes, x, y):
    print("non clade stuff", y)
    pallette = [
        "#f2a210",
        "#f53312",
        "#117283",
    ]
    for index, series_tup in enumerate(y.items()):
        label, value = series_tup
        axes.plot(x, value, label=label, color=pallette[index])
    y_ticks = list(map(lambda x: f"{x:.0f}%", np.arange(110, step=10)))
    plt.yticks(
        ticks=np.arange(110, step=10), labels=y_ticks,
    )
    axes.legend(ncol=3, loc="lower center", bbox_to_anchor=(0.5, -0.5))
    # Make some space on the right side for the extra y-axis.
    fig.subplots_adjust(right=0.9)
    fig.savefig("./updated_curation/output/random.png", dpi=300)


def generate_clade_progression(file: str, out: str):
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
        "O": 0,
        "S": 1,
        "L": 2,
        "V": 3,
        "G": 4,
        "G+S477X": 5,
        "GH": 6,
        "GH+S477X": 7,
        "GR": 8,
        "GR+S477X": 9,
    }

    labels = OrderedDict(
        [
            ("O", "O"),
            ("S", "S"),
            ("L", "L"),
            ("V", "V"),
            ("G", "G"),
            ("Gn", "G+S477X"),
            ("GH", "GH"),
            ("GHn", "GH+S477X"),
            ("GR", "GR"),
            ("GRn", "GR+S477X"),
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
    fig, axes = plot_stacked_area(
        clade_df_percentage.index,
        sorted(labels.values(), key=lambda x: insert_order.get(x, 0)),
        values,
        out,
    )
    # plot_overlap_line_graph(
    #     fig, axes.twinx(), mutation_df_percentage.index, mutation_df_percentage
    # )


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

    def get_index(index: int) -> List[int]:
        index -= 1
        if index >= 3:
            index -= 3
            return [1, index]
        return [0, index]

    labels = OrderedDict(
        [
            ("O", "O"),
            ("S", "S"),
            ("L", "L"),
            ("V", "V"),
            ("G", "G"),
            ("Gn", "G+S477X"),
            ("GH", "GH"),
            ("GHn", "GH+S477X"),
            ("GR", "GR"),
            ("GRn", "GR+S477X"),
        ]
    )
    colors = {
        "G": "#ffcccc",
        "G+S477X": "#ffcccc",
        "GH": "#f4b183",
        "GH+S477X": "#f4b183",
        "GR": "#ff7c80",
        "GR+S477X": "#ff7c80",
        "L": "#d9d9d9",
        "O": "#808080",
        "S": "#70ad47",
        "V": "#ff99ff",
    }
    to_hatch = ["G+S477X", "GH+S477X", "GR+S477X"]
    indices = get_index(index)
    # print(values)
    values.rename(index=labels, inplace=True)
    # print("renamed", values)
    normalize = values.values.sum(axis=0) != 0
    # print(normalize)

    axs[indices[0], indices[1]].set_title(label=title, fontdict={"fontsize": 10})
    pie = axs[indices[0], indices[1]].pie(
        values.values, normalize=normalize, labeldistance=None, colors=colors.values()
    )
    plotted_val = values[values > 0]

    for wedge, clade in zip(pie[0], plotted_val.index):
        if clade in to_hatch:
            wedge.set_hatch("++")

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
    continents_index = [2, 3, 5, 6, 1, 4]
    dfs = make_df(file, continents)
    fig, axs = plt.subplots(2, 3)
    patches = {}
    for i, item in enumerate(dfs.items()):
        title, chart_labels = item
        graph_index = continents_index[i]
        patches.update(make_pie_chart(graph_index, chart_labels, title, axs))
    axs[1, 0].legend(
        handles=list(patches.values()),
        labels=list(patches.keys()),
        ncol=5,
        loc="lower center",
        bbox_to_anchor=(1.7, -0.4),
    )
    fig.savefig(out, dpi=300)


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
    return parser.parse_args()


def main():
    args = _parse_args()
    if args.gc:
        generate_clade_progression(args.c_file, args.c_out)
    if args.gg:
        generate_geoclade_progression(args.g_file, args.g_out)


if __name__ == "__main__":
    main()
# r"./data/input_sample/clade_progression.tsv"
