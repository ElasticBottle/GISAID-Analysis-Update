#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import seaborn as sns

# Your x and y axis
x = range(1, 6)
y = [[10, 4, 6, 5, 3], [12, 2, 7, 10, 1], [8, 18, 5, 7, 6]]

# use a known color palette (see..)
pal = sns.color_palette("Set1")
plt.stackplot(x, y, labels=["A", "B", "C"], colors=pal, alpha=0.4)
plt.legend(loc="upper right")
plt.show()

# create your palette
pal = ["#9b59b6", "#e74c3c", "#34495e", "#2ecc71"]
plt.stackplot(x, y, labels=["A", "B", "C"], colors=pal, alpha=0.4)
plt.legend(loc="upper center")
# plt.show()

# Make data
data = pd.DataFrame(
    {
        "group_A": [1, 4, 6, 8, 9],
        "group_B": [2, 24, 7, 10, 12],
        "group_C": [2, 8, 5, 10, 6],
    },
    index=range(1, 6),
)

# We need to transform the data from raw data to percentage (fraction)
data_perc = data.divide(data.sum(axis=1), axis=0)

# Make the plot
plt.stackplot(
    range(1, 6),
    data_perc["group_A"],
    data_perc["group_B"],
    data_perc["group_C"],
    labels=["A", "B", "C"],
)
plt.legend(loc="upper left")
plt.margins(0, 0)
plt.title("100 % stacked area chart")
plt.show()

#%%
df = pd.read_csv(
    r"./data/input_sample/clade_progression.tsv", sep="\t", skiprows=1, index_col=0,
)
indiv_df = df.loc[:"Cumulative:", :].iloc[:-1]
indiv_df = indiv_df.astype(float)
hor_sum = indiv_df.sum(axis=1)
df_percentage = indiv_df.div(hor_sum, axis=0) * 100

# cum_df = df.loc["Month":, :].iloc[1:]
# cum_df = cum_df.astype(float)

values = [None] * 7
insert_order = {
    "O": 0,
    "S": 1,
    "L": 2,
    "V": 3,
    "G": 4,
    "GH": 5,
    "GR": 6,
}
for label, value in df_percentage.items():
    values[insert_order.get(label, 0)] = value
pal = ["#f291ee"]
plt.stackplot(
    df_percentage.index,
    *values,
    labels=sorted(df.columns, key=lambda x: insert_order.get(x, 0)),
)
y_ticks = list(map(lambda x: f"{x:.0f}%", np.arange(110, step=10)))
plt.yticks(
    ticks=np.arange(110, step=10), labels=y_ticks,
)
plt.margins(0.0)
plt.legend(ncol=7, loc="lower center", bbox_to_anchor=(0.5, -0.3))
