import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter

df = pd.read_excel('time.xlsx')

# Getting the values for each category
Me = []
Metis = []
mt_metis = []
ParMetis = []
chaco = []
kahip = []

for i in range(5):
    index = i * 6
    data = df.iloc[index, 2:17].values.tolist()
    Me.extend([float(val) for val in data])
    data = df.iloc[index + 1, 2:17].values.tolist()
    Metis.extend([float(val) for val in data])
    data = df.iloc[index + 2, 2:17].values.tolist()
    mt_metis.extend([float(val) for val in data])
    data = df.iloc[index + 3, 2:17].values.tolist()
    ParMetis.extend([float(val) for val in data])
    data = df.iloc[index + 4, 2:17].values.tolist()
    chaco.extend([float(val) for val in data])
    data = df.iloc[index + 5, 2:17].values.tolist()
    kahip.extend([float(val) for val in data])

def format_func(value, tick_number):
    if value >= 1e4:
        return f"${{{'{:.0e}'.format(value).replace('e+0', 'x10^')}}}$"
    else:
        return f"${{{int(value):,}}}$"

# Creating a 3x5 grid layout
fig, axes = plt.subplots(nrows=3, ncols=5, figsize=(17, 8.5), sharex=True)

# Bar chart titles and labels
titles = [
    "vas_stokes_4M", "rgg_n_2_24_s0", "Queen_4147", "stokes", "HV15R", 
    "Cube_Coup_dt6", "Flan_1565", "ML_Geer", "dielFilterV3real", "hugebubbles",
    "Geo_1438", "nv2", "ldoor", "ss", "Emilia_923"]
legends = ['Metis', 'mt-metis', 'ParMetis', 'Chaco', 'KaHIP', 'Hunyuan']
partitions = ['8', '32', '64', '128', '256']

# Setting the bar width
bar_width = 0.14

colors = ['#2800d7', '#1982c4', '#52a675', '#8ac926', '#ffca3a', '#f05006']

handles = []

# Plotting each vertical bar chart
for i, ax in enumerate(axes.flatten()):
    a = []
    b = []
    c = []
    d = []
    e = []
    f = []

    for j in range(5):
        a.append(float(Me[i * 5 + j]))
        b.append(float(Metis[i * 5 + j]))
        c.append(float(mt_metis[i * 5 + j]))
        d.append(float(ParMetis[i * 5 + j]))
        e.append(float(chaco[i * 5 + j]))
        f.append(float(kahip[i * 5 + j]))

    index = np.arange(len(partitions))

    for j in range(6):
        handles.append(ax.bar(0, 0, color=colors[j], label=legends[j], zorder=0))

    ax.bar(index, b, bar_width, color=colors[0], label='Metis', edgecolor='black', zorder=2)
    ax.bar(index + bar_width, c, bar_width, color=colors[1], label='mt-metis', edgecolor='black', zorder=2)
    ax.bar(index + bar_width * 2, d, bar_width, color=colors[2], label='ParMetis', edgecolor='black', zorder=2)
    ax.bar(index + bar_width * 3, e, bar_width, color=colors[3], label='Chaco', edgecolor='black', zorder=2)
    ax.bar(index + bar_width * 4, f, bar_width, color=colors[4], label='KaHIP', edgecolor='black', zorder=2)
    ax.bar(index + bar_width * 5, a, bar_width, color=colors[5], label='Hunyuan', edgecolor='black', zorder=2)
    
    ax.grid(axis='y', linestyle='--', alpha=0.9, zorder=1)

    ax.set_title(titles[i],fontsize=16)
    formatter = ScalarFormatter(useMathText=True)
    ax.yaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(format_func))
    ax.yaxis.get_offset_text().set_position((0.2, 1.0))
    ax.set_xticks(index + bar_width * 2.5)  # Adding ticks
    ax.set_xticklabels(partitions, fontsize=12)  # Setting tick labels
    ax.tick_params(axis='y', labelsize=12)

# Adjusting the layout of the subplots
fig.subplots_adjust(hspace=0.4, wspace=0.1)

# Adding a shared x-axis label
fig.text(0.5, 0.03, 'Number of Partitions', ha='center', fontsize=18)

# 添加共享的y轴标签
fig.text(0.0, 0.5, 'Speedup', ha='center', va='center', fontsize=18, rotation=90)

fig.legend(handles=handles, labels=legends, loc='upper center', bbox_to_anchor=(0.5, 1.0), fontsize=14, handlelength=1, handletextpad=1, ncol=6)

# Setting the overall layout
plt.tight_layout(rect=[0, 0.05, 1, 0.95])

# Displaying the chart
plt.savefig("speedup.pdf", dpi=400, format="pdf", bbox_inches='tight', pad_inches=0.0)