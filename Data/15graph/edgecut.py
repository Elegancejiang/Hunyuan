import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter

df = pd.read_excel('edge.xlsx')

# 获取指定单元格范围的数值
Me=[]
Metis=[]
mt_metis=[]
ParMetis=[]
chaco=[]
kahip=[]

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

# 创建一个3行5列的布局
fig, axes = plt.subplots(nrows=3, ncols=5, figsize=(17, 8.5), sharey=True)

# 柱状图标题和标签
titles = [
    "vas_stokes_4M", "rgg_n_2_24_s0", "Queen_4147", "stokes", "HV15R", 
    "Cube_Coup_dt6", "Flan_1565", "ML_Geer", "dielFilterV3real", "hugebubbles",
    "Geo_1438", "nv2", "ldoor", "ss", "Emilia_923"]
legends = ['Metis', 'mt-metis', 'ParMetis', 'Chaco', 'KaHIP', 'Hunyuan']
partitions = ['8', '32', '64', '128', '256']

# 设置柱子宽度
bar_width = 0.14

color1 = '#2800d7'  # Metis
color2 = '#1982c4'  # mt-metis
color3 = '#52a675'  # ParMetis
color4 = '#8ac926'  # Chaco
color5 = '#ffca3a'  # KaHIP
color6 = '#f05006'  # Hunyuan

# 绘制每个横向柱状图
for i, ax in enumerate(axes.flatten()):
    a = []
    b = []
    c = []
    d = []
    e = []
    f = []

    for j in range(5):
        a.append(float(Me[j * 15 + i]))
        b.append(float(Metis[j * 15 + i]))
        c.append(float(mt_metis[j * 15 + i]))
        d.append(float(ParMetis[j * 15 + i]))
        e.append(float(chaco[j * 15 + i]))
        f.append(float(kahip[j * 15 + i]))

    index = np.arange(len(a))

    ax.barh(index, b, bar_width, color=color1, label='Metis', edgecolor='black', zorder=2)
    ax.barh(index + bar_width, c, bar_width, color=color2, label='mt-metis', edgecolor='black', zorder=2)
    ax.barh(index + bar_width * 2, d, bar_width, color=color3, label='ParMetis', edgecolor='black', zorder=2)
    ax.barh(index + bar_width * 3, e, bar_width, color=color4, label='Chaco', edgecolor='black', zorder=2)
    ax.barh(index + bar_width * 4, f, bar_width, color=color5, label='KaHIP', edgecolor='black', zorder=2)
    ax.barh(index + bar_width * 5, a, bar_width, color=color6, label='Hunyuan', edgecolor='black', zorder=2)
    
    ax.legend(loc='lower right', bbox_to_anchor=(1.02, -0.03), ncol=1, fontsize=9.5, handlelength=0.5, handletextpad=0.5)

    ax.grid(axis='x', linestyle='--', alpha=0.9, zorder=1)

    ax.set_title(titles[i], fontsize=16)
    formatter = ScalarFormatter(useMathText=True)
    ax.xaxis.set_major_formatter(formatter)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
    ax.xaxis.get_offset_text().set_position((1.0, 0.2))
    ax.set_yticks(index + bar_width * 2.5)  # 添加刻度
    ax.set_yticklabels(partitions, fontsize=12)  # 设置刻度标签
    ax.tick_params(axis='x', labelsize=12)

# 调整子图的布局
fig.subplots_adjust(hspace=0.4, wspace=0.1)

# 添加共享的x轴标签
fig.text(0.5, 0.03, 'Edge-cut Values', ha='center', fontsize=18)

# 添加共享的y轴标签
fig.text(0.0, 0.5, 'Number of Partitions', ha='center', va='center', fontsize=18, rotation=90)

# 设置整体布局
plt.tight_layout(rect=[0, 0.05, 1, 0.95])

# 显示图形
plt.savefig("edgecut.pdf", dpi=400, format="pdf", bbox_inches='tight', pad_inches=0.0)