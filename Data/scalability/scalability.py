import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df0 = pd.read_excel('20240321metis_32.xlsx',sheet_name='Ave')
df1 = pd.read_excel('20240306_32.xlsx',sheet_name='Ave')
df2 = pd.read_excel('20240306_32.xlsx',sheet_name='Min')
df3 = pd.read_excel('20240401_32.xlsx',sheet_name='Ave')
df4 = pd.read_excel('20240401_32.xlsx',sheet_name='Min')
df0 = df1

all_avemetis = df0.iloc[0:15, 4].values.tolist()  # 注意索引从 0 开始，第 5 列的索引为 4
all_ave4090 = df1.iloc[0:15, 4].values.tolist()  # 注意索引从 0 开始，第 5 列的索引为 4
all_min4090 = df2.iloc[0:15, 4].values.tolist()  # 注意索引从 0 开始，第 5 列的索引为 4
all_ave4070 = df3.iloc[0:15, 4].values.tolist()  # 注意索引从 0 开始，第 5 列的索引为 4
all_min4070 = df4.iloc[0:15, 4].values.tolist()  # 注意索引从 0 开始，第 5 列的索引为 4

coarsen_avemetis = df0.iloc[0:15, 5].values.tolist()
coarsen_ave4090 = df1.iloc[0:15, 5].values.tolist()
coarsen_min4090 = df2.iloc[0:15, 5].values.tolist()
coarsen_ave4070 = df3.iloc[0:15, 5].values.tolist()
coarsen_min4070 = df4.iloc[0:15, 5].values.tolist()

uncoarsen_avemetis = df0.iloc[0:15, 7].values.tolist()
uncoarsen_ave4090 = df1.iloc[0:15, 7].values.tolist()
uncoarsen_min4090 = df2.iloc[0:15, 7].values.tolist()
uncoarsen_ave4070 = df3.iloc[0:15, 7].values.tolist()
uncoarsen_min4070 = df4.iloc[0:15, 7].values.tolist()

for i in range(len(all_avemetis)):
    coarsen_ave4070[i] += uncoarsen_ave4070[i]
    coarsen_ave4090[i] += uncoarsen_ave4090[i]
    coarsen_min4070[i] += uncoarsen_min4070[i]
    coarsen_min4090[i] += uncoarsen_min4090[i]

    coarsen_ave4070[i] = coarsen_ave4070[i] / coarsen_ave4090[i]
    coarsen_ave4090[i] = coarsen_ave4090[i] / coarsen_ave4090[i]
    coarsen_min4070[i] = coarsen_min4070[i] / coarsen_min4090[i]
    coarsen_min4090[i] = coarsen_min4090[i] / coarsen_min4090[i]

    uncoarsen_ave4070[i] = uncoarsen_ave4070[i] / uncoarsen_ave4090[i]
    uncoarsen_ave4090[i] = uncoarsen_ave4090[i] / uncoarsen_ave4090[i]
    uncoarsen_min4070[i] = uncoarsen_min4070[i] / uncoarsen_min4090[i]
    uncoarsen_min4090[i] = uncoarsen_min4090[i] / uncoarsen_min4090[i]

labels = ['vas_stokes_4M', 'rgg_n_2_24_s0', 'Queen_4147', 'stokes', 'HV15R', 
          'Cube_Coup_dt6', 'Flan_1565', 'ML_Geer', 'dielFilterV3real', 'hugebubbles', 
          'Geo_1438', 'nv2', 'ldoor','ss', 'Emilia_923']

# 设置图形大小
fig, ax = plt.subplots(figsize=(7, 4.5))

# 设置柱子宽度
bar_width = 0.35

# 设置柱子位置
index = np.arange(len(labels))

color1 = '#eaa56e'  # GTX-4090
color2 = '#a1c4e0'  # GTX-4070Ti
color3 = '#a1c4e0'

alpha1 = 1

ax.bar(index, coarsen_ave4090, bar_width, color=color1, alpha=alpha1, label='GTX 4090', edgecolor='black', zorder=2)
ax.bar(index + bar_width, coarsen_ave4070, bar_width, color=color2, alpha=alpha1, label='GTX 4070Ti', edgecolor='black', zorder=2)

# 在“stokes”的“GTX-4070Ti”柱子上添加斜线阴影
stokes_index = labels.index('stokes')
ax.bar(stokes_index + bar_width, 2, bar_width, color=color3, hatch='////', alpha=alpha1, label='GTX 4070Ti  (Out of Memory)', edgecolor='black', zorder=2)

# 添加横坐标网格线 zorder=1 代表图层
plt.grid(axis='y', linestyle='--', alpha=0.9, zorder=1)

# 设置标签和标题

ax.set_ylabel('Execution Time Relative to 4090', fontsize=12)

# 设置x轴刻度标签
ax.set_xticks(index + bar_width * 0.5)
ax.set_xticklabels(labels, rotation=45, ha='right')

# 添加图例
legend_labels = ['GTX 4090', 'GTX 4070Ti', 'GTX 4070Ti (Out of Memory)']
legend = ax.legend(labels=legend_labels, ncol=3, bbox_to_anchor=(0.5, 1.14), prop={'size': 10}, handlelength=1, handletextpad=1, loc='upper center')

# 设置整体布局
plt.tight_layout(rect=[0, 0.05, 1, 0.95])

# 显示图形
plt.savefig("scalability.pdf", dpi=400, format="pdf", bbox_inches='tight', pad_inches=0.0)