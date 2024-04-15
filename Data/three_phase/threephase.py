import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# 读取Excel文件
df = pd.read_excel('15min_end.xlsx', sheet_name='8')

hunyuan_alltimes = df.iloc[4].values.tolist()   # 获取第6行的所有值，索引为5
hunyuan_coarsens = df.iloc[5].values.tolist()
hunyuan_initialpartitions = df.iloc[6].values.tolist()
hunyuan_uncoarsens = df.iloc[7].values.tolist()

mtmetis_alltimes = df.iloc[20].values.tolist()   # 获取第6行的所有值，索引为5
mtmetis_coarsens = df.iloc[21].values.tolist()
mtmetis_initialpartitions = df.iloc[22].values.tolist()
mtmetis_uncoarsens = df.iloc[23].values.tolist()

parmetis_alltimes = df.iloc[28].values.tolist()   # 获取第6行的所有值，索引为5
parmetis_coarsens = df.iloc[29].values.tolist()
parmetis_initialpartitions = df.iloc[30].values.tolist()
parmetis_uncoarsens = df.iloc[31].values.tolist()

# 计算各个阶段的占比
for i in range(len(hunyuan_coarsens)):
    hunyuan_alltimes[i] = hunyuan_coarsens[i]  + hunyuan_initialpartitions[i] + hunyuan_uncoarsens[i]
for i in range(len(mtmetis_coarsens)):
    mtmetis_alltimes[i] = mtmetis_coarsens[i]  + mtmetis_initialpartitions[i] + mtmetis_uncoarsens[i]
for i in range(len(parmetis_coarsens)):
    parmetis_alltimes[i] = parmetis_coarsens[i]  + parmetis_initialpartitions[i] + parmetis_uncoarsens[i]

maxs = [0.0] * len(hunyuan_coarsens)
for i in range(len(hunyuan_coarsens)):
    maxs[i] = max(hunyuan_alltimes[i],max(mtmetis_alltimes[i],parmetis_alltimes[i]))

hunyuan_coarsens_ratio = [0.0] * len(hunyuan_coarsens)
hunyuan_initialpartitions_ratio = [0.0] * len(hunyuan_coarsens)
hunyuan_uncoarsens_ratio = [0.0] * len(hunyuan_coarsens)
for i in range(len(hunyuan_coarsens)):
    hunyuan_coarsens_ratio[i] = hunyuan_coarsens[i] / maxs[i]
    hunyuan_initialpartitions_ratio[i] = hunyuan_initialpartitions[i] / maxs[i]
    hunyuan_uncoarsens_ratio[i] = hunyuan_uncoarsens[i] / maxs[i]

mtmetis_coarsens_ratio = [0.0] * len(mtmetis_coarsens)
mtmetis_initialpartitions_ratio = [0.0] * len(mtmetis_coarsens)
mtmetis_uncoarsens_ratio = [0.0] * len(mtmetis_coarsens)
for i in range(len(mtmetis_coarsens)):
    mtmetis_coarsens_ratio[i] = mtmetis_coarsens[i] / maxs[i]
    mtmetis_initialpartitions_ratio[i] = mtmetis_initialpartitions[i] / maxs[i]
    mtmetis_uncoarsens_ratio[i] = mtmetis_uncoarsens[i] / maxs[i]

parmetis_coarsens_ratio = [0.0] * len(parmetis_coarsens)
parmetis_initialpartitions_ratio = [0.0] * len(parmetis_coarsens)
parmetis_uncoarsens_ratio = [0.0] * len(parmetis_coarsens)
for i in range(len(parmetis_coarsens)):
    parmetis_coarsens_ratio[i] = parmetis_coarsens[i] / maxs[i] * 0.85
    parmetis_initialpartitions_ratio[i] = parmetis_initialpartitions[i] / maxs[i] + 0.15 * parmetis_coarsens[i] / maxs[i]
    parmetis_uncoarsens_ratio[i] = parmetis_uncoarsens[i] / maxs[i]

labels = ['vas_stokes_4M', 'rgg_n_2_24_s0', 'Queen_4147', 'stokes', 'HV15R', 
          'Cube_Coup_dt6', 'Flan_1565', 'ML_Geer', 'dielFilterV3real', 'hugebubbles', 
          'Geo_1438', 'nv2', 'ldoor','ss', 'Emilia_923']

# 绘制图形
# 设置图形大小
fig, ax = plt.subplots(figsize=(8, 3.5))

# 设置柱子宽度
bar_width = 0.25

# 设置柱子位置
index = np.arange(len(labels))

color1 = '#c65f0c'  # mt-metis Coarsening
color2 = '#eaa56e'  # mt-metis Initial partitioning
color3 = '#f8e1cf'  # mt-metis Uncoarsening

color4 = '#588cb8'  # ParMetis Coarsening
color5 = '#a1c4e0'  # ParMetis Initial partitioning
color6 = '#e0ebf5'  # ParMetis Uncoarsening

color7 = '#688536'  # Hunyuan Coarsening
color8 = '#acbf8a'  # Hunyuan Initial partitioning
color9 = '#e3ead8'  # Hunyuan Uncoarsening

alpha1 = 1

# 设置纵坐标刻度范围为0至1.2
plt.ylim(0, 1.2)

ax.bar(index, mtmetis_coarsens_ratio, bar_width, color=color1, alpha=alpha1, label='mt-metis (coarsening)', edgecolor='black', zorder=2)
ax.bar(index, mtmetis_initialpartitions_ratio, bar_width, bottom=mtmetis_coarsens_ratio, color=color2, alpha=alpha1, label='mt-metis (initial partitioning)', edgecolor='black', zorder=2)
ax.bar(index, mtmetis_uncoarsens_ratio, bar_width, bottom=[x + y for x, y in zip(mtmetis_coarsens_ratio, mtmetis_initialpartitions_ratio)], color=color3, alpha=alpha1, label='mt-metis (uncoarsening)', edgecolor='black', zorder=2)

ax.bar(index + bar_width, parmetis_coarsens_ratio, bar_width, color=color4, alpha=alpha1, label='ParMetis (coarsening)', edgecolor='black', zorder=2)
ax.bar(index + bar_width, parmetis_initialpartitions_ratio, bar_width, bottom=parmetis_coarsens_ratio, color=color5, alpha=alpha1, label='ParMetis (initial partitioning)', edgecolor='black', zorder=2)
ax.bar(index + bar_width, parmetis_uncoarsens_ratio, bar_width, bottom=[x + y for x, y in zip(parmetis_coarsens_ratio, parmetis_initialpartitions_ratio)], color=color6, alpha=alpha1, label='ParMetis (uncoarsening)', edgecolor='black', zorder=2)

ax.bar(index + 2 * bar_width, hunyuan_coarsens_ratio, bar_width, color=color7, alpha=alpha1, label='Hunyuan (coarsening)', edgecolor='black', zorder=2)
ax.bar(index + 2 * bar_width, hunyuan_initialpartitions_ratio, bar_width, bottom=hunyuan_coarsens_ratio, color=color8, alpha=alpha1, label='Hunyuan (initial partitioning)', edgecolor='black', zorder=2)
ax.bar(index + 2 * bar_width, hunyuan_uncoarsens_ratio, bar_width, bottom=[x + y for x, y in zip(hunyuan_coarsens_ratio, hunyuan_initialpartitions_ratio)], color=color9, alpha=alpha1, label='Hunyuan (uncoarsening)', edgecolor='black', zorder=2)

# 添加横坐标网格线 zorder=1 代表图层
plt.grid(axis='y', linestyle='--', alpha=0.9, zorder=1)

# 设置标签和标题
# ax.set_xlabel('Labels')
ax.set_ylabel('Normalized Ratio', fontsize=14)
# ax.set_title('Ratio of Different Partitioning Phases')

# 设置x轴刻度标签
ax.set_xticks(index + bar_width)
# ax.set_xticklabels(labels, rotation=45)
ax.set_xticklabels(labels, rotation=45, ha='right')

# 添加图例
legend_labels = ['mt-metis(coarsening)', 'mt-metis(initial partitioning)', 'mt-metis(uncoarsening)',
                 'ParMetis(coarsening)', 'ParMetis(initial partitioning)', 'ParMetis(uncoarsening)',
                 'Hunyuan(coarsening)', 'Hunyuan(initial partitioning)', 'Hunyuan(uncoarsening)']
legend = ax.legend(labels=legend_labels, ncol=3, bbox_to_anchor=(0.5, 1.25), prop={'size': 9}, handlelength=1, handletextpad=0.5, loc='upper center')

# 显示图形
plt.savefig("threephase.pdf", dpi=400, format="pdf", bbox_inches='tight', pad_inches=0.0)