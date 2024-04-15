import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# 从Excel文件读取数据
df1 = pd.read_excel('metis_32.xlsx', skiprows=1)
df2 = pd.read_excel('mt-metis_32.xlsx', skiprows=1)  # 假设数据从第二行开始
df3 = pd.read_excel('parmetis_32.xlsx', skiprows=1)
df4 = pd.read_excel('chaco_32.xlsx', skiprows=1)
df5 = pd.read_excel('kahip_32.xlsx', skiprows=1)

# 创建一个2行5列的布局
fig, axes = plt.subplots(nrows=2, ncols=5, figsize=(12, 4), sharey=True)

# 散点图标题和标签
labels = ['Metis vs Hunyuan', 'mt-metis vs Hunyuan', 'ParMetis vs Hunyuan', 'Chaco vs Hunyuan', 'KaHIP vs Hunyuan']
legends = ['Metis', 'mt-metis', 'ParMetis', 'Chaco', 'KaHIP']

color1 = '#2800d7'  # Metis
color2 = '#1982c4'  # mt-metis
color3 = '#52a675'  # ParMetis
color4 = '#8ac926'  # Chaco
color5 = '#ffca3a'  # KaHIP
color6 = '#f05006'  # Hunyuan

# 绘制每个散点图
for i, ax in enumerate(axes.flatten()):

    j = i % 5
    if j == 0:
        df = df1
        col = color1
    elif j == 1:
        df = df2
        col = color2
    elif j == 2:
        df = df3
        col = color3
    elif j == 3:
        df = df4
        col = color4
    elif j == 4:
        df = df5
        col = color5
    # 提取所需列的数据
    nvtxs = df.iloc[:, 1].tolist()  # 第2列数据为顶点数
    nedges = df.iloc[:, 2].tolist()  # 第3列数据为边数
    mt_edgecut = df.iloc[:, 3].tolist()  # 第4列数据为mt-metis边切值
    mt_alltimes = df.iloc[:, 4].tolist()  # 第5列数据为mt-metis总时间
    # me_edgecut = df.iloc[:, 9].tolist()  # 第10列数据为hunyuan边切值
    me_edgecut = df.iloc[:, 15].tolist()  # 第10列数据为hunyuan边切值
    me_alltimes = df.iloc[:, 10].tolist()  # 第11列数据为hunyuan总时间

    sort_idx = np.argsort(nvtxs)
    nvtxs = np.array(nvtxs)[sort_idx].tolist()
    nedges = np.array(nedges)[sort_idx].tolist()
    mt_edgecut = np.array(mt_edgecut)[sort_idx].tolist()
    mt_alltimes = np.array(mt_alltimes)[sort_idx].tolist()
    me_edgecut = np.array(me_edgecut)[sort_idx].tolist()
    me_alltimes = np.array(me_alltimes)[sort_idx].tolist()

    if i < 5:
        x = np.log10(np.array(nvtxs).astype(float))
        y1 = np.log10(mt_edgecut)
        y2 = np.log10(me_edgecut)
        ax.scatter(x, y2, color=color6, label='Hunyuan', s=3, zorder=2)
        ax.scatter(x, y1, color=col, label=legends[i], s=3, zorder=2)
        # 设置图例位置和对齐方式
        ax.legend(loc='upper center', bbox_to_anchor=(0.82, 0.35), ncol=1, prop={'size': 8}, handlelength=0.5, handletextpad=0.5)
    else:
        x = np.log10(nvtxs)
        y1 = np.log10(mt_alltimes)
        y2 = np.log10(me_alltimes)
        ax.scatter(x, y2, color=color6, label='Hunyuan', s=3, zorder=2)
        ax.scatter(x, y1, color=col, label=legends[i-5], s=3, zorder=2)
        # 设置图例位置和对齐方式
        ax.legend(loc='upper center', bbox_to_anchor=(0.82, 0.35), ncol=1, prop={'size': 8}, handlelength=0.5, handletextpad=0.5)
    
    if i == 0:
        ax.set_ylabel("Log10 Edgecut", labelpad=5)  # 根据需要调整labelpad的值
    elif i == 5:
        ax.set_ylabel("Log10 Total time", labelpad=5)  # 根据需要调整labelpad的值
    
    if i < 5:
        # 调整标题位置
        ax.title.set_position([0.5, 2.5])
        # 设置标题居中
        ax.title.set_ha("center")
        ax.set_title(labels[i-5])

# 设置整体布局
plt.tight_layout(rect=[0, 0.05, 1, 0.95])

# 添加整体标题
# plt.suptitle('2行5列散点图布局', fontsize=16, fontweight='bold')

# 添加共享的x轴标签
fig.text(0.5, 0.02, 'Log10 Number of Vertices', ha='center', fontsize=12)

# 显示图形
plt.savefig("scatter.pdf", dpi=400, format="pdf", bbox_inches='tight', pad_inches=0.0)