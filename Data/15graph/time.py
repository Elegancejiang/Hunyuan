import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import MaxNLocator

x = np.arange(5)
total_width, n = 0.5, 5
width = total_width / n
x = x - (total_width - width) * 5

x_label=[8,32,64,128,256]

df = pd.read_excel('time.xlsx')

#读取数据
Me=[]
Metis=[]
mt_metis=[]
ParMetis=[]
chaco=[]
kahip=[]

for i in range(3):
    for j in range(5):
        index = i * 5 + j
        Me.append(float(df.iloc[0, index + 2]))
        Me.append(float(df.iloc[6, index + 2]))
        Me.append(float(df.iloc[12, index + 2]))
        Me.append(float(df.iloc[18, index + 2]))
        Me.append(float(df.iloc[24, index + 2]))

        Metis.append(float(df.iloc[1, index + 2]))
        Metis.append(float(df.iloc[7, index + 2]))
        Metis.append(float(df.iloc[13, index + 2]))
        Metis.append(float(df.iloc[19, index + 2]))
        Metis.append(float(df.iloc[25, index + 2]))
   
        mt_metis.append(float(df.iloc[2, index + 2]))
        mt_metis.append(float(df.iloc[8, index + 2]))
        mt_metis.append(float(df.iloc[14, index + 2]))
        mt_metis.append(float(df.iloc[20, index + 2]))
        mt_metis.append(float(df.iloc[26, index + 2]))
       
        ParMetis.append(float(df.iloc[3, index + 2]))
        ParMetis.append(float(df.iloc[9, index + 2]))
        ParMetis.append(float(df.iloc[15, index + 2]))
        ParMetis.append(float(df.iloc[21, index + 2]))
        ParMetis.append(float(df.iloc[27, index + 2]))

        chaco.append(float(df.iloc[4, index + 2]))
        chaco.append(float(df.iloc[10, index + 2]))
        chaco.append(float(df.iloc[16, index + 2]))
        chaco.append(float(df.iloc[22, index + 2]))
        chaco.append(float(df.iloc[28, index + 2]))

        kahip.append(float(df.iloc[5, index + 2]))
        kahip.append(float(df.iloc[11, index + 2]))
        kahip.append(float(df.iloc[17, index + 2]))
        kahip.append(float(df.iloc[23, index + 2]))
        kahip.append(float(df.iloc[29, index + 2]))

def draw_step_fig(row_idx,y_Me, y_Metis, y_mtmetis, y_ParMetis, y_chaco,  y_kahip , mat_name, fig):
    width = 0.1
    ax = fig.add_subplot(3, 5, row_idx + 1)  

    ax.bar(x, y_Metis, width=width, label='Metis', color='#2800d7')
    ax.bar(x + width, y_mtmetis, width=width, label='mt-metis', color='#1982c4')
    ax.bar(x + 2*width, y_ParMetis, width=width, label='ParMetis', color='#52a675')
    ax.bar(x + 3*width, y_chaco, width=width, label='Chaco', color='#8ac926')
    ax.bar(x + 4*width, y_kahip, width=width, label='KaHIP', color='#ffca3a')
    ax.bar(x + 5*width, y_Me, width=width, label='Hunyuan', color='#f05006')

    plt.yticks(fontsize = 13)
    plt.xticks(x + 3 * width, x_label, fontsize=15)  # 设置横坐标刻度和标签

    plt.title(mat_name, fontsize=22, pad=2)
    if row_idx == 0:
        ax.legend(handletextpad=0.2, handleheight=0.5, columnspacing=1.5, labelspacing=0.2, handlelength=0.8,
                  loc='upper center', bbox_to_anchor=(0.55, 1.1), bbox_transform=plt.gcf().transFigure,
                  ncol=9, fancybox=False, shadow=False, fontsize=20, frameon=False)

    if index == 2:
        plt.ylabel('Speedup', fontsize=25)
        ax.yaxis.set_label_coords(-2.35, -0.7)  

    bwith=1.5
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=6))  

fig = plt.figure(figsize=(28, 12))  # 调整整个图形的大小

mat_num = 15

titles = [
    "vas_stokes_4M", "rgg_n_2_24_s0", "Queen_4147", "stokes", "HV15R", 
    "Cube_Coup_dt6", "Flan_1565", "ML_Geer", "dielFilterV3real", "hugebubbles-00020",
    "Geo_1438", "nv2", "ldoor", "ss", "Emilia_923"]

for index in range(mat_num):
    a = []
    b = []
    c = []
    d = []
    e = []
    f = []

    # 将原来的索引方式改为按列读取数据
    for j in range(5):
        a.append(float(Me[index * 5 + j]))
        b.append(float(Metis[index * 5 + j]))
        c.append(float(mt_metis[index * 5 + j]))
        d.append(float(ParMetis[index* 5 + j]))
        e.append(float(chaco[index * 5 + j]))
        f.append(float(kahip[index * 5 + j]))

    draw_step_fig(index, a, b, c, d, e, f, titles[index], fig)
fig.text(0.55, 0.14, 'Number of partitions', ha='center', fontsize=25) 


plt.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.9, wspace=0.1, hspace=0.3)
plt.gcf().subplots_adjust(left=0.15, top=1, bottom=0.2)

plt.show()
# fig.savefig('cut2.eps', dpi=300, format='eps', bbox_inches='tight', pad_inches=0.1)

# 显示图形
plt.savefig("time.pdf", dpi=400, format="pdf", bbox_inches='tight', pad_inches=0.0)


