#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd
from matplotlib.ticker import FuncFormatter
from collections import Counter
from matplotlib.ticker import FuncFormatter, MaxNLocator 

x = np.arange(5)
total_width, n = 0.5, 5
width = total_width / n
x = x - (total_width - width) * 5

x_label=[8,32,64,128,256]

df = pd.read_excel(r"D:\HuaweiMoveData\Users\huawei\Desktop\better.xlsx")



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
def draw_step_fig(row_idx, y_Me, y_Metis, y_mtmetis, y_ParMetis, y_chaco, y_kahip, mat_name, fig):
    width = 0.3
    ax = fig.add_subplot(3, 5, row_idx + 1)  # 调整子图布局为4行5列

    ax.plot(x, y_Metis, marker='o', markersize=2, linewidth=2, label='Metis', color='#2800d7')
    ax.plot(x, y_mtmetis, marker='o', markersize=2, linewidth=2, label='mt-metis', color='#1982c4')
    ax.plot(x, y_ParMetis, marker='o', markersize=2, linewidth=2, label='ParMetis', color='#52a675')
    ax.plot(x, y_chaco, marker='o', markersize=2, linewidth=2, label='Chaco', color='#8ac926')
    ax.plot(x, y_kahip, marker='o', markersize=2, linewidth=2, label='KaHIP', color='#ffca3a')
    ax.plot(x, y_Me, marker='o', markersize=2, linewidth=2, label='Hunyuan', color='#f05006')

    plt.yticks(fontsize=12)
    plt.xticks(x + 3 * width-0.9 , x_label, fontsize=16)

    abs_data = [abs(value) for value in y_Me + y_Metis + y_mtmetis + y_ParMetis + y_chaco + y_kahip]
    median_value = np.median(abs_data)
    print(median_value)  # 打印 median_value
    if median_value > 1e7:
        data = 1e7
        exponent = '7'
    elif median_value > 1e6:
        data = 1e6
        exponent = '6'
    elif median_value > 1e5:
        data = 1e5
        exponent = '5'
    else:
        data = 1e4
        exponent = '4'

    def custom_formatter(x, pos):
        return '{:.2f}'.format(x / data)

    ax.yaxis.set_major_formatter(FuncFormatter(custom_formatter))
    ax.text(0, 1.02, '$\mathdefault{×1e' + exponent + '}$', transform=ax.transAxes,fontsize=11)

    plt.title(mat_name, fontsize=20, pad=2)
    legend = ax.legend(handletextpad=0.2, handleheight=1.5, columnspacing=1.5, labelspacing=0.2, handlelength=0.8,
                       loc='upper left', fontsize=11, frameon=True)
    legend.get_frame().set_edgecolor('grey')
    legend.get_frame().set_linewidth(1.5)


    if index == 2:
        plt.ylabel('Edge-cut', fontsize=25)
        ax.yaxis.set_label_coords(-2.45, -0.7)  # 将y轴标签调整到指定位置，紧贴边界
    bwith=1.5
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)


    ax.yaxis.set_major_locator(MaxNLocator(nbins=6))  # 设置纵坐标最大刻度为5个








fig = plt.figure(figsize=(28, 12))


titles = [
    "vas_stokes_4M", "rgg_n_2_24_s0", "Queen_4147", "stokes", "HV15R", 
     "Cube_Coup_dt6", "Flan_1565", "ML_Geer", "dielFilterV3real",
    "hugebubbles-00020", "Geo_1438", "nv2", "ldoor", "ss", "Emilia_923",
    
]

mat_num = 15

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
        d.append(float(ParMetis[index * 5 + j]))
        e.append(float(chaco[index * 5 + j]))
        f.append(float(kahip[index * 5 + j]))

    draw_step_fig(index, a, b, c, d, e, f, titles[index], fig)
    
fig.text(0.55, 0.14, 'Number of partitions', ha='center', fontsize=25)  # 在整个图的底部中央添加横坐标解释

plt.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.9, wspace=0.15, hspace=0.3)
plt.gcf().subplots_adjust(left=0.15, top=1, bottom=0.2)


plt.show()
fig.savefig('cut2.eps', dpi=300, format='eps', bbox_inches='tight', pad_inches=0.1)

