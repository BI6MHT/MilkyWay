# -*- coding: utf-8 -*-
"""
Created on Fri Jan 14 23:32:00 2022

@author: BI6MHT
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 

# 太阳距离银心的距离和旋转速度
R0 = 8.5 #kpc
V0 = 220 #km/s

# 读取文件
table = pd.read_csv('MilkyWay.csv')
#返回行数
rows_num = table.shape[0]

# Outfile

R_list = []
V_list = []

for i in range(0,rows_num):
    # 读取第一行的内容
    myRow = table.iloc[i,:]
    # 去掉赤经赤纬
    myRow = myRow[2:len(myRow)]
    # 读取银经银维
    GLON = myRow[0]
    GLAT = myRow[1]
    
    # 银经0到90，270到360采用的算法
    if GLON<=90 or GLON>=270:
        # 获取最大的速度
        vmax = max(myRow[2:len(myRow)])
        # 计算距离和旋转速度
        R = R0*np.sin(GLON*np.pi/180)
        V = vmax + V0*np.sin(GLON*np.pi/180)
        R_list.append(R)
        V_list.append(V)

R_array = np.array(R_list)
V_array = np.array(V_list)

# 按距离排序，得到索引
Order = R_array.argsort()

# 利用索引对速度和距离进行排序
V_array = V_array[Order]
R_array = R_array[Order]

# 用5次多项式拟合，输出系数从高到低
Coe = np.polyfit(R_array, V_array, 5)
# 生成多项式
Pol = np.poly1d(Coe) 
# 得到结果
V_Smooth = Pol(R_array)

plt.rcParams['font.sans-serif']=['SimHei'];
plt.rcParams['axes.unicode_minus'] = False;

plt.xlim(0,10)
plt.ylim(0,300)
plt.title('银河旋转曲线')
plt.xlabel('到银心的距离[kpc]') 
plt.ylabel('旋转速度[km/s]') 
plt.scatter(R_array,V_array,marker='x',color='b')
plt.plot(R_array, V_Smooth,'r')
plt.plot()
plt.savefig('rocurve.png')