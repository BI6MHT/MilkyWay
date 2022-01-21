# -*- coding: utf-8 -*-
"""
Created on Sat Jan 15 13:58:29 2022

@author: BI6MHT
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cmath

# 太阳距离银心的距离和旋转速度
R0 = 8.5 #kpc
V0 = 220 #km/s

# 读取文件
table = pd.read_csv('MilkyWay.csv')

#返回行数
rows_num = table.shape[0]
#print(rows_num)
# Outfile

R_list = []
V_list = []
fig = plt.figure()
ax = fig.add_subplot(111)

for i in range(0,rows_num):
    # 读取第i行的内容
    myRow = table.iloc[i,:]
    #print(myRow)
    # 去掉赤经赤纬
    myRow = myRow[2:len(myRow)]
    #print(myRow)
    # 读取银经银维
    GLON = myRow[0]
    GLAT = myRow[1]
    
    # 除去最后一行的最大径向速度
    vels = myRow[2:len(myRow)-1]
    #print(vels)
    theta = GLON -90 + 180;
    for VR in vels:
        
        # 如果为NaN或者None，则跳过本次循环
        if pd.isnull(VR):
            continue
        
        # 求出峰值速度对应的到银心的距离
        R = R0*V0*np.sin(GLON*np.pi/180)/(V0*np.sin(GLON*np.pi/180)+VR)
       
        # 计算分子云的x,y坐标轴
        rp = cmath.sqrt(R**2-R0**2*(np.sin(GLON*np.pi/180))**2)+R0*np.cos(GLON*np.pi/180)
        rm = -cmath.sqrt(R**2-R0**2*(np.sin(GLON*np.pi/180))**2)+R0*np.cos(GLON*np.pi/180)
        
        
        plt.rcParams['font.sans-serif']=['SimHei'];
        plt.rcParams['axes.unicode_minus'] = False;
        
        # 如果有两者都为实数，显然有rm<rp，故只有一根的情况下，肯定是rm.real<0
        # rp.real<0，则直接没有根了
        if rp.real>0:
               
            # rp.real>0的时候，rm.real<0，则说明两个都为实数
            # 因为负数开平方根只能开出来虚数，对实部的正负无影响
            # 故如果是复数的话，两者的实部部分正负是相同的。
            if rm.real<0:
                r = rp.real
                x = r*np.cos(theta*np.pi/180)
                y = -R0 + r*np.sin(theta*np.pi/180)
                # 作图
                plt.scatter(x,y,marker='x',c='black')
            
            # rp.real>0，且rm.real>=0，两者又都为实数
            # 故有两个解存在
            elif rp.imag == 0 and  rm.imag == 0:
                # 发现二解
                # 处理这些的最简单的方法是把它们都画出来，但是用另一种颜色 '
                r = rp.real;
                x = r*np.cos(theta*np.pi/180)
                y = -R0+ r*np.sin(theta*np.pi/180)
                # 画图
                plt.scatter(x,y,marker='x',c='r')
                r = rm.real
                x = r*np.cos(theta*np.pi/180)
                y = -R0 + r*np.sin(theta*np.pi/180)
                # 画图
                plt.scatter(x,y,marker='x',c='b')
        #print([GLON,VR,rp,rm])

# 画出银河和太阳的位置
plt.text(0, 0, 'C')

plt.text(0, -8.5, 'Sun')

plt.text(-15, 15, 'Q1')

plt.text(-15, -15, 'Q2')

plt.text(15, -15, 'Q3')

plt.text(15, 15, 'Q4')

# 设置画图极限，即0到10kpc，0到300km/s


plt.xlim(-25,25)
plt.ylim(-25,25)
# 设置为正方形
ax.set_aspect('equal', adjustable='box')
plt.title('银河旋臂图')
plt.xlabel('[kpc]') 
plt.ylabel('[kpc]') 
plt.savefig('map.png')