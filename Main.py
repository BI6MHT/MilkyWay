# -*- coding: utf-8 -*-
"""
Created on Mon Jan 17 16:22:39 2022

@author: BI6MHT
"""

# 批量处理

import os
import datetime
import pandas as pd
from MilkyWay import MilkyWay


obs_Dir = 'Data/'
obs_Site = '110.81285 32.59175 40'

# 读取目录下的文件
obs_Sources = os.listdir(obs_Dir)


"""
根据文件创建时间返回UTC时间
具体的是以修改时间还是创建时间为观测时间，还是要看看文件的时间属性的
如果要以修改时间为观测时间，可把LocalTime = os.stat(obs_Dir+obs_Source+'.txt').st_ctime中的st_ctime改为st_mtime
"""
def get_ObsTime(obs_Source):
    # 获取文件创建时间，类型为float
    LocalTime = os.stat(obs_Dir+obs_Source+'.txt').st_ctime
    # 获取文件的UTC时间，类为'datetime.datetime'
    UTCTime = datetime.datetime.utcfromtimestamp(LocalTime)
    # 将UTC时间转为字符串
    obs_Time = UTCTime.strftime("%Y-%m-%d %H:%M:%S")
    
    return obs_Time

# 用于列表排序，即根据第2个元素进行排序
def takeSecond(elem):
        return elem[1]

List = []

for i in obs_Sources:
    
    # 去除txt后缀的文件名，得到观测目标的天球坐标
    obs_Source = i.split('.',1)[0]
    
    # 得到观测时间
    obs_Time = get_ObsTime(obs_Source)
    Hydro = MilkyWay(obs_Dir,obs_Site,obs_Source,obs_Time)
    
    # 只画出银经0°到90°，270°到360°的频谱
    # if 0<=Hydro.obs_GalLon<=90 or 270<=Hydro.obs_GalLon<=360:
    Hydro.get_Plot()
    child_List = [Hydro.obs_Source,Hydro.obs_GalLon,Hydro.obs_GalLat]
    rv_Array = Hydro.rv_Array.tolist()
    child_List = child_List+rv_Array
    List.append(child_List)

# 根据递增的银河经度进行排序
List = sorted(List,key=takeSecond)

# 储存到excel表格中
# 注意，如果再次运行储存csv操作，请确保上次的MilkyWay文件没有被打开
List_Excel = pd.DataFrame(List)
List_Excel.to_csv('MilkyWay.csv')
