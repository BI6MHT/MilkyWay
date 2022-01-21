# -*- coding: utf-8 -*-
"""
Created on Mon Jan 17 09:19:39 2022

@author: BI6MHT
"""

from astropy.coordinates import SkyCoord,Angle,EarthLocation
from astropy.time import Time
from astropy.constants import c,k_B
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt

# 对曲线进行平滑
#import statsmodels.api as sm
#lowess = sm.nonparametric.lowess

# 寻求极大值；对曲线进行平滑
from scipy import signal

# 单个天体的获取

class MilkyWay():
    

    """
    初始化类的属性
    
    obs_Dir:观测文件所在根目录
    
    obs_Site:观测者在地球上所在位置，例如
    
    obs_Source:观测源的赤经赤纬，例如'12h10m0s 20d'
    
    obs_Time:观测的UTC时间，例如'2021-10-30 13:00:43'
    
    smooth_Window(默认值31): 即Savitzky-Golay滤波器的窗口长度，值越小，曲线越贴近真实曲线；值越大，平滑效果越厉害（备注：该值必须为正奇整数）
    
    smooth_k(默认值3): 值越大，曲线越贴近真实曲线；值越小，曲线平滑越厉害。另外，当k值较大时，受窗口长度限制，拟合会出现问题，高频曲线会变成直线。
    
    """
    def __init__(self,obs_Dir,obs_Site,obs_Source,obs_Time,smooth_Window=21,smooth_k=3):
        
        self.obs_Freq = 1420.4057517667
        
        self.obs_Dir = obs_Dir
        self.obs_Site = obs_Site
        self.obs_Source = obs_Source
        self.obs_Time = obs_Time
        self.smooth_Window = smooth_Window
        self.smooth_k = smooth_k
        self.obs_GalLon,self.obs_GalLat = self.get_GalCoord()
        
        # 获取频谱参数，第一列为频率(MHz)；第二列是强度，默认单位为W
        self.spec = np.loadtxt(self.obs_Dir+self.obs_Source+'.txt',dtype=float,skiprows=1)
        
        # 取-150到150km/s之间的径向速度rv，以及对应的等效温度(平滑和未平滑)
        self.rv, self.equiTemp_Smooth,self.equiTemp = self.get_Smooth()
        
        # 得到峰值点以最大径向速度那一点
        # 按峰值降序排列
        # 最大速度那一点对应末端元素
        self.rv_Array,self.T_Array = self.get_Peaks()
    
    """
    获取源所在的银经和银维
    返回类型：list类型
    返回结果：银经 银维
    返回例子：'54.2663 0.224174'
    
    """
    def get_GalCoord(self):
        
        Source = SkyCoord(self.obs_Source,frame='icrs')
        Gal = Source.galactic.to_string()
        Gal_lon = float(Gal.split(' ')[0])
        Gal_lat = float(Gal.split(' ')[1])
        return Gal_lon,Gal_lat
            
    """
    根据多普勒效应，获取不同频率对应的修正到日心的径向速度
    返回类型：一维的numpy.ndarray
    返回结果：不同频率对应的径向速度
    返回例子：[1098.31559636 ... -1005.95448058 -1008.01546988]
    """
    def get_CorrectVel(self):
        
        obs_Time = Time(self.obs_Time)
        
        obs_SiteLon = Angle(self.obs_Site.split(' ')[0],unit=u.deg)
        obs_SiteLat = Angle(self.obs_Site.split(' ')[1],unit=u.deg)
        obs_SiteHei = float(self.obs_Site.split(' ')[2])*u.m
        loc = EarthLocation.from_geodetic(obs_SiteLon,obs_SiteLat,obs_SiteHei)
        
        # 观测天体源所在的赤经赤纬，坐标框架采用icrs，即天球坐标系坐标
        Source = SkyCoord(self.obs_Source,frame='icrs')
        
        # 读取频谱文件，获取观测频点
        
        spec_freq = self.spec[:,0]
        
        # 地球上观测者和中性氢在视线上的相对速度，即径向速度
        # 根据多普勒效应，不同频点对应的不同径向速度
        rv = (self.obs_Freq-spec_freq)*c/self.obs_Freq
         
        # 径向速度vr是观测者和中性氢在视线上的相对速度
        # 现在我们要求得以太阳系质心为原点时（或者大概认为是在太阳中心观测），太阳中心和中性氢在两者连线上的相对速度
        # 即进行坐标转换，因而要求出一个速度修正项
        vcorr = Source.radial_velocity_correction(kind='barycentric',obstime=obs_Time, location=loc)
        
        # 加上速度修正项以后，便得到了在太阳上观测中性氢的径向速度
        # m/s转成km/s,然后取其数值
        rv = (rv + vcorr + rv*vcorr/c)
        rv = rv.to(u.km/u.s).to_value()
        
        # 显然，频率较大(蓝移)时，rv为负；频率较小（红移）时，rv为正
        # 所以目前rv是的排序是由正到负的，不妨让rv逆序排列
        rv = rv[::-1]
        
        return rv
    
    """
    获取频率强度(电平)对应的等效温度
    返回类型：一维的numpy.ndarray
    返回结果：不同径向速度（频率）对应的等效温度
    返回例子：[1.31519478 1.26451874 ... 1.40664323 1.37319388 1.22032352]
    """
    def get_EquiTemp(self):
       
        spec_freq = self.spec[:,0]
        spec_level = self.spec[:,1]
        
        # 观测带宽，由于MHz，故乘以10^6
        spec_width = (spec_freq[-1]-spec_freq[0])*np.power(10,6)
       
        # 根据奈奎斯特噪声或者说热噪声，有P=kBT，其中P为功率，k为玻尔兹曼常数，B为带宽，T为温度
        # 即我们可以将功率等效为温度
        # 由于等效出来的数量级过大，故可乘以10^(-27)以减小数量级
        equiTemp = (spec_level/k_B/spec_width)*pow(10,-14)
        
        # 由于径向速度rv是逆序排列的，故另EquiTemp也逆向排序，以和rv一一对应
        equiTemp = equiTemp[::-1]
        
        # 由于等式中出现从astropy中调用的k_B常量，故equiTemp是 <class 'astropy.units.quantity.Quantity'>
        # 取其数值
        equiTemp = equiTemp.to_value()

        return equiTemp
    
    """
    获取-150到150km/s之间的径向速度rv以及对应的等效温度T，并进行平滑处理
    返回类型：一维的numpy.ndarray，一维的numpy.ndarray
    返回结果：径向速度，等效温度
    返回例子：[-148.58293111 ... 146.13853902],[1.31519478 ... 1.22032352]
    """
    def get_VT(self):
        
        rv = self.get_CorrectVel()
        equiTemp = self.get_EquiTemp()
        
        rv_index=np.where((rv>-150)&(rv<150))[0]
        rv_start = rv_index[0]
        rv_end = rv_index[-1]
        rv_150 = rv[rv_start:rv_end]
        equiTemp_150 = equiTemp[rv_start:rv_end]

        return rv_150,equiTemp_150
    

    # def get_Smooth(self):
    #     rv, equiTemp = self.get_VT()
    #     # 对数值进行平滑化
    #     # 其中的1/20代表平滑程度，可自行调整
    #     T_Smooth = lowess(equiTemp,rv,frac=self.Smooth)[:,1]
    #     return rv,T_Smooth,equiTemp
    
     
    """
    对 径向速度rv vs 等效温度equiTemp 的曲线进行平滑处理
    其中vr在-150km/s到150km/s之间
    返回类型：一维的numpy.ndarray，一维的numpy.ndarray，一维的numpy.ndarray
    返回结果：径向速度，平滑等效温度，未平滑等效温度
    返回例子：[-148.58293111 ... 146.13853902],[1.31519478 ... 1.22032352],[1.31519478 ... 1.22032352]
    """
    def get_Smooth(self):
        rv, equiTemp = self.get_VT()
        # 对数值进行平滑
        
        T_Smooth = signal.savgol_filter(equiTemp,self.smooth_Window,self.smooth_k)
        return rv,T_Smooth,equiTemp
    
    """
    获取频谱图
    """
    def get_Plot(self):
        plt.figure()        
        
        #设置使用的字体为支持中文的字体
        plt.rcParams['font.sans-serif']=['SimHei'];
        plt.rcParams['axes.unicode_minus'] = False;
        
        # 未平滑
        plt.subplot(2,1,1)
        plt.plot(self.rv,self.equiTemp)
        plt.title(self.obs_Source + ' ' + self.obs_Time)
        
        # 平滑
        plt.subplot(2,1,2)
        plt.plot(self.rv,self.equiTemp_Smooth);
        plt.xlabel('中性氢分子团相对于日心的径向速度vr(km/s)')
        plt.ylabel('等效温度T(K)')
        
        # 画出峰值，以及最大径向速度那一点
        for i in range(0,len(self.T_Array)):
            plt.scatter(self.rv_Array[i],self.T_Array[i],s=25,c='r') 
            plt.text(self.rv_Array[i],self.T_Array[i],\
                 str(np.around(self.rv_Array[i],3))+','+str(np.around(self.T_Array[i],3)),\
                 fontdict={'fontsize':8});
        
        # 以银经命名，银经保留一位小数
        plt.savefig('Images/'+str(np.around(self.obs_GalLon,1))+'.png')
        

        plt.clf() # 清图
        plt.cla() # 清坐标轴
        plt.close() # 关窗口
    
    """
    寻求最大峰以及右边的峰，即求极值问题
    返回类型：一维的numpy.ndarray，一维的numpy.ndarray
    返回结果：峰的径向速度以及最大径向速度，对应的等效温度
    返回例子：[-148.58293111 ... 146.13853902],[1.31519478 ... 1.22032352]
    """
    def get_Peaks(self):
        
        # 观察频谱规律后，对于银经0到90°，270°到360°，寻峰算法如下
        
        # 求取最高峰所在位置
        max_Index = np.argmax(self.equiTemp_Smooth)
        # 峰值间的最小间隔点数，两峰之间的间隔需大于此
        # 此处设两峰之间的间隔需大于30km/s
        lim_Distance = np.ceil((30*len(self.equiTemp_Smooth)/300))
        # 寻峰，即求极值点
        # 此处只求最大峰及其左边的峰
        peaks =signal.find_peaks(self.equiTemp_Smooth[0:max_Index+2],distance = lim_Distance)[0]
        
        # 找出极大值对应的径向速度和等效温度
        T_Array = np.array([])
        rv_Array = np.array([])
        for peak in peaks:
            T_Array= np.append(T_Array,self.equiTemp_Smooth[peak])
            rv_Array = np.append(rv_Array,self.rv[peak])
        
        # 按峰值降序排序，得到索引
        Order = T_Array.argsort()[::-1]
        # 利用索引对径向速度和等效温度进行排序
        rv_Array = rv_Array[Order]
        T_Array = T_Array[Order]
        
        # 在尾端加入最大径向速度
        mark_Vel,mark_T = self.get_MaxVel()
        
        rv_Array = np.append(rv_Array,mark_Vel)
        T_Array = np.append(T_Array,mark_T)
        
        
        return rv_Array,T_Array
    
    """
    获取氢谱线峰的最大径向速度，以及对应的等效温度
    峰是弥散开来，可以求右边最高峰衰减到一定值处对应的速度
    返回类型： float
    返回结果：最大径向速度，等效温度
    """
    def get_MaxVel(self):
        
        # 此处以(右边最高峰)衰减到(右边最高峰)的右边频谱积分平值为截止标准
        max_Index = np.argmax(self.equiTemp_Smooth)
        mark_Index = self.get_RightFirst(self.equiTemp_Smooth,max_Index)
        mark_Vel = self.rv[mark_Index]
        mark_T =  self.equiTemp_Smooth[mark_Index]
        
        return mark_Vel,mark_T
    
    """
    返回第一个小于右边频谱积分平值的索引
    返回类型： int
    返回结果：索引
    """

    def get_RightFirst(self,Arr,Index):
        
        # 边缘峰的右边频谱的积分平均值
        # 不能用[Index:-1]，可以实际测试一下
        spec_Ave = np.sum(Arr[Index:len(Arr)])/len(Arr[Index:len(Arr)])
        # 寻求第一个小于spec_Ave对应的索引
        for i in range(Index,len(Arr)):
            i+1
            if Arr[i] <= spec_Ave:
                break
        return i
    