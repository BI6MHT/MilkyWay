# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 01:13:02 2022

@author: BI6MHT
"""
from astropy.coordinates import SkyCoord 
import astropy.units  as u

# 无关文件
# 获取银纬0°时，每隔5°对应赤经赤纬
# 0°和360°是一样的

for i in range(0,91,5):
    
    Coord = str(i)+'d'+' '+'0d'
    target = SkyCoord(Coord, frame="galactic",unit = u.deg)
    Ra = target.icrs.ra.to_string(u.hour)
    Dec = target.icrs.dec.to_string(u.deg)
    print(str(i)+' '+Ra+' '+Dec)
    
for i in range(270,360,5):
    Coord = str(i)+'d'+' '+'0d'
    target = SkyCoord(Coord, frame="galactic", unit = u.deg)
    Ra = target.icrs.ra.to_string(u.hour)
    Dec = target.icrs.dec.to_string(u.deg)
    print(str(i)+' '+Ra+' '+Dec)
    