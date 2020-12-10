#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 23:00:30 2020

@author: felipe
"""

# importar librerias
import os
import geopandas
from matplotlib import pyplot as plt
import rasterio
import numpy as np
import pandas as pd
from cdo import Cdo
import gc
import cv2


def leap_year(year, calendar='standard'):
    """Determine if year is a leap year"""
    leap = False
    if ((calendar in ['standard', 'gregorian',
        'proleptic_gregorian', 'julian']) and
        (year % 4 == 0)):
        leap = True
        if ((calendar == 'proleptic_gregorian') and
            (year % 100 == 0) and
            (year % 400 != 0)):
            leap = False
        elif ((calendar in ['standard', 'gregorian']) and
                 (year % 100 == 0) and (year % 400 != 0) and
                 (year < 1583)):
            leap = False
    return leap

dpm = {'noleap': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       '365_day': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'standard': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'gregorian': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'proleptic_gregorian': [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       'all_leap': [0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       '366_day': [0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       '360_day': [0, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30]}

def get_dpm(time, calendar='standard'):
    """
    return a array of days per month corresponding to the months provided in `months`
    """
    month_length = np.zeros(len(time), dtype=np.int)

    cal_days = dpm[calendar]

    for i, (month, year) in enumerate(zip(time.month, time.year)):
        month_length[i] = cal_days[month]
        if leap_year(year, calendar=calendar):
            month_length[i] += 1
    return month_length

#%% Acotar Xarray de la cuenca
       
os.chdir('/home/carlos/Downloads')
import rioxarray
import xarray
from shapely.geometry import mapping

cdo = Cdo()

cob_glacial = 'ls8_sr_2013_2020_glaciar (1).nc'
cob_glacial = 'ls8_toa_2013_2020_glaciar.nc'
glaciares = xarray.open_dataset(cob_glacial)
nirswir = glaciares['nir_swir']
swir1 = glaciares['swir1']
#coordenadas
nirswir.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=True)
nirswir.rio.write_crs("epsg:32719", inplace=True)

############# Read shapefile and first feature
Cuenca_Shape = geopandas.read_file('Rio_Olivares.shp', crs="epsg:32719")
catastro_DGA = geopandas.read_file('glaciares.shp', crs="epsg:32719")

tiempo = pd.read_csv('./time_2013_2020.csv',index_col = 0, parse_dates = True)
tiempo = pd.read_csv('./time_2013_2020_toa.csv',index_col = 0, parse_dates = True)
img_nubosas_sombra = np.zeros((len(tiempo.index)))
areas_glaciares = np.zeros((len(tiempo.index)))
flags = np.zeros((len(tiempo.index)))
   
############# Read NetCDF variables
default = 2
dicc = {1 : default, 2: default, 3: default ,4: default ,5: default , 6: 2.2, 7: 2.1, 8: 3, 
        9: 1.9, 10: default, 11: 1.9, 12: 1.5}
#for i in range(round(len(tiempo[0:41]))): 
for i in range(round(len(tiempo))): 


    nc_ua  = nirswir[i,:,:]

    a = ndimage.interpolation.zoom(nc_ua,.5)  # (20, 40)


#    swir1_flag = swir1[i,:,:].values[swir1[i,:,:].values == -9999]
#    flags[i] = len(swir1_flag)
#    if len(swir1_flag) > 1e2:
#        areas_glaciares[i] = np.nan
#        continue
#    else:

    thres = dicc[tiempo.index[i].month]
#        thres = 2
    nc_ua = nc_ua.where(((nc_ua >= thres) & (nc_ua <= 3e10)), other=np.nan) 

############# Calculate mask

#        clipped = nc_ua.rio.clip(Cuenca_Shape.geometry.apply(mapping), Cuenca_Shape.crs, drop=False)
    clipped = nc_ua.rio.clip(catastro_DGA.geometry.apply(mapping), catastro_DGA.crs, drop=False)

#numero pixeles glaciar dentro de la cuenca
    area_glaciar = np.count_nonzero(~np.isnan(clipped))*900/1e6
    areas_glaciares[i] = area_glaciar
    gc.collect()
    del nc_ua
        
plt.close("all")
plot_areas = pd.DataFrame(rolling_windows.mean())
plot_areas = plot_areas.resample('YS').mean()
rolling_windows = plot_areas.rolling(3, min_periods=1)
plot_areas_mm = rolling_windows.mean().plot()

plt.plot(tiempo.index.strftime("%b %Y"), areas_glaciares)
plt.xticks(rotation=90)
plt.grid()
