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
from scipy.ndimage import gaussian_filter


#%% Acotar Xarray de la cuenca
       
os.chdir('/home/carlos/Downloads')
import rioxarray
import xarray
from shapely.geometry import mapping

#%%
def main():
    #%%
    cdo = Cdo()
    
    cob_glacial = 'ls8_sr_2013_2020_glaciar (1).nc'
    cob_glacial = 'ls8_toa_2013_2020_glaciar.nc'
    glaciares = xarray.open_dataset(cob_glacial)
    nirswir = glaciares['nir_swir']
    times = pd.read_csv('./time_2013_2020_toa.csv',index_col = 0, parse_dates = True).index
    # assign times
    
    nirswir = nirswir.assign_coords(time = times)
        
    #coordenadas
    nirswir.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=True)
    nirswir.rio.write_crs("epsg:32719", inplace=True)
    
    ############# Read shapefile and first feature
    cobertura_nival = geopandas.read_file('glaciares.shp', crs="epsg:32719")
    
#    img_nubosas_sombra = np.zeros((len(tiempo.index)))
    areas_glaciares = np.zeros((len(times)))
#    flags = np.zeros((len(tiempo.index)))
       
    ############# Read NetCDF variables
    default = 2
    dicc = {1 : default, 2: default, 3: default ,4: default ,5: default , 6: 2.2, 7: 2.1, 8: 3, 
            9: 1.9, 10: default, 11: 1.9, 12: 1.5}
    #for i in range(round(len(tiempo[0:41]))): 
    for i in range(round(len(times))): 
    
    
        nc_ua  = nirswir[i,:,:]
            
#        result = gaussian_filter(nc_ua, sigma=2)
#        nc_ua.values = result  
        
        thres = dicc[times[i].month]
        thres = 5
        nc_ua = nc_ua.where(((nc_ua >= thres)), other=np.nan) 
        
    ############# Calculate mask
    
        clipped = nc_ua.rio.clip(cobertura_nival.geometry.apply(mapping), cobertura_nival.crs, drop=False)
    
    #numero pixeles glaciar dentro de la cuenca
        area_glaciar = np.count_nonzero(~np.isnan(clipped))*900/1e6
        areas_glaciares[i] = area_glaciar
        gc.collect()
        del nc_ua
    #    del a
  #%%          
    plt.close("all")
    plot_areas = pd.DataFrame(areas_glaciares, index = times)
    plot_areas.columns = ['Area']
    plot_areas = plot_areas[plot_areas['Area'] > 20]
    ax = plot_areas.rolling(15).mean().plot(rot = 90)
#     xlabel = plot_areas.index.strftime("%b %Y").to_list(), rot = 90
#    ax.set_xticks(plot_areas.index, times.strftime("%b %Y").to_list(), rotation = 90)
#    rolling_windows = plot_areas.rolling(15, min_periods=1)
#    plot_areas_mm = rolling_windows.mean().plot()
    ax.set_ylabel('Cobertura glaciar ($km^2$)')
    ax.set_xlabel('')
    ax.set_ylim(bottom =0)
    ax.grid()