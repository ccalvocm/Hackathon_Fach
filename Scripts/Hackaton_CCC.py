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


#%% Acotar Xarray de la cuenca
       
os.chdir('/home/carlos/Downloads')
import rioxarray
import xarray
from shapely.geometry import mapping

cdo = Cdo()

cob_glacial = 'ls8_sr_2013_2020_glaciar (1).nc'
glaciares = xarray.open_dataset(cob_glacial)
nirswir = glaciares['nir_swir']
nirswir.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=True)
nirswir.rio.write_crs("epsg:32719", inplace=True)

############# Read shapefile and first feature
Cuenca_Shape = geopandas.read_file('Rio_Olivares.shp', crs="epsg:32719")


tiempo = pd.read_csv('./time_2013_2020.csv',index_col = 0, parse_dates = True)
img_nubosas_sombra = np.zeros((len(tiempo.index)))
areas_glaciares = np.zeros((len(tiempo.index)))
   
############# Read NetCDF variables

for i in range(len(tiempo)): 
    nc_ua  = nirswir[i,:,:]
    nc_ua = nc_ua.where(nc_ua >= 2, other=np.nan) 

    ############# Calculate mask
    
    clipped = nc_ua.rio.clip(Cuenca_Shape.geometry.apply(mapping), Cuenca_Shape.crs, drop=False)
    
    #numero pixeles glaciar dentro de la cuenca
    area_glaciar = np.count_nonzero(~np.isnan(clipped))*900/1e6
    areas_glaciares[i] = area_glaciar

