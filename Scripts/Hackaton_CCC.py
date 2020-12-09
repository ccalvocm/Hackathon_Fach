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
import rasterio.mask
from netCDF4 import Dataset
import numpy as np
import pandas as pd
import xarray as xr
import fiona
from numba import jit
from cdo import Cdo

cdo = Cdo()


@jit(nopython=True, nogil=True)
def distance(x1, y1, x2, y2):
    """
    Calculate distance from (x1,y1) to (x2,y2)
    """
    return ((x1-x2)**2 + (y1-y2)**2)**0.5

@jit(nopython=True, nogil=True)
def point_is_on_line(x, y, x1, y1, x2, y2):
    """
    Check whether point (x,y) is on line (x1,y1) to (x2,y2)
    """

    d1 = distance(x,  y,  x1, y1)
    d2 = distance(x,  y,  x2, y2)
    d3 = distance(x1, y1, x2, y2)

    eps = 1e-12
    return np.abs((d1+d2)-d3) < eps

@jit(nopython=True, nogil=True)
def is_left(xp, yp, x0, y0, x1, y1):
    """
    Check whether point (xp,yp) is left of line segment ((x0,y0) to (x1,y1))
    returns:  >0 if left of line, 0 if on line, <0 if right of line
    """

    return (x1-x0) * (yp-y0) - (xp-x0) * (y1-y0)

@jit(nopython=True, nogil=True)
def is_inside(xp, yp, x_set, y_set, size):
    """
    Given location (xp,yp) and set of line segments (x_set, y_set), determine
    whether (xp,yp) is inside polygon.
    """

    # First simple check on bounds
    if (xp < x_set.min() or xp > x_set.max() or yp < y_set.min() or yp > y_set.max()):
        return False

    wn = 0
    for i in range(size-1):

        # Second check: see if point exactly on line segment:
        if point_is_on_line(xp, yp, x_set[i], y_set[i], x_set[i+1], y_set[i+1]):
            return False

        # Calculate winding number
        if (y_set[i] <= yp):
            if (y_set[i+1] > yp):
                if (is_left(xp, yp, x_set[i], y_set[i], x_set[i+1], y_set[i+1]) > 0):
                    wn += 1
        else:
            if (y_set[i+1] <= yp):
                if (is_left(xp, yp, x_set[i], y_set[i], x_set[i+1], y_set[i+1]) < 0):
                    wn -= 1

    if wn == 0:
        return False
    else:
        return True

@jit(nopython=True, nogil=True)
def calc_mask(mask, lon, lat, shp_lon, shp_lat):
    """
    Calculate mask of grid points which are inside `shp_lon, shp_lat`
    """

    for j in range(lat.size):    
        for i in range(lon.size):
            if is_inside(lon[i], lat[j], shp_lon, shp_lat, shp_lon.size):
                mask[j,i] = True

#%% Acotar DEM de la cuenca

# Cargar shp de la cuenca del río Olivares

#nc_file = './Olivares_DEM_mask_1.nc'
#mascara
#mask_cuenca = Dataset(nc_file, mode='r')

#cob_glacial = 'ls8_sr_2013_2013_.nc'
                
                
os.chdir('./home/carlos/Downloads')
import rioxarray
import xarray
from shapely.geometry import mapping

test = xarray.open_dataset('ls8_sr_2013_2020_glaciar (1).nc')
MSWEP_monthly2 = test['nir_swir']
MSWEP_monthly2.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=True)
MSWEP_monthly2.rio.write_crs("epsg:32719", inplace=True)
Africa_Shape = geopandas.read_file('Rio_Olivares.shp', crs="epsg:32719")

clipped = MSWEP_monthly2.rio.clip(Africa_Shape.geometry.apply(mapping), Africa_Shape.crs, drop=False)


cob_glacial = 'ls8_sr_2013_2020_glaciar (1).nc'
glaciares = Dataset(cob_glacial, mode = 'r')
nirswir = glaciares['nir_swir']
tiempo = pd.read_csv('./time_2013_2020.csv',index_col = 0, parse_dates = True)
img_nubosas_sombra = np.zeros((len(tiempo.index)))
areas_glaciares = np.zeros((len(tiempo.index)))
   
############# Read shapefile and first feature
fc = fiona.open("Rio_Olivares.shp")
feature = next(iter(fc))

############# Extract array of lat/lon coordinates:
coords = feature['geometry']['coordinates'][0]
shp_lon = np.array(coords)[:,0]
shp_lat = np.array(coords)[:,1]

############# Read NetCDF variables, shifting the longitudes
nc_lon = glaciares['x'][:]
nc_lat = glaciares['y'][:][:]

for i in range(3): 
    nc_ua  = nirswir[i,:,:]
    nc_ua[nc_ua < 20] = np.nan

    ############# Calculate mask
    mask = np.zeros_like(nc_ua, dtype=bool)
    calc_mask(mask, nc_lon, nc_lat, shp_lon, shp_lat)
    
    ############# Mask the data array
    nc_ua_masked = np.ma.masked_where(~mask, nc_ua)
    
#    plt.pcolormesh(nc_lon, nc_lat, nc_ua_masked, vmin = 2, vmax = 30)
    ############# area pixeles dentro de la cuenca############# 
    #area_cuenca = 900*np.count_nonzero(mask)/1_000_000
    
    #numero pixeles glaciar dentro de la cuenca
    area_glaciar = np.count_nonzero(~np.isnan(nc_ua_masked))*900/1e6
    areas_glaciares[i] = area_glaciar
#plt.imshow(mask_cuenca['Band1'])
# cargar DEM ALOS-PALSAR
# path = '../Etapa 1 y 2/DEM/DEM Alos 5a a 8a mar.jp2'

# print(type(geom))
# with rasterio.open(path, 'r') as src:
#     out_image, out_transform = rasterio.mask.mask(src, geom, crop=True)
#     out_meta = src.meta

# out_meta.update({"driver": "GTiff",
#                  "height": out_image.shape[1],
#                  "width": out_image.shape[2],
#                  "transform": out_transform})

# with rasterio.open("Olivares_DEM.tif", "w", **out_meta) as dest:
#     dest.write(out_image)


# print bounding box and buffer in WGS84

basin = basin.to_crs('EPSG:4326')
bbox = basin.bounds
print(bbox)

center = [bbox['maxx'] + bbox['minx'], bbox['maxy'] + bbox['miny']]
print('centro: ', center[0]/2, center[1]/2)

buffer = [abs(bbox['maxx'] - bbox['minx']), abs(bbox['maxy'] - bbox['miny'])]
print(buffer)

#%% Raster plotting
import rasterio
from matplotlib import pyplot as plt

l8_2019_ro = rasterio.open('/home/carlos/Downloads/ls8_sr_2019.tiff')
#l8_2019_ro = rasterio.open('./scripts/ls8_sr_2019.tiff')
plt.show(l8_2019_ro)

# %%
