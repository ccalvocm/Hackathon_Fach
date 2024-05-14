import os
import osgeo
import ee
import pandas as pd
import datetime
import requests
import zipfile
import xarray as xr
import rioxarray as rxr
import geopandas as gpd
import numpy as np
import json
import matplotlib.pyplot as plt
import geemap

# load landsat 8 imagecollection using earth engine api
def load_landsat8_imagecollection(start_date, end_date):
    landsat8 = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR').filterDate(start_date, end_date)
    return landsat8

def filter_landsat8_clouds_bitwise(image):
    cloud_shadow_bit_mask = 1 << 3
    clouds_bit_mask = 1 << 5
    qa = image.select('pixel_qa')
    mask = qa.bitwiseAnd(cloud_shadow_bit_mask).eq(0).And(qa.bitwiseAnd(clouds_bit_mask).eq(0))
    return image.updateMask(mask)

def calculate_nir(image):
    return image.normalizedDifference(['B5', 'B4'])

def calcularte_ndgi_landsat8(image):
    return image.normalizedDifference(['B3', 'B5'])

# calcular NIR/SWIR from landsat 8 surface reflectance
def calculate_nir_swir_landsat8(image):
    # get NIR from image
    nir = image.select('B5')
    # adjust scale and offset
    nir = nir.multiply(2.75e-05).add(-0.2)
    # get SWIR from image
    # adjust scale and offset
    # mask swir where it is different from 0, else it is 1e-8
    swir = image.select('B6').multiply(2.75e-05).add(-0.2).where(image.select('B6').neq(0),1e-8)
    return nir.divide(swir)

def filter_nir_swir(image):
    # retrun values greater than 3, else 0
    return image.gt(3).where(image.gt(3),0)

def fixMultipoly(geo):
    if 'MultiPolygon' in geo.geometry.iloc[0].geom_type:
        geo2=geo.explode()
        geo3=gpd.GeoDataFrame([],geometry=geo2.geometry,crs='4326')
        geo3['area']=''
        geo3['area']=geo3.to_crs(epsg='32719').apply(lambda x: x['geometry'].area,
                                                        axis=1)
        geo3=geo3.sort_values('area',ascending=False)
        geo3=geo3.iloc[0].drop('area')
        geo3=gpd.GeoDataFrame(pd.DataFrame(geo3).T)
        return geo3.buffer(0)
    else:
        return geo.buffer(0)

def gdf2FeatureCollection(gs):
    features = []
    for i in range(gs.shape[0]):
        geom = gs.iloc[i:i+1,:] 
        geom=fixMultipoly(geom)
        jsonDict = json.loads(geom.to_json())
        x=np.array([x[0] for x in jsonDict['features'][0]['geometry']['coordinates'][0]])
        y=np.array([x[1] for x in jsonDict['features'][0]['geometry']['coordinates'][0]])
        cords = np.dstack((x[:],y[:])).tolist()
        # g=ee.Geometry.Polygon(cords).bounds()
        g=ee.Geometry.Polygon(cords)
        # feature = ee.Feature(g,{'name':self.gdf.loc[i,col].astype(str)})
        feature = ee.Feature(g)
        features.append(feature)
    return ee.FeatureCollection(features)

def load_glacier_shapefile(path=os.path.join('..','geoData','glaciaresOlivares.shp')):
    glacier_shapefile = gpd.read_file(path)
    return glacier_shapefile

def calculate_area(image):
    # calculate sum of total area in EPSG 32719 of non zero pixels of image
    return image.mask(image).reduceRegion(reducer=ee.Reducer.sum(),
                                          geometry=image.geometry(),
                                          scale=30)


def main():
    ee.Authenticate()
    ee.Initialize()

    # load glacier shapefile
    glacier_shapefile = load_glacier_shapefile().to_crs(epsg='4326')

    # convert to featurecollection
    glacier_fc = gdf2FeatureCollection(glacier_shapefile)

    # load landsat 8 imagecollection and filterbounds
    landsat8 = load_landsat8_imagecollection('2015-01-01', '2024-05-02')

    # filter bounds
    landsat8 = landsat8.filterBounds(glacier_fc).map(lambda image: image\
    .clip(glacier_fc)).map(lambda image: image.set('system:time_start', \
    image.get('system:time_start'))).map(lambda image: filter_landsat8_clouds_bitwise(image))

    # calculate ndvi
    nir_swir = landsat8.map(lambda image: calculate_nir_swir_landsat8(image))\
        .map(lambda image: filter_nir_swir(image))

    # calculate area by date
    area = nir_swir.map(calculate_area)

    # plot first image from nir_swir using geemap
    nir_swir_first = nir_swir.mean()
    # plot using geemap
    Map = geemap.Map()
    Map.addLayer(nir_swir_first, {'min': 4, 'max': 1000, 'palette': ['blue', 'white', 'green']})
    Map




