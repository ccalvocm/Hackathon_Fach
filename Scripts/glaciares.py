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
import ee

# load landsat 8 imagecollection using earth engine api
def load_landsat8_imagecollection(start_date, end_date):
    landsat8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterDate(start_date, end_date)
    return landsat8

def filter_landsat8_clouds(image):
    # get QA_PIXEL band from image
    qa = image.select('QA_PIXEL')
    # get cloud and cirrus bits
    cloud = 1 << 5
    cirrus = 1 << 6
    # get cloud and cirrus masks
    mask = qa.bitwiseAnd(cloud).eq(0).And(qa.bitwiseAnd(cirrus).eq(0))
    # return image with mask applied
    return image.updateMask(mask)

def filter_landsat8_clouds_bitwise(image):
    cloud_shadow_bit_mask = 1 << 3
    clouds_bit_mask = 1 << 5
    qa = image.select('QA_PIXEL')
    mask = qa.bitwiseAnd(cloud_shadow_bit_mask).eq(0).And(qa.bitwiseAnd(clouds_bit_mask).eq(0))

    # # If more than 40864234 pixels are masked, return everything masked
    masked_pixels = mask.reduceRegion(reducer=ee.Reducer.sum(), geometry=mask.geometry(), scale=30,
                                      maxPixels=1e13)
    # return 0 if masked_pixels QA_PIXEL is greater than 20% of the total pixels of the image
        
    return ee.Algorithms.If(ee.Number(masked_pixels.get('QA_PIXEL')).gt(408642.34),
                            image.updateMask(mask),ee.Image(0)
                            )

# calcular NIR/SWIR from landsat 8 surface reflectance
def calculate_nir_swir_landsat8(image):
    # get NIR from image
    nir = image.select('SR_B5')
    # adjust scale and offset
    nir = nir.multiply(2.75e-05).add(-0.2)
    # get SWIR from image
    # adjust scale and offset
    # mask swir where it is different from 0, else it is 1e-8
    swir = image.select('SR_B6').multiply(2.75e-05).add(-0.2).where(image.select('SR_B6').neq(0),1e-8)
    # return nir/swir renamed as nir_swir
    return nir.divide(swir).rename('nir_swir')

def filter_nir_swir(image):
    # retrun values greater than 3, else 0
    return image.gt(3).where(image.gt(3),0)

def mask_nir_swir(image):
    # mask values grater than 3 for further computations of area
    return image.updateMask(image.gt(3))

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
    glacier_shapefile = gpd.GeoDataFrame([],geometry=gpd.read_file(path).buffer(0))
    return glacier_shapefile

def calculate_area(image):
    # calculate sum of total area in EPSG 32719 of non zero pixels of image
    return image.reduceRegion(reducer=ee.Reducer.sum(),
                               geometry=image.geometry(),
                                 scale=30).getInfo()
def set_time_start(image):
    return image.set('system:time_start', image.get('system:time_start'))

# calculate ndvi
def calculate_area(image):
    # calculate sum of total area in EPSG 32719 of non zero pixels of image
    return ee.Image(image.reduceRegion(reducer=ee.Reducer.sum(),
                                geometry=glacier_fc.geometry(),
                                scale=30).getInfo())
def clip_image(image):
    return image.clip(glacier_fc)

def main():
    ee.Authenticate()
    ee.Initialize()

    # load glacier shapefile
    glacier_shapefile = load_glacier_shapefile().to_crs(epsg='4326')

    # convert to featurecollection
    glacier_fc = gdf2FeatureCollection(glacier_shapefile)

    # load landsat 8 imagecollection and filterbounds
    yrs=['2013','2014','2015','2016','2017','2018','2019','2020','2021',
    '2022','2023']

    df=pd.DataFrame(index=yrs,columns=['area'])

    for yr in yrs:
        landsat8 = load_landsat8_imagecollection(yr+'-12-01', 
        str(int(yr)+1)+'-03-31')

        # Filter out images with missing values
        landsat8 = landsat8.filterBounds(glacier_fc).map(clip_image)\
            .map(set_time_start).map(filter_landsat8_clouds)

        # Check if bands SR_B5 and SR_B6 exist in the Landsat 8 image
        nir_swir = landsat8.filter(ee.Filter.listContains('system:band_names', 
        'SR_B5')).filter(ee.Filter.listContains('system:band_names', 
        'SR_B6')).map(calculate_nir_swir_landsat8).map(mask_nir_swir)

        # Calculate the mean of NIR_SWIR
        nir_swir_mean = nir_swir.mean()

        maskedArea = ee.Image.pixelArea().mask(nir_swir_mean.gt(3))
        area = maskedArea.reduceRegion(reducer=ee.Reducer.sum(), geometry=glacier_fc.geometry(), scale=30, maxPixels=1e13).get('area').getInfo()
        print(area)
        df.loc[yr, 'area'] = area


    # def calculate_yearly_mean(image_collection):
    #     def calculate_mean(year):
    #         start_date = ee.Date.fromYMD(year, 1, 1)
    #         end_date = ee.Date.fromYMD(year, 12, 31)
    #         yearly_images = image_collection.filterDate(start_date, end_date)
    #         return yearly_images.mean().set('year', year)
        
    #     years = ee.List.sequence(2019, 2021)  # Update the range of years if needed
    #     yearly_means = years.map(calculate_mean)
        
    #     return ee.ImageCollection.fromImages(yearly_means)

    # Filter out images with missing values
    nir_swir_filtered = nir_swir.filter(ee.Filter.notNull(['nir_swir']))

    # # Calculate yearly mean of filtered image collection
    # nir_swir_yearly_mean = calculate_yearly_mean(nir_swir_filtered)




    # Filter out images with missing values
    # nir_swir_filtered = nir_swir_yearly_mean.filter(ee.Filter.notNull(['nir_swir']))

    # Compute area for each non-masked pixel in filtered image collection nir_swir
    # area_list = nir_swir_filtered.map(calculate_area)

    # Convert to pandas dataframe
    # df = pd.DataFrame.from_records(area_list.getInfo())

    # # Print the dataframe
    # print(df)

    map = geemap.Map()
    map.addLayer(nir_swir_mean)
    map

    # area_list = nir_swir.map(calculate_area)

    # print max value for nir_swir.mean()
    print(nir_swir.mean().reduceRegion(reducer=ee.Reducer.mean(), 
                                        scale=30,
                                        maxPixels=1e13,
                                        geometry=glacier_fc).getInfo())
    map=geemap.Map()
    map.addLayer(nir_swir.mean())
    map

