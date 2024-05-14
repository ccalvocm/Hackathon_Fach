import os
import ee
import pandas as pd
import datetime
import requests
import zipfile
import xarray as xr
import osgeo
import rioxarray as rxr
import geopandas as gpd
import numpy as np
import json

# def time_index_from_filenames(filenames,first=14,last=4):
def time_index_from_filenames(filenames,first=41,last=33):
    '''helper function to create a pandas DatetimeIndex
       Filename example: 20150520_0164.tif'''
    # return pandas DatetimeIndex from a list of filenames with dates in format YYYYMMDD
    return pd.DatetimeIndex([pd.Timestamp(f[-first:-last]) for f in filenames])

def files2NETCDF(folderName,band='sr_',format='.tif'):
    filenames = os.listdir(folderName)
    files = sorted([x for x in filenames if (format in x) & ('xml' not in x) & (band in x)])
    # time = xr.Variable('time', time_index_from_filenames(files,41,33))
    time = xr.Variable('time', time_index_from_filenames(files,14,4))

    chunks = {'x': 5490, 'y': 5490, 'band': 1}
    da = xr.concat([xr.open_rasterio(os.path.join(folderName, f), chunks=chunks) for f in files], dim=time)

    # change DataArray values based on the time daysinmonth
    daysinmonth = da.time.dt.days_in_month
    da = da / daysinmonth / 24

    # set xarray da projection to epsg=32719
    da.rio.write_crs("EPSG:4326", inplace=True)
    da.to_netcdf(os.path.join('.', band+'.nc'))

    # open netcdf with xr
    path=r'/Users/farrospide/Library/CloudStorage/OneDrive-ciren.cl/Proyectos_RH_2024/ODEPA_EROSION/3_Desarrollo/0_Datos/pp/IMERG/prIMERG_2000_2023_mmmes.nc'
    debug = xr.open_dataset(path)

    # change data variable name precip by pp from debug DataArray
    debug = debug.rename_vars({'precip': 'pp'})

    # assign to debug as an attribute units: mm/mes
    debug.pp.attrs['units'] = 'mm/mes'

    # plot debug values from time index 2002-01-01
    debug.pp.sel(time='2002-01-01').plot()

    da.to_netcdf(path.replace('.nc','_mmhr.nc'))

    debug = xr.open_dataset(path)
    debug.pp.attrs['units'] = 'mm/hr'
    daysinmonth = debug.time.dt.days_in_month
    data = debug['pp'].values / daysinmonth.values[:, None, None, None]/24
    debug_new = xr.DataArray(data, coords=debug.coords, dims=debug.dims)
    debug_new.attrs = debug.attrs
    debug_new['pp'].attrs = debug['pp'].attrs

def pythonDate2eeDate(date=datetime.datetime):
    eeDate = date.strftime("%Y-%m-%d")
    return eeDate

def outFolder(folderName):
    try:
        os.mkdir(os.path.join('.', folderName))
    except:
        error = 'directorio ya existe'
    return None

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

def load_glacier_shapefile(path=os.path.join('.','geoData','glaciaresOlivares.shp')):
    glacier_shapefile = gpd.GeoDataFrame([],geometry=gpd.read_file(path).buffer(0))
    return glacier_shapefile

# calcular NIR/SWIR from landsat 8 surface reflectance
def calculate_nir_swir_landsat8(image):
    # get NIR from image
    nir = image.select('B5')
    # adjust scale and offset
    nir = nir
    # get SWIR from image
    # adjust scale and offset
    # mask swir where it is different from 0, else it is 1e-8
    swir = image.select('B6').where(image.select('B6').neq(0),1e-8)
    # return nir/swir renamed as nir_swir
    return nir.divide(swir).rename('nir_swir')

def download(collection, folderName):
    # login()

    # 2. define python dates
    fechaIni = datetime.datetime(2015, 1, 1)
    fechaFin = datetime.datetime(2016, 12, 1)
    # print(fechaIni, fechaFin)
    # print(pythonDate2eeDate(fechaFin))

    # 3. define the collection
    # collection = "ECMWF/ERA5_LAND/DAILY_AGGR"
    imageCollection = ee.ImageCollection(collection)

    # 4. define the bounding box and filter accordingly
    glacier_shapefile = load_glacier_shapefile().to_crs(epsg='4326')
    rectangle = gdf2FeatureCollection(glacier_shapefile)

    # llx, lly, urx, ury = -77,-58, -65, -17
    # rectangle = ee.Geometry.Rectangle([llx, lly, urx, ury])
    imageCollectionRegion = imageCollection.filterBounds(rectangle)

    # 5. select the band
    # band = "evaporation_from_open_water_surfaces_excluding_oceans_sum"
    # map calculate_nir_swir_landsat8 over landsat8
    nir_swir = imageCollectionRegion.map(calculate_nir_swir_landsat8)

    imageCollectionBand = nir_swir.select('nir_swir')
    f = open("log.txt", "w")

    # 6. output folder
    outFolder(folderName)
    rutaDl = os.path.join(folderName)

    # 7. select date
    # iterate over imageCollectionBand dates
    
    for feature in imageCollectionBand.getInfo()['features']:
        date = feature['properties']['system:time_start']
        print(f"Attempting download image for {date}")
        try:
            imageCollectionDate = imageCollectionBand.filterDate(date).min().clip(rectangle)

            url = imageCollectionDate.getDownloadURL({
                'scale': 30,
                'crs': 'EPSG:32719',
                'format': 'GEO_TIFF',
                'region': rectangle})
            print(url)
            # print(f"Attempting download image for {date} OK!")
            response = requests.get(url, stream=True)
            filePath = os.path.join(rutaDl, f"tmax{date}.tif")

            # zip
            with open(os.path.join(rutaDl, band+'_'+date+'.tif'), "wb") as fd:
                for chunk in response.iter_content(chunk_size=1024):
                    fd.write(chunk)
            fd.close()

        except:
            # print(f"Attempting download image for {date} FAILED!")
            f.write(f"{date}\n")


if __name__ == '__main__':
    collection = "LANDSAT/LC08/C02/T1_TOA"
    band = "B5"  # total_precipitation_sum, temperature_2m, 
    folderName = "l8"
    outFolder(folderName=folderName)
    
    # 1. authenticate
    ee.Authenticate()
    ee.Initialize()
    download(collection=collection, folderName=folderName)
    files2NETCDF(folderName=folderName)

# open .nc file
debug=xr.open_dataset(os.path.join("/", "Users", "farrospide", "Downloads",'ERA5', 'snowDepth_netcdf.nc'))
# plot
debug.band[0].plot()

import rioxarray as rxr
debug=rxr.open_rasterio(os.path.join(folderName,'ERA5', 'netcdf.nc'))
debug.plot()

import geemap
ndwiViz = {'min': 0.5, 'max': 10, 'palette': ['00FFFF', '0000FF']}
Map=geemap.Map()
Map.addLayer(imageCollectionDate, ndwiViz, 'pp TerraClimate 1990-01-01')
Map