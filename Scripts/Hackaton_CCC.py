#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 23:00:30 2020

@author: felipe
"""

# importar librerias


#%%
def main(folder = './', nc_cloud = 'ls8_toa_2013_2020_icesnow.nc' , nc_glaciar = 'ls8_toa_2013_2020_glaciar.nc', df_times = './time_2013_2020_toa.csv'):
    #%%
    
    import geopandas
    from matplotlib import pyplot as plt
    import rasterio
    import numpy as np
    import pandas as pd
    from cdo import Cdo
    import gc
    import rioxarray
    import xarray
    from shapely.geometry import mapping

    cob_nieve = folder + '//'+ nc_cloud
    cob_glacial = folder + '//'+ nc_glaciar
    nieve = xarray.open_dataset(cob_nieve)['pixel_qa']
    nirswir = xarray.open_dataset(cob_glacial)['nir_swir']
    times = pd.read_csv(folder + '//'+ df_times,index_col = 0, parse_dates = True).index
    # assign times
    
    nirswir = nirswir.assign_coords(time = times)
    nieve = nieve.assign_coords(time = times)  
    clip = xarray.DataArray(np.zeros((170,1679,1416)), coords=nirswir.coords)
    clip = clip.assign_coords(time = times)  


    #coordenadas
    nirswir.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=True)
    nirswir.rio.write_crs("epsg:32719", inplace=True)

    nieve.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=True)
    nieve.rio.write_crs("epsg:32719", inplace=True)    
    ############# Read shapefile and first feature
    cobertura_nival = geopandas.read_file(folder + '//'+'glaciares.shp', crs="epsg:32719")
    
    areas_glaciares = np.zeros((len(times)))
       
    ############# Read NetCDF variables
    default = 2
    dicc = {1 : default, 2: default, 3: default ,4: default ,5: default , 6: 2.2, 7: 2.1, 8: 3, 
            9: 1.9, 10: default, 11: 1.9, 12: 1.5}
    #for i in range(round(len(tiempo[0:41]))): 
    for i in range(round(len(times))): 
    
    
        nc_ua  = nirswir[i,:,:]
        snow = nieve[i,:,:]
                   
        thres = dicc[times[i].month]
        thres = 3
        
        #Filtrado de nubes y relleno
        
        nc_ua = nc_ua.where(~snow.isin([352, 368, 416, 432, 480, 864, 880, 928, 944, 992]), other=np.nan) 
        try:
            nc_ua = nc_ua.interpolate_na(dim="x", method="nearest", fill_value="extrapolate")
        except:
            areas_glaciares[i] = 0
            print('muchas nubes')
            clip[i] = np.nan
            gc.collect()
            del nc_ua
            del snow
            continue
        nc_ua = nc_ua.where(snow.isin([336, 368, 400, 432, 848, 880, 912, 944, 1352, 324, 388, 836, 900, 1348]), other=np.nan) 
    
        # filtrado de glaciares
        nc_ua = nc_ua.where(((nc_ua >= thres)), other=np.nan) 

        
    # Calculate mask cuenca
    
        clipped = nc_ua.rio.clip(cobertura_nival.geometry.apply(mapping), cobertura_nival.crs, drop=False)
    
    #numero pixeles glaciar dentro de la cuenca

        clip[i] =   clipped.values   
        area_glaciar = np.count_nonzero(~np.isnan(clipped))*900/1e6
        areas_glaciares[i] = area_glaciar
        gc.collect()
        del nc_ua
        del clipped

    nirswir.close()
    nieve.close()
  #%% Graficar
          

    plt.close("all")
    plot_areas = pd.DataFrame(areas_glaciares, index = times)
    plot_areas.columns = ['Area']
    plot_areas = plot_areas[plot_areas['Area'] > 20]
    retroceso = (plot_areas.resample('Y').mean().iloc[-1]-plot_areas.resample('Y').mean().iloc[0])/plot_areas.resample('Y').mean().iloc[0]
    ax = plot_areas.rolling(15).mean().plot(rot = 90)
    ax.set_ylabel('Cobertura glaciar ($km^2$)')
    ax.set_xlabel('')
    ax.set_ylim(bottom =0)
    ax.grid()
    
    media_first = clip.sel(time=slice('2013-01-01', '2013-12-31'))
    media_first_mean = media_first.resample(time = '1Y').mean()[0]
    media_first_mean.to_netcdf('cobertura_glaciar_media_2013.nc')
    media_last =  clip.sel(time=slice('2020-01-01', '2020-08-21'))
    media_last_mean = media_last.resample(time = '1Y').mean()[-1]
    media_last_mean.to_netcdf('cobertura_glaciar_media_2020.nc')

    
    x = media_last_mean.x.values
    y = media_last_mean.y.values
    
#   
    plt.figure(figsize = (16,19)) 
    plt.pcolormesh(x, y, media_first_mean.values, cmap = 'afmhot', vmin = 3, vmax = 200, linewidths = 4, label='y0')

    plt.pcolormesh(x, y, media_last_mean.values, cmap = 'cool', vmin = 3, vmax = 200, linewidths = 0.1, label='y1')
        
    plt.xlabel('Este (m)')
    plt.ylabel('Norte (m)')

def Temperatura(df_folder = '/home/carlos/Downloads/modis_temp.xlsx'):
    df_t = pd.read_excel(df_folder, index_col = 1, parse_dates = True)
    df_t['LST'].resample('YS').mean().plot()
    plt.ylabel('Temperatura media de la superficie terrestre (°C)')
    
#%%
if __name__ == '__main__':
    main()