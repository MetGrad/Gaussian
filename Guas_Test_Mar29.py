import math
import xarray as xr
import cartopy.crs as ccrs
import matplotlib as mpl
from matplotlib import pyplot as plt 
import matplotlib.pylab as plt
import matplotlib.path as mpath
from matplotlib.pyplot import figure
import numpy as np
import cartopy.feature as cfeature
import pyproj
import utm
import pandas as pd 
import traceback 
from pathlib import Path  
import os  
import numpy.ma as ma
import netCDF4 as nc

from scipy.ndimage import gaussian_filter
import numpy as np

#Open HYCOM 25
ds_25 = xr.open_mfdataset('/Users/Anna/Desktop/MSMET/Thesis/WRFprepData/March2018HYCOM/020_archv.2018_084_*_3z.nc', combine = 'by_coords', concat_dim = 'time')
ds_25.to_netcdf('test_25.nc')
ds25 = xr.open_dataset('test_25.nc')

#test domain
latbounds = [27, 28]
lonbounds = [-89, -88]

#real domain 
#latbounds = [24.5, 31]
#lonbounds = [-92, -87.5]

lats = ds25.Latitude.values[:]
lons = ds25.Longitude.values[:]

#set upper & lower bounds
latli = np.argmin(np.abs(lats-latbounds[0]))
latui = np.argmin(np.abs(lats-latbounds[1]))
lonli = np.argmin(np.abs(lons-lonbounds[0]))
lonui = np.argmin(np.abs(lons-lonbounds[1]))

u_subset_25=ds25.variables['u'][:,0,latli:latui, lonli:lonui]
v_subset_25=ds25.variables['v'][:,0,latli:latui, lonli:lonui]

lat_array=ds25.Latitude.values[latli:latui]
lon_array=ds25.Longitude.values[lonli:lonui]

meanU_25 = u_subset_25.mean('MT')
meanV_25 = v_subset_25.mean('MT')

#Plot U25 current 
mapcrs = ccrs.PlateCarree()
fig = plt.figure(figsize=(8, 4), dpi=150)
ax = fig.add_subplot(111, projection=mapcrs)
ax.coastlines()
cs = ax.contourf(lon_array, lat_array, meanU_25)
#ax.contour(lon_array, lat_array, meanU_25, colors = 'white', linewidths = 1, transform=ccrs.PlateCarree())
PCM=ax.get_children()[2]
plt.colorbar(cs, ax=ax, label = 'Speed m/s')
gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, color = 'black')
plt.tight_layout()
gl.top_labels = False
gl.ylabels_right = False
plt.title('Test Domain:\nHYCOM Mean u-Current for March 25, 2018')
plt.show()

#Gaussian Filter
#HYCOM: lat_array (Y), lon_array (X), meanU_25
#DS: goodUarr25_1-10, lat25_1, lon25_1s

#print('meanU25 shape: ', meanU_25.shape)
HYCOM_lat = lat_array
HYCOM_lon = lon_array
#print('HYCOM lat shape: ', HYCOM_lat.shape)
#print('HYCOM lon shape:', HYCOM_lon.shape)

#Define new grid: 1km (1000m)
#match WRF refernce lat/lon, increment by partial degree 
###### !!!!! wrf lat lon NEED TO UPDATE THESE !!!!! ######

#REAL CASE: 
#lat_0 =  24.5
#lon_0 = -87.5
#lat_max = 31
#lon_max = -92

#TEST CASE
lat_0 =  27
lon_0 = -89
lat_max = 28
lon_max = -88 

#min, max, then increment and fill lat & then lon (Y=lat, X=lon)
#0.008 degrees = 1km 
WRF_lat = np.arange(lat_0, lat_max, 0.008)
WRF_lon = np.arange(lon_0, lon_max, 0.008)
WRF_latLen = len(WRF_lat)
WRF_lonLen = len(WRF_lon)
#print('WRF lon shape:', WRF_lon.shape)
#print('WRF lat shape:', WRF_lat.shape)
#print('WRF lon: ', WRF_lon)
#make current arrays (u & v) that match lat & lon dimensions, fill with nans 
u25_WRF  = np.zeros((WRF_latLen, WRF_lonLen))   #(813, 688)
v25_WRF = np.zeros((WRF_latLen, WRF_lonLen))  
u25_WRF[:] = -999
v25_WRF[:] = -999
#at equator, 1km = 0.008 degrees, round up to 0.01
#DS sigma??? = 0.03 # 0.015*2   
deg_to_km = 1/0.008
sigma_HYCOM = 0.02 # 0.015*2 #For HYCOM 
sigma_HYCOM_km = sigma_HYCOM*deg_to_km


for y in range (0, WRF_latLen):  
    print(' y :', y)
    for x in range (0, WRF_lonLen):
        #print('x:', x)
        sum_weight_HYCOM = 0
        sum_product_u25HYCOM = 0
       
        #define area of weighting (lat,lon point +- 3 sigma in lat/lon but round up)
        #+- half size 
        #small_box = sigma_HYCOM/(HYCOM_lat[1] - HYCOM_lat[0]) + 1 #if rounds down, 1 gives enough space 
        #print(HYCOM_lat[1] - HYCOM_lat[0])
        #print('small box: ', small_box)
        
        #find all HYCOM indeces that are within WRF smallbox 
        #set bounds in lat/lon
        min_lon = max(WRF_lon[x] - 3*sigma_HYCOM, WRF_lon[0])
        max_lon = min(WRF_lon[x] + 3*sigma_HYCOM, WRF_lon[WRF_lonLen -1])
        min_lat = max(WRF_lat[y] - 3*sigma_HYCOM, WRF_lat[0])   
        max_lat = min(WRF_lat[y] + 3*sigma_HYCOM, WRF_lat[WRF_latLen -1])       

        good_lon = [] 
        good_lat = []
        #print('min lat: ', min_lat, 'WRFlat_y: ',WRF_lat[y], 'max lat: ', max_lat)
        #print('min lon: ', min_lon, 'WRFlon_x: ',WRF_lon[x], 'max lon: ', max_lon)
        #print('HYCOMlon: ',HYCOM_lon)
        
        #for q in HYCOM_lon: 
        #    if min_lon <= q <= max_lon:
        #        good_lon.append(q)
        #for r in HYCOM_lat:
        #    if min_lat <= r <= max_lat:
        #        good_lat.append(r)
        
        lon_indices = np.logical_and(HYCOM_lon >= min_lon, HYCOM_lon <= max_lon)  
        lat_indices = np.logical_and(HYCOM_lat >= min_lat, HYCOM_lat <= max_lat)
        #print('shapes:', min_lon.shape, max_lon.shape)
        print('lon_indicies: ',lon_indices,'lat_indicies: ', lat_indices)
        #print('min lon: ', min_lon)
        #print('max lon: ', max_lon)
        #print('min lat: ', min_lat)
        #print('max lat: ', max_lat)
        
        good_lon = HYCOM_lon[lon_indices]
        good_lat = HYCOM_lat[lat_indices]       
        #print('good lon: ',good_lon,  'good lat:', good_lat)
        #print('good lat:', good_lat)
        #print('good lon:', good_lon)
        
        #print('good lat shape: ', good_lat.shape, 'good lon shape: ', good_lon.shape)
    
        #get lengths of good arrays 
        #good_lon = np.array(good_lon)
        #good_lat = np.array(good_lat)
        good_lon_len = len(good_lon)
        good_lat_len = len(good_lat)
        #print('good lat len: ',good_lat_len, 'good lon len: ', good_lon_len )
        
        #loop through values of good arrays 
        #print(' y :',y)
        for y_small in range(0, good_lat_len):
            for x_small in range(0, good_lon_len):         
                
                #calcualte the distance, 1km  
                #print('x, y :', x, y)
                #print('x_small, y_small:', x_small, y_small)
                #if x_small >0 and y_small >0:
                #distance=(deg_to_km)*np.sqrt((HYCOM_lat[y_small]-WRF_lat[y])**2+((HYCOM_lon[x_small]-WRF_lon[x])*np.cos(WRF_lat[y]*(np.pi/180.0)))**2)
                distance=(deg_to_km)*np.sqrt((good_lat[y_small]-WRF_lat[y])**2+((good_lon[x_small]-WRF_lon[x])*np.cos(WRF_lat[y]*(np.pi/180.0)))**2)
                #print('distance: ', distance)
                #print('3sigma', 3*sigma_HYCOM)
                #print(y, x, y_small, x_small, distance, good_lat[y_small], WRF_lat[y], good_lon[x_small], WRF_lon[x])
                #determine weights (u & v) - get an array of poitns of the smaller domain
                if distance <= 3 * sigma_HYCOM_km: 
                    weight_HYCOM = 0.2*(1 / (sigma_HYCOM_km * np.sqrt(2 * np.pi))) * np.exp(-(distance**2) / (2*sigma_HYCOM_km**2))
                else: 
                    weight_HYCOM = 0 
                    
                #sum wieghts for HYCOM u & v
                    #print('weight_HYCOM: ', weight_HYCOM)
                sum_weight_HYCOM = sum_weight_HYCOM + weight_HYCOM
                #print('sum_weight_HYCOM: ', sum_weight_HYCOM)

        #if sum of weights != 0, sum product of weights* current value for HYCOM 
        #0.2*HYCOM weights 
                                                                                        #we want HYCOM indices 
                sum_product_u25HYCOM = sum_product_u25HYCOM + weight_HYCOM * meanU_25[HYCOM_lat[y_small], HYCOM_lon[x_small]]
                    #print('weight hycom, meanu25 @ index: ', weight_HYCOM, meanU_25[y_small, x_small])
               
        #end small box loop 
    #outside of loop: sum of products/sum of weights (u & v)
    #if sum_weight_HYCOM != 0: 
        #print('sum product, sum weight: ', sum_product_u25HYCOM, sum_weight_HYCOM)
        u25_WRF[y,x] = sum_product_u25HYCOM/sum_weight_HYCOM 
    #else: 
    #    u25_WRF[y,x] = 0 
#end loop

print(u25_WRF.shape)
print(u25_WRF)

#plot smoothed 25 HYCOM (u and v) 

mapcrs = ccrs.PlateCarree()
fig = plt.figure(figsize=(8, 4), dpi=150)
ax = fig.add_subplot(111, projection=mapcrs)
ax.coastlines()
cs = ax.contourf(WRF_lon, WRF_lat, u25_WRF)
#ax.contour(WRF_lon, WRF_lat, u25_WRF, colors = 'white', linewidths = 1, transform=ccrs.PlateCarree())
PCM=ax.get_children()[2]
plt.colorbar(cs, ax=ax, label = 'Speed m/s')
gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, color = 'black')
plt.tight_layout()
gl.top_labels = False
gl.ylabels_right = False
plt.title('HYCOM u-Current for March 25, 2018\non WRF Grid')
plt.show()