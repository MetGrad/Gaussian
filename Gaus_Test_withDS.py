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
import glob
from scipy.ndimage import gaussian_filter
import numpy as np

#Open HYCOM 25
ds_25 = xr.open_mfdataset('/Users/Anna/Desktop/MSMET/Thesis/WRFprepData/March2018HYCOM/020_archv.2018_084_*_3z.nc', combine = 'by_coords', concat_dim = 'time')
ds_25.to_netcdf('test_25.nc')
ds25 = xr.open_dataset('test_25.nc')

#test domain
latbounds = [26, 27]
lonbounds = [-89.5, -88.5]

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

SST_subset_25=ds25.variables['water_temp'][:,0,latli:latui, lonli:lonui]
#SST_subset_26=ds26.variables['water_temp'][:,0,latli:latui, lonli:lonui]

u_subset_25=ds25.variables['u'][:,0,latli:latui, lonli:lonui]
#v_subset_25=ds25.variables['v'][:,0,latli:latui, lonli:lonui]
#u_subset_26=ds26.variables['u'][:,0,latli:latui, lonli:lonui]
#v_subset_26=ds26.variables['v'][:,0,latli:latui, lonli:lonui]

lat_array=ds25.Latitude.values[latli:latui]
lon_array=ds25.Longitude.values[lonli:lonui]

meanU_25 = u_subset_25.mean('MT')
#meanV_25 = v_subset_25.mean('MT')
#meanU_26 = u_subset_26.mean('MT')
#meanV_26 = v_subset_26.mean('MT')

#meanSST_25 = SST_subset_25.mean('MT')
#meanSST_26 = SST_subset_26.mean('MT')

#March 25 Dopplerscat
ds25_6 = xr.open_dataset('/Users/Anna/Desktop/MSMET/Thesis/JPL_March22/20180325_084330_0380-0436_line06.L2.nc')

speed25_6 = ds25_6.wind_speed.values
dir25_6 = ds25_6.wind_dir.values
dir_math25_6=90-dir25_6*(np.pi/180)
uWind25_6=speed25_6*np.cos(dir_math25_6*(180/np.pi))
lat25_6=ds25_6.latitude.values
lon25_6=ds25_6.longitude.values
flag25_6=ds25_6.flag.values
u_cur25_6=ds25_6.u_current.values
v_cur25_6=ds25_6.v_current.values

          
lonlen25_6=len(lon25_6)
latlen25_6=len(lat25_6)
goodUarr25_6=np.zeros((latlen25_6,lonlen25_6))
goodUarr25_6[:]=np.nan
#goodVarr25_6=np.zeros((latlen25_6,lonlen25_6))
#goodVarr25_6[:]=np.nan
#goodWindUarr25_6=np.zeros((latlen25_6,lonlen25_6))
#goodWindUarr25_6[:]=np.nan
#goodWindVarr25_6=np.zeros((latlen25_6,lonlen25_6))
#goodWindVarr25_6[:]=np.nan
for x in range (len(lon25_6)):
    for y in range (len(lat25_6)):
        if flag25_6[y,x] == 0:
            goodUarr25_6[y,x] = u_cur25_6[y,x] 

outerDS25_6cu = []
count_latDS25_6cu = 0 
for x in lat25_6:
    count_lonDS25_6cu = 0
    for y in lon25_6: 
        innerDS25_6cu= []
        if not (np.isnan(goodUarr25_6[count_latDS25_6cu][count_lonDS25_6cu].data)): 
            innerDS25_6cu.append(x)
            innerDS25_6cu.append(y)
            #inner_UK.append(meanU_25[count_lat][count_lon].data)
            innerDS25_6cu.append(goodUarr25_6[count_latDS25_6cu][count_lonDS25_6cu])
            outerDS25_6cu.append(innerDS25_6cu)
        count_lonDS25_6cu = count_lonDS25_6cu + 1
    count_latDS25_6cu = count_latDS25_6cu + 1  
DS25_6_uCurrent = np.array(outerDS25_6cu)

lat25_6 = DS25_6_uCurrent[:,0]
lon25_6 = DS25_6_uCurrent[:,1]
#uWind25_6 = DS25_6_uWind[:,2]
#vWind25_6 = DS25_6_vWind[:,2]
uCurrent25_6 = DS25_6_uCurrent[:,2]

#Gaussian Filter
#HYCOM: lat_array (Y), lon_array (X), meanU_25
#DS: goodUarr25_1-10, lat25_1, lon25_1s

HYCOM_lat = lat_array
HYCOM_lon = lon_array

#Define new grid: 1km (1000m)
#match WRF refernce lat/lon, increment by partial degree 
###### !!!!! wrf lat lon NEED TO UPDATE THESE !!!!! ######

#REAL CASE: 
#lat_0 =  24.5
#lon_0 = -92
#lat_max = 31
#lon_max = -87.5

#TEST CASE
lat_0 =  26
lon_0 = -89.5
lat_max = 27
lon_max = -88.5 

#min, max, then increment and fill lat & then lon (Y=lat, X=lon)
#0.008 degrees = 1km 
WRF_lat = np.arange(lat_0, lat_max, 0.00898315)
WRF_lon = np.arange(lon_0, lon_max, 0.00898315)
WRF_latLen = len(WRF_lat)
WRF_lonLen = len(WRF_lon)

#make current arrays (u & v) that match lat & lon dimensions, fill with nans 
u25_WRF  = np.zeros((WRF_latLen, WRF_lonLen))   #(813, 688)
#v25_WRF = np.zeros((WRF_latLen, WRF_lonLen))  
u25_WRF[:] = -999
#v25_WRF[:] = -999
#u26_WRF  = np.zeros((WRF_latLen, WRF_lonLen))   #(813, 688)
#v26_WRF = np.zeros((WRF_latLen, WRF_lonLen))  
#u26_WRF[:] = -999
#v26_WRF[:] = -999
#sst25_WRF  = np.zeros((WRF_latLen, WRF_lonLen))   #(813, 688)
#sst26_WRF = np.zeros((WRF_latLen, WRF_lonLen))  
#sst25_WRF[:] = -999
#sst26_WRF[:] = -999

#at equator, 1deg = 111.31949077920639km, so 1km*(1deg/111.31949077920639km) = 0.00898315
#Great Circle Calculator 
#http://edwilliams.org/gccalc.htm

#DS sigma??? = 0.03 # 0.015*2   
km_per_deg= 1/0.00898315
sigma_HYCOM = 0.02 # 0.015*2 #For HYCOM 
sigma_HYCOM_km = sigma_HYCOM*km_per_deg
sigma_DS = 0.5
sigma_DS_km = sigma_DS*km_per_deg


for y in range (0, WRF_latLen):  
    #print(' y :', y)
    for x in range (0, WRF_lonLen):
        #print('x:', x)
        sum_weight_HYCOM = 0
        sum_product_u25HYCOM = 0
        #sum_product_v25HYCOM = 0
        #sum_product_u26HYCOM = 0
        #sum_product_v26HYCOM = 0
        #sum_product_sst25HYCOM = 0
        #sum_product_sst26HYCOM = 0 
       
        #define area of weighting (lat,lon point +- 3 sigma in lat/lon but round up)
        #+- half size 
        #find all HYCOM indeces that are within WRF smallbox 
        #set bounds in lat/lon
        min_lon = max(WRF_lon[x] - 3*sigma_HYCOM, WRF_lon[0])
        max_lon = min(WRF_lon[x] + 3*sigma_HYCOM, WRF_lon[WRF_lonLen -1])
        min_lat = max(WRF_lat[y] - 3*sigma_HYCOM, WRF_lat[0])   
        max_lat = min(WRF_lat[y] + 3*sigma_HYCOM, WRF_lat[WRF_latLen -1])       

        good_lon = [] 
        good_lat = []
        inter_like = []
        #print('min lat: ', min_lat, 'WRFlat_y: ',WRF_lat[y], 'max lat: ', max_lat)
        #print('min lon: ', min_lon, 'WRFlon_x: ',WRF_lon[x], 'max lon: ', max_lon)
        #print('HYCOMlon: ',HYCOM_lon)
        
        #HYCOM lat/lon comparison 
        test1_In = np.where(HYCOM_lon >= min_lon)
        test2_In = np.where(HYCOM_lon <= max_lon)
        inter_like_lon = np.intersect1d(test1_In, test2_In)
        
        test1_In = np.where(HYCOM_lat >= min_lat)
        test2_In = np.where(HYCOM_lat <= max_lat)
        inter_like_lat = np.intersect1d(test1_In, test2_In)   
        
        #DS lat/lon comparison
        DS25_test6_1 = np.where(lon25_6 >= min_lon)
        DS25_test6_2 = np.where(lon25_6 <= max_lon)
        DS25_inter_lon6 = np.intersect1d(DS25_test6_1, DS25_test6_2)
        
        DS25_test6_1 = np.where(lat25_6 >= min_lat)
        DS25_test6_2 = np.where(lat25_6 <= max_lat)
        DS25_inter_lat6 = np.intersect1d(DS25_test6_1, DS25_test6_2)
        #np.where((HYCOM_lon >= min_lon) and (HYCOM_lon <= max_lon), HYCOM_lon)
        #np.where((HYCOM_lat >= min_lat) and (HYCOM_lat <= max_lat))
        
        #HYCOM index
        lon_indices = inter_like_lon
        lat_indices = inter_like_lat
        #DS index
        DS25_lon_In6 = DS25_inter_lon6
        DS25_lat_In6 = DS25_inter_lat6
        
        #lon_indices = np.logical_and(HYCOM_lon >= min_lon, HYCOM_lon <= max_lon)  
        #lat_indices = np.logical_and(HYCOM_lat >= min_lat, HYCOM_lat <= max_lat)
        #print('shapes:', min_lon.shape, max_lon.shape)
        #print('lon_indicies: ',lon_indices,'\nlat_indicies: ', lat_indices)
        
        #Good HYCOM index
        good_lon = HYCOM_lon[lon_indices]
        good_lat = HYCOM_lat[lat_indices]       
        #print('good lon: ',good_lon,  'good lat:', good_lat)
        #print('good lat shape: ', good_lat.shape, 'good lon shape: ', good_lon.shape)
        
        #Good DS index
        DS25_good_lon6 = lon25_6[DS25_lon_In6]
        DS25_good_lat6 = lat25_6[DS25_lat_In6]
    
        #get lengths of good arrays HYCOM 
        good_lon_len = len(good_lon)
        good_lat_len = len(good_lat)
        #print('good lat len: ',good_lat_len, 'good lon len: ', good_lon_len )
        #DS good lat/lon length 
        DSgood_lon_len6 = len(DS25_good_lon6)
        DSgood_lat_len6 = len(DS25_good_lat6)        
        
        
        #loop through values of good arrays 
        for y_small in range(0, good_lat_len):
            for x_small in range(0, good_lon_len):         
                
                #calcualte the distance, 1km  
                #print('x, y :', x, y)
                #print('x_small, y_small:', x_small, y_small)
                distance=(km_per_deg)*np.sqrt((good_lat[y_small]-WRF_lat[y])**2+((good_lon[x_small]-WRF_lon[x])*np.cos(WRF_lat[y]*(np.pi/180.0)))**2)
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
                sum_product_u25HYCOM = sum_product_u25HYCOM + weight_HYCOM * meanU_25[lat_indices[y_small], lon_indices[x_small]]
                #sum_product_v25HYCOM = sum_product_v25HYCOM + weight_HYCOM * meanV_25[lat_indices[y_small], lon_indices[x_small]]
                #sum_product_u26HYCOM = sum_product_u26HYCOM + weight_HYCOM * meanU_26[lat_indices[y_small], lon_indices[x_small]]
                #sum_product_v26HYCOM = sum_product_v26HYCOM + weight_HYCOM * meanV_26[lat_indices[y_small], lon_indices[x_small]]
                #sum_product_sst25HYCOM = sum_product_sst25HYCOM + weight_HYCOM * meanSST_25[lat_indices[y_small], lon_indices[x_small]]
                #sum_product_sst26HYCOM = sum_product_sst26HYCOM + weight_HYCOM * meanSST_26[lat_indices[y_small], lon_indices[x_small]]

    #outside of loop: sum of products/sum of weights
        u25_WRF[y,x] = sum_product_u25HYCOM/sum_weight_HYCOM 
        #v25_WRF[y,x] = sum_product_v25HYCOM/sum_weight_HYCOM 
        #u26_WRF[y,x] = sum_product_u26HYCOM/sum_weight_HYCOM 
        #v26_WRF[y,x] = sum_product_v26HYCOM/sum_weight_HYCOM
        #sst25_WRF[y,x] = sum_product_sst25HYCOM/sum_weight_HYCOM
        #sst26_WRF[y,x] = sum_product_sst26HYCOM/sum_weight_HYCOM
#end loop
