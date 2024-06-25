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
import pygrib

#Open HYCOM 25
ds_25 = xr.open_mfdataset('/Users/Anna/Desktop/MSMET/Thesis/WRFprepData/March2018HYCOM/020_archv.2018_084_*_3z.nc', combine = 'by_coords', concat_dim = 'time')
ds_25.to_netcdf('test_25.nc')
ds25 = xr.open_dataset('test_25.nc')

#ds_26 = xr.open_mfdataset('/Users/Anna/Desktop/MSMET/Thesis/WRFprepData/March2018HYCOM/020_archv.2018_085_*_3z.nc', combine = 'by_coords', concat_dim = 'time')
#ds_26.to_netcdf('test_26.nc')
#ds26 = xr.open_dataset('test_26.nc')

#test domain
latbounds = [26.5, 27]
lonbounds = [-91, -90.5]

#real domain 
#latbounds = [24.5, 29]
#lonbounds = [-93, -87.5]

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