from reader_SG800 import geo2grid
import netCDF4 as nc
import math
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import xarray as xr
import pyproj
import numpy as np
samples_folder = 'E:/fromIngrid/samplefiles_reruns/'
samples_prefix = 'samplesNSEW_'
samples_year = 2017
samples_month = 2

filename = samples_folder + samples_prefix + str(samples_year) + "{:02d}".format(samples_month) + '.nc'
nc_file = xr.open_dataset(filename, decode_times=False)

grid_mapping = nc_file['grid_mapping'].attrs


fig, axs = plt.subplots(figsize=(12, 8), layout='constrained')
axs.pcolormesh(nc_file['temperature'][0,0,:,:], cmap=plt.get_cmap('hot'))
breakpoint()

#geo2grid(-42, -57, case='get_xy')


proj_daymet = "+proj=lcc +lat_0=42.5 +lon_0=-100 +lat_1=25 +lat_2=60 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" #my custom CRS




crs = pyproj.CRS.from_cf(grid_mapping)
proj4 = crs.to_proj4()
p1 = pyproj.Proj(proj4)
x = 100
y=100
x_0=-175200
y_0 = 1484800
gm = nc_file['grid_mapping']
x0 = gm.false_easting
y0 =  gm.false_northing
dx = dy_m = gm.horizontal_resolution
lat_0 = gm.latitude_of_projection_origin
lon_0 = gm.longitude_of_projection_origin
lat_ts= gm.standard_parallel

proj_str = '+proj=lcc'
#proj4string = f'lcc +lat_0={lat_0} +lon_0={lon_0} +lat_1={lat_1} +lat_2={lat_1} +a={radius} +b={radius} +x_0={x0} +y_0={y0}'
proj4string = f'stere +lat_0={lat_0} +lon_0={lon_0} +lat_ts={lat_ts} +a={radius} +b={radius} +x_0={x0} +y_0={y0}'
proj = pyproj.Proj(proj=proj4string)
lon, lat = proj(x, y, inverse=True)
p1(x_0,y_0,inverse=True)
breakpoint()


