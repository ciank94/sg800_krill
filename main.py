from reader_SG800 import SGReader, geo2grid
from SG_krill import Krill
import netCDF4 as nc
import datetime
import math
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import xarray as xr
import pyproj
import numpy as np
samples_folder = 'E:/fromIngrid/samplefiles_reruns/'
samples_prefix = 'samplesNSEW_'


# model parameters:
n = 10000
x_min = 110.0
x_max = 810.4
y_min = 210.0
y_max = 610.5

# time parameters:
init_year = 2017
init_month = 2
dt = datetime.timedelta(hours=12)


# reader that adds info from current file
reader_SG = SGReader(samples_folder, samples_prefix, init_month, init_year)

for i in range(0, 10):

    if i == 0:
        k = Krill(reader_SG.nc_file)
        k.init_krill(n=n, x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max)

    k.step_krill(dt, reader_SG)
    reader_SG.update_time(dt)
    k.plot_currents()


breakpoint()









#todo: add functionality to read to read multiple files










#todo: initialise case of krill if first instance to get grid information



breakpoint()








k.plot_init()





grid_mapping = nc_file['grid_mapping'].attrs



breakpoint()

#geo2grid(-42, -57, case='get_xy')





breakpoint()


