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
local_folder = 'C:/Users/ciank/PycharmProjects/sinmod/sg800_krill/'
trajectory_folder = local_folder + 'trajectory/'

# time parameters: #todo: have a dateinit variable that calls the correct file
init_year = 2017
init_month = 5

# model parameters:
n = 10000
x_min = 210.0
x_max = 610.4
y_min = 210.0
y_max = 610.5

# time parameters
duration_days = datetime.timedelta(days=30)  # simulation duration in days;
minutes = 4
dt = datetime.timedelta(hours=minutes/60)
save_step = datetime.timedelta(hours=4)  # save time step
time_threshold = datetime.timedelta(seconds=0)
simulation_steps = duration_days/dt
save_number = duration_days/save_step

# save_file_name
save_file_prefix = 'trajectory_d' + str(duration_days.days)

# reader that adds info from current file
reader_SG = SGReader(samples_folder, samples_prefix, init_month, init_year)

# start of simulation:
time_counter = save_step
save_counter = -1
for i in np.arange(0,simulation_steps,1):
    time_counter -= dt
    if i == 0:
        k = Krill(reader_SG.nc_file)
        k.init_krill(n=n, x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max)
        k.init_netcdf(trajectory_folder=trajectory_folder, n=n, save_n=save_number, save_file_prefix=save_file_prefix)

    k.step_krill(dt, reader_SG)
    reader_SG.update_time(dt)
    # save simulation
    if time_counter <= time_threshold:
        save_counter += 1
        time_counter = save_step
        print('saving ' + str(save_counter + 1) + ' of ' + str(save_number))
        k.save_step(save_counter=save_counter, current_datetime=reader_SG.current_datetime)

k.trajectory_file.close()
