import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import num2date, date2num
from plotting import PlotData
remote = False

# file directories:
local_directory = 'C:/Users/ciank/PycharmProjects/sinmod/sg800_krill/'
if remote:
    remote_folder = 'A:/Cian_sinmod/sg800_krill/'
    trajectory_folder = remote_folder + 'trajectory/'
else:
    trajectory_folder = local_directory + 'trajectory/'
save_folder = local_directory + 'figures/'

# initialise plotting class
pl = PlotData(save_folder, trajectory_folder)


# read the file that needs reading: #todo: functions that call different files, determines number of plots and saves if time
time_prefix = 'trajectory_20170301_d40'
pl.read_trajectory_file(file_prefix=time_prefix)
#breakpoint()

pl.plot_growth(skip_t=1, kk=0)

#pl.plot_anomalies(skip_t=1, kk=0)

#pl.plot_depth_region(skip_t=1, kk=6)
# plot depths:
#pl.plot_depths(skip_t=1, kk=6)

# plot temperature, speeds etc.
#pl.plot_temp()

# first plots: environment and trajectories over time;
#pl.plot_currents()
#pl.plot_trajectory_color(skip_n=10, skip_t=3, kk=0)
#pl.plot_retention()

# dom_pathways plots:
#pl.plot_dom_pathways(skip_t=1, kk=slice(0,1,1))



breakpoint()


