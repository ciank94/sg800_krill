import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import num2date, date2num
from plotting import PlotData

local_directory = 'C:/Users/ciank/PycharmProjects/sinmod/sg800_krill/'
save_folder = local_directory + 'figures/'
trajectory_folder = local_directory + 'trajectory/'
pl = PlotData(save_folder, trajectory_folder)

# read the file that needs reading: #todo: functions that call different files, determines number of plots and saves if time
time_prefix= 'sim_201705_d30'
pl.read_trajectory_file(file_prefix=time_prefix)

# first plots: environment and trajectories over time;
#pl.plot_currents()
#pl.plot_trajectory_color(skip_n=5, skip_t=1, kk=0)

# dom_pathways plots:
pl.plot_dom_pathways(skip_t=1, kk=0)



breakpoint()


