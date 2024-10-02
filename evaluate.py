import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import num2date, date2num
from plotting import PlotData

local_directory = 'C:/Users/ciank/PycharmProjects/sinmod/sg800_krill/'
save_folder = local_directory + 'figures/'
trajectory_folder = local_directory + 'trajectory/'
pl = PlotData(save_folder, trajectory_folder)


pl.read_trajectory_file(file_prefix='trajectory_d30')
pl.plot_trajectory_color(skip_n=5, skip_t=1)

#pl.plot_currents()
breakpoint()


