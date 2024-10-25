from reader_SG800 import SGReader
from SG_krill import Krill
import datetime
import numpy as np
import logging
import time
import matplotlib.pyplot as plt

# add a tag for this experiment:
exp_tag = 'sens2dep'

# switches for changing the type of simulation:
remote = False  # remote or local server
test = True  # load data less frequently
light_mapping = False # mapping between SINMOD and ERA5 grid
bio_mapping = False  # mapping between SINMOD and PISCES grid
temp_beh = False  # vertical response to temperature gradient
dvm_beh = False  # vertical response to light conditions (dvm_beh overrides temp_beh)
feed_beh = True  # f parameter defining feeding

# time parameters:
date_init = datetime.datetime(2017, 1, 1, 1, 0)
duration_days = datetime.timedelta(days=40)  # simulation duration in days;
minutes = 120
dt = datetime.timedelta(hours=minutes/60)
save_step = datetime.timedelta(hours=4)  # save time step
time_threshold = datetime.timedelta(seconds=0)
simulation_steps = duration_days/dt
save_number = duration_days/save_step

# model parameters:
n = 1600  # particle number
N = 1  # ensemble members;
x_min = 210.0  # coordinates for initialization; todo: make initialization
x_max = 610.4
y_min = 210.0
y_max = 610.5

# Furnish reader with initial parameters:
reader_SG = SGReader(remote, test, light_mapping, bio_mapping, temp_beh, dvm_beh, feed_beh, exp_tag)  # initialise switches for simulation
reader_SG.file_explorer()  # define paths to files in this folder
reader_SG.init_time(duration_days.days, date_init)  # reader that adds time info from current file
reader_SG.read_bio()  # read biology file
reader_SG.read_light()
reader_SG.log_init(n, N, x_min, x_max, y_min, y_max, dt, save_step, simulation_steps, save_number)  # log init to file

# start of simulation:
time_counter = save_step
save_counter = -1
start_time = time.time()
for i in np.arange(0, simulation_steps, 1):
    time_counter -= dt
    if i == 0:
        reader_SG.logger.info('Initialising krill')
        k = Krill(reader_SG.nc_file)
        if reader_SG.temp_beh:
            k.init_temp_response(reader_SG)
        k.init_krill(reader_SG)
        k.init_netcdf(reader_SG)

    # function that steps krill forward in time
    k.step_krill(dt, reader_SG)

    # update the current time in simulation
    reader_SG.update_time(dt)

    # save simulation
    if time_counter <= time_threshold:
        save_counter += 1
        time_counter = save_step
        reader_SG.logger.info('saving ' + str(save_counter + 1) + ' of ' + str(reader_SG.save_number))
        k.save_step(save_counter, reader_SG)
        # k.plot_init(kk=0)  # plot an example ensemble member

# finished simulation;
python_time = time.time() - start_time
reader_SG.logger.info(f"Simulation took {python_time/(60*60):.2f} hours.")
k.trajectory_file.close()
