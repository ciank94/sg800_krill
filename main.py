from reader_SG800 import SGReader
from SG_krill import Krill
import datetime
import numpy as np
import logging
import time
import matplotlib.pyplot as plt

# switches for changing the type of simulation:
remote = False  # remote or local server
test = True  # load data less frequently
bio_mapping = False  # if true do mapping between SINMOD and PISCES grid
temp_beh = False # add response to the temperature gradient
dvm_beh = False # move particles based on light conditions
feed_beh = True # particles have an f parameter defining feeding;

# time parameters:
date_init = datetime.datetime(2017, 3, 1, 1, 0)
duration_days = datetime.timedelta(days=10)  # simulation duration in days;
minutes = 4
dt = datetime.timedelta(hours=minutes/60)
save_step = datetime.timedelta(hours=4)  # save time step
time_threshold = datetime.timedelta(seconds=0)
simulation_steps = duration_days/dt
save_number = duration_days/save_step

# model parameters:
n = 400  # particle number
N = 4 # Ensemble member;
x_min = 210.0  # coordinates for initialization; todo: make initialization
x_max = 610.4
y_min = 210.0
y_max = 610.5

# Furnish reader with initial conditions:
reader_SG = SGReader(remote, test, bio_mapping, temp_beh, dvm_beh, feed_beh) # initialise switches for simulation
reader_SG.file_explorer()  # define paths to files in this folder
reader_SG.init_time(duration_days.days, date_init) # reader that adds time info from current file
reader_SG.read_bio() # read biology file
reader_SG.log_init(n, N, x_min, x_max, y_min, y_max, dt, save_step, simulation_steps, save_number) # log init to file

# start of simulation:
time_counter = save_step
save_counter = -1
start_time = time.time()
for i in np.arange(0, simulation_steps, 1):
    time_counter -= dt
    if i == 0:
        reader_SG.logger.info('Initialising krill')
        k = Krill(reader_SG.nc_file)
        k.init_temp_response(N=N)
        k.init_krill(N=N, n=n, x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max)
        k.init_netcdf(trajectory_folder=reader_SG.trajectory_folder, N=N, n=n,
                      save_n=save_number, save_file=reader_SG.save_file)


    k.step_krill(dt, reader_SG)
    #k.plot_init(kk=0) # plot an example ensemble member

    # update the current time in simulation
    reader_SG.update_time(dt)

    # save simulation
    if time_counter <= time_threshold:
        save_counter += 1
        time_counter = save_step
        reader_SG.logger.info('saving ' + str(save_counter + 1) + ' of ' + str(save_number))
        k.save_step(save_counter=save_counter, current_datetime=reader_SG.current_datetime)

python_time = time.time() - start_time
reader_SG.logger.info(f"Simulation took {python_time/(60*60):.2f} hours.")
k.trajectory_file.close()
