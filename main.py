from reader_SG800 import SGReader
from SG_krill import Krill
import datetime
import numpy as np
import matplotlib.pyplot as plt
samples_folder = 'E:/fromIngrid/samplefiles_reruns/'
samples_prefix = 'samplesNSEW_'
local_folder = 'C:/Users/ciank/PycharmProjects/sinmod/sg800_krill/'
trajectory_folder = local_folder + 'trajectory/'

counter = -1
range_vals = np.arange(0, 1, 0.001)
phi = np.zeros(range_vals.shape[0])
max_sp = 0.012
v_lower = 0.4
v_upper = v_lower + 0.2
mu, sigma = 0, 0.06
rv = np.random.normal(mu, sigma, range_vals.shape[0])
for x in range_vals:
    counter = counter + 1
    if (x >= v_lower) & (x <= v_upper):
        phi[counter] = 0
    elif (x > v_upper):
        phi[counter] = ((x - v_upper) ** 2) * rv[counter]
    elif (x < v_lower):
        phi[counter] = ((x - v_lower) ** 2) * rv[counter]

    if phi[counter] > max_sp:
        phi[counter] = max_sp
    if phi[counter] < -max_sp:
        phi[counter] = -max_sp


plt.plot(range_vals, phi)
plt.xlim([0, 1])
plt.show()
breakpoint()


# time parameters: #todo: have a date_init variable that calls the correct file
date_init = datetime.datetime(2017, 5, 31, 21, 0)
duration_days = datetime.timedelta(days=10)  # simulation duration in days;
minutes = 4
dt = datetime.timedelta(hours=minutes/60)
save_step = datetime.timedelta(hours=4)  # save time step
time_threshold = datetime.timedelta(seconds=0)
simulation_steps = duration_days/dt
save_number = duration_days/save_step

# model parameters:
n = 10000  # particle number
N = 1 # Ensemble member;
x_min = 210.0  # coordinates for initialization;
x_max = 610.4
y_min = 210.0
y_max = 610.5

# reader that adds info from current file
reader_SG = SGReader(samples_folder, samples_prefix, duration_days.days, date_init)

# start of simulation:
time_counter = save_step
save_counter = -1
for i in np.arange(0, simulation_steps, 1):
    time_counter -= dt
    if i == 0:
        k = Krill(reader_SG.nc_file)
        k.init_krill(N=N, n=n, x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max)
        k.init_netcdf(trajectory_folder=trajectory_folder, N=N, n=n,
                      save_n=save_number, save_file=reader_SG.save_file)

    k.step_krill(dt, reader_SG)
    #k.plot_init(kk=0) # plot an example ensemble member

    # update the current time in simulation
    reader_SG.update_time(dt)

    # save simulation
    if time_counter <= time_threshold:
        save_counter += 1
        time_counter = save_step
        print('saving ' + str(save_counter + 1) + ' of ' + str(save_number))
        k.save_step(save_counter=save_counter, current_datetime=reader_SG.current_datetime)

    if reader_SG.current_datetime.month != reader_SG.file_month:
        reader_SG.nc_file.close()
        reader_SG.update_file_month(samples_folder=samples_folder, samples_prefix=samples_prefix)

k.trajectory_file.close()
