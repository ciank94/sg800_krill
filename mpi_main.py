from reader_SG800 import SGReader
from SG_krill import Krill
import datetime
import numpy as np
import logging
import time
import matplotlib.pyplot as plt
from mpi4py import MPI


world_comm = MPI.COMM_WORLD
world_size = world_comm.Get_size()
my_rank = world_comm.Get_rank()

# namelist for experiment setup:
#namelist_path = '/cluster/projects/nn9828k/Cian_sinmod/sg800_krill/'
#namelist_path = 'C:/Users/ciank/PycharmProjects/sinmod/sg800_krill/'
namelist_path = '/cluster/work/ciank/sg800_krill/'

if my_rank == 0:
    release_day_increment = 0
else:
    release_day_increment = my_rank*5

# initialise reader based on namelist configuration
reader_SG = SGReader(namelist_path, release_day_increment)

# process:
reader_SG.logger.warning(str(reader_SG.init_datetime) + ' for my_rank = ' + str(my_rank))

# start of simulation:
time_counter = reader_SG.save_step
save_counter = -1
start_time = time.time()
for i in np.arange(0, reader_SG.simulation_steps, 1):
    time_counter -= reader_SG.dt
    if i == 0:
        reader_SG.logger.info('Initialising krill')
        k = Krill(reader_SG.nc_file)
        if reader_SG.temp_beh:
            k.init_temp_response(reader_SG)
        k.init_krill(reader_SG)
        k.init_netcdf(reader_SG)

    # function that steps krill forward in time
    k.step_krill(reader_SG)

    # update the current time in simulation
    reader_SG.update_time()

    # save simulation
    if time_counter <= reader_SG.time_threshold:
        save_counter += 1
        time_counter = reader_SG.save_step
        reader_SG.logger.info('saving ' + str(save_counter + 1) + ' of ' + str(reader_SG.save_number))
        k.save_step(save_counter, reader_SG)
        # k.plot_init(kk=0)  # plot an example ensemble member

# finished simulation;
k.trajectory_file.close()
python_time = time.time() - start_time
reader_SG.logger.info(f"Simulation took {python_time/(60*60):.2f} hours.")

