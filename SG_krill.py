from netCDF4 import num2date, date2num
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import os


class Krill:
    def __init__(self, nc_file):
        self.i_max = nc_file['xc'].shape[0]
        self.j_max = nc_file['yc'].shape[0]
        self.lay_depths = nc_file['LayerDepths']
        self.depth = np.array(nc_file['depth'])
        self.depth[self.depth>10000] = np.nan
        self.res = 800
        self.load_counter = -1
        return

    def init_krill(self, N, n, x_min, x_max, y_min, y_max):
        self.x = np.zeros([n, N])
        self.y = np.zeros([n, N])
        step_xy = (n / np.sqrt(n))

        for kk in range(0, N, 1):
            counter = -1
            for ii in np.arange(0, step_xy, 1):
                x = x_min + ((x_max - x_min) * (((ii + 1) / step_xy)))
                for jj in np.arange(0, step_xy, 1):
                    y = y_min + ((y_max - y_min) * (((jj + 1) / step_xy)))
                    counter = counter + 1
                    if (x > x_max) | (y > y_max):
                        breakpoint()
                    if counter >= self.x.shape[0]:
                        breakpoint()
                    if self.depth[np.floor(y).astype(int), np.floor(x).astype(int)] > 0:
                        self.x[counter, kk] = x
                        self.y[counter, kk] = y
                    else:
                        self.x[counter, kk] = np.nan
                        self.y[counter, kk] = np.nan
        print('Initialised krill with x and y values')
        return

    def step_krill(self, dt_datetime, reader_SG):
        time_index = reader_SG.time_index
        nc_file = reader_SG.nc_file
        dt = dt_datetime.seconds
        if time_index > self.load_counter:
            self.u_east = nc_file['u_east'][time_index, :, :, :]
            self.v_north = nc_file['v_north'][time_index, :, :, :]
            self.load_counter = time_index
            print('loading variables')
        layer_ii = 0
        for kk in range(0, self.x.shape[1]):
            for ii in range(0, self.x.shape[0]):
                if self.x[ii, kk] > 0:
                    u = self.u_east[layer_ii, np.floor(self.y[ii, kk]).astype(int), np.floor(self.x[ii, kk]).astype(int)]
                    v = self.v_north[layer_ii, np.floor(self.y[ii, kk]).astype(int), np.floor(self.x[ii, kk]).astype(int)]
                    self.x[ii, kk] = self.x[ii, kk] + ((dt * u) / self.res)
                    self.y[ii, kk] = self.y[ii, kk] + ((dt * v) / self.res)
                    if ((self.x[ii, kk] > self.i_max) | (self.y[ii, kk] > self.j_max) | (self.y[ii, kk] < 1)
                            | (self.x[ii, kk] < 1)):
                        self.x[ii, kk] = np.nan
                        self.y[ii, kk] = np.nan
                    elif (np.isnan
                        (self.depth[np.floor(self.y[ii, kk]).astype(int), np.floor(self.x[ii, kk]).astype(int)])):
                        self.x[ii, kk] = np.nan
                        self.y[ii, kk] = np.nan

        return

    def save_step(self, save_counter, current_datetime):
        self.trajectory_file['xp'][:,:,save_counter] = self.x
        self.trajectory_file['yp'][:,:,save_counter] = self.y
        self.trajectory_file['time'][save_counter] = date2num(current_datetime,self.trajectory_file['time'].unit)
        return

    def init_netcdf(self, trajectory_folder, N, n, save_n, save_file):
        trajectory_filename = trajectory_folder + save_file

        self.trajectory_file = nc.Dataset(trajectory_filename, mode='w')
        dimension_key_dict = {'trajectory': n, 'ensemble': N, 'obs': save_n}

        for dimension in dimension_key_dict:
            self.trajectory_file.createDimension(dimension, dimension_key_dict[dimension])

        variable_key_dict = {'xp': {'datatype': 'f4', 'dimensions': ('trajectory', 'ensemble', 'obs'),
                                           'description': 'x position of particle'},
                             'yp': {'datatype': 'f4', 'dimensions': ('trajectory', 'ensemble', 'obs'),
                                    'description': 'y position of particle'},
                             'time': {'datatype': 'f4', 'dimensions': ('obs',),
                                      'description': 'datetime of particle'}
                             }

        for variable in variable_key_dict:
            self.trajectory_file.createVariable(variable, variable_key_dict[variable]['datatype'],
                                                variable_key_dict[variable]['dimensions'])
            self.trajectory_file[variable].description = variable_key_dict[variable]['description']

        time_unit_out = "seconds since 2014-04-01 00:00:00"
        self.trajectory_file['time'].setncattr('unit', time_unit_out)
        print('Initialising trajectory file: ' + trajectory_filename)
        return


    def plot_init(self, kk):
        fig, axs = plt.subplots(figsize=(12, 8), layout='constrained')
        axs.pcolormesh(self.depth, cmap=plt.get_cmap('viridis'))
        axs.scatter(self.x[:,kk], self.y[:,kk], c='r', s=1)
        plt.show()
        return

    def plot_currents(self):
        fig, axs = plt.subplots(figsize=(12, 8), layout='constrained')
        mag_uv0 = np.sqrt(np.square(self.u_east) + np.square(self.v_north))
        d_map = axs.pcolormesh(mag_uv0, cmap=plt.get_cmap('viridis'), vmin=0)
        sk = 10
        xx = np.arange(0, self.i_max,1)
        yy = np.arange(0, self.j_max,1)
        [ux, yx] = np.meshgrid(xx,yy)
        axs.quiver(ux[::sk,::sk], yx[::sk,::sk], self.u_east[::sk,::sk],self.v_north[::sk,::sk], edgecolors='k', alpha=0.5, linewidths=0.05)
        axs.scatter(self.x, self.y, s=0.4, c='r')
        fig.colorbar(d_map)
        plt.show()
        return