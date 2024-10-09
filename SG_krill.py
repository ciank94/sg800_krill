from netCDF4 import num2date, date2num
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import logging
from numba import njit
import os


class Krill:
    def __init__(self, nc_file):
        self.i_max = nc_file['xc'].shape[0]
        self.j_max = nc_file['yc'].shape[0]
        self.lay_depths = np.cumsum(np.array(nc_file['LayerDepths'][:]))
        self.depth = np.array(nc_file['depth'])
        self.depth[self.depth>10000] = np.nan
        self.res = 800
        self.load_counter = -1
        self.init_temp = False
        self.logger = logging.getLogger(__name__)
        return

    def init_temp_response(self, N):
        sqN_int = np.sqrt(N).astype(int)
        t_scale = np.linspace(-1.5, 1.5, sqN_int)
        w_scale = np.linspace(1, np.sqrt(N), sqN_int)
        [t_vals, w_vals] = np.meshgrid(t_scale, w_scale)
        self.t_min = np.reshape(t_vals, [N])
        self.t_max = self.t_min + 3
        self.w_max = np.reshape(w_vals, [N])*0.06
        self.logger.warning('min temps: ' + str(self.t_min))
        self.logger.warning('max temps: ' + str(self.t_max))
        self.logger.warning('max w: ' + str(self.w_max))
        self.init_temp = True
        return

    def init_krill(self, N, n, x_min, x_max, y_min, y_max):
        self.x = np.zeros([n, N])
        self.y = np.zeros([n, N])
        self.z = np.ones([n, N])*75
        step_xy = (n / np.sqrt(n))

        for kk in range(0, N, 1):
            counter = -1
            for ii in np.arange(0, step_xy, 1):
                x = x_min + ((x_max - x_min) * (((ii + 1) / step_xy)))
                for jj in np.arange(0, step_xy, 1):
                    y = y_min + ((y_max - y_min) * (((jj + 1) / step_xy)))
                    counter = counter + 1
                    #if (x > x_max) | (y > y_max):
                        #breakpoint()
                    #if counter >= self.x.shape[0]:
                        #breakpoint()
                    if self.depth[np.floor(y).astype(int), np.floor(x).astype(int)] > 0:
                        self.x[counter, kk] = x
                        self.y[counter, kk] = y
                        self.z[counter, kk] = self.z[counter, kk]
                    else:
                        self.x[counter, kk] = np.nan
                        self.y[counter, kk] = np.nan
                        self.z[counter, kk] = np.nan
        print('Initialised krill with x and y values')
        return

    def step_krill_test(self, dt_datetime, reader_SG):
        time_index = reader_SG.time_index
        nc_file = reader_SG.nc_file
        dt = dt_datetime.seconds
        if time_index != self.load_counter:
            self.u_east = np.array(nc_file['u_east'][time_index, :, :, :])
            self.v_north = np.array(nc_file['v_north'][time_index, :, :, :])
            self.temp = np.array(nc_file['temperature'][time_index, :, :, :])
            self.temp[self.temp <- 3000]=np.nan
            self.v_north[self.v_north < -3000] = np.nan
            self.u_east[self.u_east < -3000] = np.nan
            self.load_counter = time_index
            print('loading variables')
        self.t_save = np.ones([self.x.shape[0], self.x.shape[1]]) * -32000
        self.w_save = np.ones([self.x.shape[0], self.x.shape[1]]) * -32000
        N = self.x.shape[1]
        n = self.x.shape[0]
        [x1, y1, z1, t_save, w_save] = loop_individuals(self.x, self.y, self.z, self.i_max, self.j_max,
                                                               self.u_east, self.v_north, self.temp, self.depth,
                                                               self.lay_depths, self.t_min, self.t_max, dt, self.res,
                                                               self.t_save, self.w_save, self.w_max, N, n)


        self.x = x1
        self.y = y1
        self.z = z1
        self.t_save = t_save
        self.w_save = w_save

        #breakpoint()

        return

    def step_krill(self, dt_datetime, reader_SG):
        time_index = reader_SG.time_index
        nc_file = reader_SG.nc_file
        dt = dt_datetime.seconds
        if time_index != self.load_counter:
            self.u_east = np.array(nc_file['u_east'][time_index, :, :, :])
            self.v_north = np.array(nc_file['v_north'][time_index, :, :, :])
            self.temp = np.array(nc_file['temperature'][time_index, :, :, :])
            self.temp[self.temp <- 3000]=np.nan
            self.v_north[self.v_north < -3000] = np.nan
            self.u_east[self.u_east < -3000] = np.nan
            self.load_counter = time_index
            print('loading variables')
        self.t_save = np.ones([self.x.shape[0], self.x.shape[1]]) * -32000
        self.w_save = np.ones([self.x.shape[0], self.x.shape[1]]) * -32000
        for kk in range(0, self.x.shape[1]):
            for ii in range(0, self.x.shape[0]):
                if self.x[ii, kk] > 0:
                    layer_ii = np.argmin((self.lay_depths - self.z[ii, kk])**2)
                    u = self.u_east[layer_ii, np.floor(self.y[ii, kk]).astype(int), np.floor(self.x[ii, kk]).astype(int)]
                    v = self.v_north[layer_ii, np.floor(self.y[ii, kk]).astype(int), np.floor(self.x[ii, kk]).astype(int)]
                    t = self.temp[layer_ii, np.floor(self.y[ii, kk]).astype(int), np.floor(self.x[ii, kk]).astype(int)] - 273.15
                    grads_t = np.gradient(self.temp[:, np.floor(self.y[ii, kk]).astype(int), np.floor(self.x[ii, kk]).astype(int)] - 273.15)
                    grad_t = grads_t[layer_ii]

                    if (np.isnan(t))|(np.isnan(u))|(np.isnan(v)):
                        self.x[ii, kk] = np.nan
                        self.y[ii, kk] = np.nan
                        self.z[ii, kk] = np.nan
                    elif (self.x[ii, kk] > self.i_max) | (self.y[ii, kk] > self.j_max) | (self.y[ii, kk] < 1)| (self.x[ii, kk] < 1):
                        self.x[ii, kk] = np.nan
                        self.y[ii, kk] = np.nan
                        self.z[ii, kk] = np.nan
                    elif np.isnan(self.depth[np.floor(self.y[ii, kk]).astype(int), np.floor(self.x[ii, kk]).astype(int)]):
                        self.x[ii, kk] = np.nan
                        self.y[ii, kk] = np.nan
                        self.z[ii, kk] = np.nan
                    else: # other options for krill
                        if (t > self.t_max[kk]) | (t < self.t_min[kk]):
                            w = self.swim_temp(t, kk, grad_t)
                        else:
                            w = 0.
                        self.x[ii, kk] = self.x[ii, kk] + ((dt * u) / self.res)
                        self.y[ii, kk] = self.y[ii, kk] + ((dt * v) / self.res)
                        self.z[ii, kk] = self.z[ii, kk] + ((dt * w) / self.res)
                        self.t_save[ii, kk] = t
                        self.w_save[ii, kk] = w
        return

    def swim_temp(self, t, kk, grad_t):
        if (grad_t == 0) | (np.isnan(grad_t)):
            scale_w = np.random.choice([-1, 1])
            w = scale_w * self.w_max[kk]
        elif t < self.t_min[kk]:
            #scale_w = ((t - self.t_min) ** 2) / (2 ** 2)
            w = np.sign(grad_t) * self.w_max[kk]
        else:
            #scale_w = (((t - self.t_max) ** 2) / (2 ** 2))
            w = -1 * np.sign(grad_t) * self.w_max[kk]
        return w

    def save_step(self, save_counter, current_datetime):
        self.trajectory_file['xp'][:, :, save_counter] = self.x
        self.trajectory_file['yp'][:, :, save_counter] = self.y
        self.trajectory_file['zp'][:, :, save_counter] = self.z
        self.trajectory_file['w'][:, :, save_counter] = self.w_save
        self.trajectory_file['temp'][:, :, save_counter] = self.t_save
        self.trajectory_file['time'][save_counter] = date2num(current_datetime, self.trajectory_file['time'].unit)
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
                             'zp': {'datatype': 'f4', 'dimensions': ('trajectory', 'ensemble', 'obs'),
                                    'description': 'z position of particle'},
                             'temp': {'datatype': 'f4', 'dimensions': ('trajectory', 'ensemble', 'obs'),
                                    'description': 'temperature of particle'},
                             'w': {'datatype': 'f4', 'dimensions': ('trajectory', 'ensemble', 'obs'),
                                      'description': 'speed of particle'},
                             'time': {'datatype': 'f4', 'dimensions': ('obs',),
                                      'description': 'datetime of particle'},
                             't_min':{'datatype': 'f4', 'dimensions': ('ensemble',),
                                      'description': 'datetime of particle'},
                             't_max': {'datatype': 'f4', 'dimensions': ('ensemble',),
                                       'description': 'datetime of particle'},
                             'w_max': {'datatype': 'f4', 'dimensions': ('ensemble',),
                                       'description': 'datetime of particle'}
                             }

        for variable in variable_key_dict:
            self.trajectory_file.createVariable(variable, variable_key_dict[variable]['datatype'],
                                                variable_key_dict[variable]['dimensions'])
            self.trajectory_file[variable].description = variable_key_dict[variable]['description']

        if self.init_temp:
            self.trajectory_file['t_min'][:] = self.t_min
            self.trajectory_file['t_max'][:] = self.t_max
            self.trajectory_file['w_max'][:] = self.w_max
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

@njit
def loop_individuals(x1, y, z, i_max, j_max, u_east, v_north, temp, depth, lay_depths, t_min, t_max,
                     dt, res, t_save, w_save, w_max, N, n):
    choice_v = np.array([-1, 1])
    for kk in range(0, N, 1):
        for ii in range(0, n, 1):
            if x1[ii, kk] > 0:
                x_int = round(x1[ii, kk])
                y_int = round(y[ii, kk])
                layer_ii = np.argmin((lay_depths - z[ii, kk]) ** 2)
                ue = u_east[layer_ii, :, :]
                vv = v_north[layer_ii, :, :]
                tt = temp[layer_ii, :, :] - 273.15
                t = tt[y_int, x_int]
                u = ue[y_int, x_int]
                v = vv[y_int, x_int]
                temp_list = temp[:, y_int, x_int] - 273.15
                grads_t = np.zeros(temp_list.shape[0])
                grads_t[1:None] = temp_list[1:None] - temp_list[0:-1]
                #grads_t = np.gradient(temp[:, y_int, x_int] - 273.15)
                grad_t = grads_t[layer_ii]

                if (np.isnan(t)) | (np.isnan(u)) | (np.isnan(v)):
                    x1[ii, kk] = np.nan
                    y[ii, kk] = np.nan
                    z[ii, kk] = np.nan
                elif (x1[ii, kk] > i_max) | (y[ii, kk] > j_max) | (y[ii, kk] < 1) | (
                        x1[ii, kk] < 1):
                    x1[ii, kk] = np.nan
                    y[ii, kk] = np.nan
                    z[ii, kk] = np.nan
                elif np.isnan(depth[y_int, x_int]):
                    x1[ii, kk] = np.nan
                    y[ii, kk] = np.nan
                    z[ii, kk] = np.nan
                else:  # other options for krill
                    if (t > t_max[kk]) | (t < t_min[kk]):
                        if (grad_t == 0) | (np.isnan(grad_t)):
                            scale_w = np.random.choice(choice_v)
                            w = scale_w * w_max[kk]
                        elif t < t_min[kk]:
                            # scale_w = ((t - self.t_min) ** 2) / (2 ** 2)
                            w = np.sign(grad_t) * w_max[kk]
                        else:
                            # scale_w = (((t - self.t_max) ** 2) / (2 ** 2))
                            w = -1 * np.sign(grad_t) * w_max[kk]
                    else:
                        w = 0.

                    x1[ii, kk] = x1[ii, kk] + ((dt * u) / res)
                    y[ii, kk] = y[ii, kk] + ((dt * v) / res)
                    z[ii, kk] = z[ii, kk] + ((dt * w) / res)
                    t_save[ii, kk] = t
                    w_save[ii, kk] = w

    return x1, y, z, t_save, w_save

