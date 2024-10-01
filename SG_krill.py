from netCDF4 import num2date
import numpy as np
import matplotlib.pyplot as plt


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

    def init_krill(self, n, x_min, x_max, y_min, y_max):
        self.x = np.zeros([n])
        self.y = np.zeros([n])
        step_xy = (n / np.sqrt(n))
        counter = -1
        for ii in np.arange(0, step_xy, 1):
            x = x_min + ((x_max - x_min)*(((ii + 1)/step_xy)))
            for jj in np.arange(0, step_xy, 1):
                y = y_min + ((y_max - y_min) * (((jj + 1)/step_xy)))
                if (x > x_max) | (y>y_max):
                    breakpoint()
                counter = counter + 1
                if self.depth[np.floor(y).astype(int), np.floor(x).astype(int)] > 0:
                    self.x[counter] = x
                    self.y[counter] = y
                else:
                    self.x[counter] = np.nan
                    self.y[counter] = np.nan

        print('Initialised krill with x and y values')
        return

    def step_krill(self, dt_datetime, reader_SG):
        time_index = reader_SG.time_index
        nc_file = reader_SG.nc_file
        dt = dt_datetime.seconds
        if time_index > self.load_counter:
            self.u_east = nc_file['u_east'][time_index, 0, :, :]
            self.v_north = nc_file['v_north'][time_index, 0, :, :]
            self.load_counter = time_index
            print('loading variables')
        for ii in range(0, self.x.shape[0]):
            if self.x[ii] > 0:
                u = self.u_east[np.floor(self.y[ii]).astype(int), np.floor(self.x[ii]).astype(int)]
                v = self.v_north[np.floor(self.y[ii]).astype(int), np.floor(self.x[ii]).astype(int)]
                self.x[ii] = self.x[ii] + ((dt*u)/self.res)
                self.y[ii] = self.y[ii] + ((dt*v)/self.res)
                if (self.x[ii] > self.i_max) | (self.y[ii] > self.j_max)|(self.y[ii] < 1)| (self.x[ii] < 1):
                    self.x[ii] = np.nan
                    self.y[ii] = np.nan
                elif np.isnan(self.depth[np.floor(self.y[ii]).astype(int), np.floor(self.x[ii]).astype(int)]):
                    self.x[ii] = np.nan
                    self.y[ii] = np.nan


        return

    def plot_init(self):
        fig, axs = plt.subplots(figsize=(12, 8), layout='constrained')
        axs.pcolormesh(self.depth, cmap=plt.get_cmap('copper'))
        axs.scatter(self.x, self.y, c='r', s=1)
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
        #axs.scatter(self.x, self.y, s=0.4, c='r')
        fig.colorbar(d_map)
        plt.show()
        return