import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import cmasher as cmr
from netCDF4 import num2date, date2num

class PlotData:
    def __init__(self, save_folder, trajectory_folder):
        # bunch of parameters and so on
        self.save_folder = save_folder
        self.trajectory_folder = trajectory_folder
        self.depth_file = self.save_folder + 'depth.npy'
        self.u_file = self.save_folder + 'u_east.npy'
        self.v_file = self.save_folder + 'v_north.npy'
        try: self.depth = np.load(self.depth_file)
        except: self.init_depth_file()
        self.i_max = self.depth.shape[1]
        self.j_max = self.depth.shape[0]
        try: self.u_east = np.load(self.u_file)
        except: self.init_current_file()
        self.v_north = np.load(self.v_file)
        self.bath_contours = np.linspace(0, 3000, 10)
        return


    def read_trajectory_file(self, file_prefix):
        trajectory_file = self.trajectory_folder + file_prefix + '.nc'
        self.df = nc.Dataset(trajectory_file, mode='r')
        return


    def plot_trajectory_color(self, skip_n, skip_t):
        fig, axs = plt.subplots(figsize=(12, 8), layout='constrained')
        #cmap = plt.get_cmap('Blues')
        cmap = cmr.ocean_r
        cmap.set_bad('gray', 1.)
        d_map = axs.pcolormesh(self.depth, cmap=cmap, alpha=0.75, vmax=4000)
        axs.contour(self.depth, levels=self.bath_contours, colors='k', alpha=0.35,
                    linewidths=1.5)
        x_times = self.df['xp'][:, ::skip_t]
        y_times = self.df['yp'][:, ::skip_t]
        #time_stamps = num2date(self.df['time'][::skip_t], self.df['time'].unit)
        t_vals = np.arange(0,skip_t*(4/24)*x_times.shape[1],skip_t*(4/24))


        t = np.arange(x_times.shape[1])
        for i in range(0, self.df['xp'].shape[0], skip_n):
            x = x_times[i]
            y = y_times[i]
            axs.plot(x, y, 'k', alpha=0.8,zorder=1, linewidth=1)
            axs.scatter(x, y, c=t, s=2.3, cmap=cmr.ember,zorder=2)
        scatter_map = axs.scatter(x, y, c=t_vals, s=2, cmap=cmr.ember, zorder=2)
        cbar1 = fig.colorbar(d_map)
        cbar1.ax.tick_params(labelsize=12)
        cbar1.ax.set_ylabel('depth (m)', loc='center', size=12, weight='bold')
        cbar=fig.colorbar(scatter_map)
        cbar.ax.tick_params(labelsize=12)
        cbar.ax.set_ylabel('time elapsed (days)', loc='center', size=12, weight='bold')
        self.save_plot(plt_name='trajectory_color')
        return

    def plot_currents(self):
        fig, axs = plt.subplots(figsize=(12, 8), layout='constrained')
        mag_uv0 = np.sqrt(np.square(self.u_east) + np.square(self.v_north))
        cmap = plt.get_cmap('viridis')
        cmap.set_bad(color='gray', alpha=0.6)
        d_map = axs.pcolormesh(mag_uv0, cmap=cmap,vmax=0.5)
        axs.contour(self.depth, levels=self.bath_contours, colors='k', alpha=0.35,
                        linewidths=1.3)
        sk = 12
        xx = np.arange(0, self.i_max,1)
        yy = np.arange(0, self.j_max,1)
        [ux, yx] = np.meshgrid(xx,yy)
        axs.quiver(ux[::sk,::sk], yx[::sk,::sk], self.u_east[::sk,::sk],self.v_north[::sk,::sk], edgecolors='k', alpha=0.7, linewidths=0.1)
        #axs.scatter(self.x, self.y, s=0.4, c='r')
        fig.colorbar(d_map)
        self.save_plot(plt_name='current_map')
        return

    def init_current_file(self):
        print('creating current file')
        samples_folder = 'E:/fromIngrid/samplefiles_reruns/'
        samples_prefix = 'samplesNSEW_'
        init_year = 2017
        init_month = 2
        filename = samples_folder + samples_prefix + str(init_year) + "{:02d}".format(init_month) + '.nc'
        df = nc.Dataset(filename)
        self.u_east = np.array(np.nanmean(df['u_east'][0:10,0,:,:], 0))
        self.v_north= np.array(np.nanmean(df['v_north'][0:10,0,:,:],0))
        #self.u_east[self.u_east<-1000] = np.nan
        #self.v_north[self.v_north < -1000] = np.nan
        np.save(self.u_file, self.u_east)
        np.save(self.v_file, self.v_north)
        return


    def init_depth_file(self):
        print('creating depth file')
        samples_folder = 'E:/fromIngrid/samplefiles_reruns/'
        samples_prefix = 'samplesNSEW_'
        init_year = 2017
        init_month = 2
        filename = samples_folder + samples_prefix + str(init_year) + "{:02d}".format(init_month) + '.nc'
        df = nc.Dataset(filename)
        self.depth = np.array(df['depth'])
        self.depth[self.depth > 10000] = np.nan
        np.save(self.depth_file, self.depth)
        return

    def save_plot(self, plt_name):
        savefile = self.save_folder + plt_name + '.png'
        print('Saving file: ' + savefile)
        plt.savefig(savefile, dpi=400)
        plt.close()
        return
