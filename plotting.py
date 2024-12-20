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
        self.depth_file = self.save_folder + "depth.npy"
        self.u_file = self.save_folder + "u_east.npy"
        self.v_file = self.save_folder + "v_north.npy"
        try:
            self.depth = np.load(self.depth_file)
        except:
            self.init_depth_file()
        self.i_max = self.depth.shape[1]
        self.j_max = self.depth.shape[0]
        try:
            self.u_east = np.load(self.u_file)
        except:
            self.init_current_file()
        self.v_north = np.load(self.v_file)
        self.bath_contours = np.linspace(0, 3000, 10)
        return

    def read_trajectory_file(self, file_prefix):
        trajectory_file = self.trajectory_folder + file_prefix + ".nc"
        self.df = nc.Dataset(trajectory_file, mode="r")
        self.dates = num2date(self.df["time"], self.df["time"].unit)
        # [np.sum(np.isnan(self.df['xp'][:, 0, i])) for i in range(0, self.df['xp'].shape[2])]
        return

    def entry_area(self, kk):
        x_min = 500
        x_max = 630
        x_slice = slice(x_min, x_max, 1)
        y_min = 475
        y_max = 625
        y_slice = slice(y_min, y_max, 1)
        skip_t = 1
        # kk=0
        fig, axs = plt.subplots(1, 1)
        counter = -1

        x_times = np.array(self.df["xp"][:, kk, ::skip_t])
        y_times = np.array(self.df["yp"][:, kk, ::skip_t])
        # exit_ind = np.zeros(x_times.shape[0])
        # for ii in range(0, x_times.shape[0]):
        #     id_ret = (x_times[ii, :] >= x_min) & (x_times[ii, :] <= x_max) & (y_times[ii, :] >= y_min) & (
        #             y_times[ii, :] <= y_max)
        #     if np.sum(id_ret) > 0:
        #         exit_ind[ii] = np.where(id_ret)[0][0]
        #
        #     print('n individuals leaving' + str(np.sum(exit_ind > 0)))

        exit_time = np.zeros(x_times.shape[1])
        for jj in range(0, x_times.shape[1]):
            id_ret = (
                (x_times[:, jj] >= x_min)
                & (x_times[:, jj] <= x_max)
                & (y_times[:, jj] >= y_min)
                & (y_times[:, jj] <= y_max)
            )
            if jj > 0:
                exit_time[jj] = np.sum(id_ret)
            x_times[id_ret, :] = np.nan
            y_times[id_ret, :] = np.nan

        axs.plot(np.cumsum(exit_time))
        axs.set_title("Initialised month = " + str(self.dates[0].month), fontsize=25)
        axs.set_xlabel("time step")
        axs.set_ylabel("cumulative number of individuals entering")
        plt.show()
        breakpoint()

    def plot_retention(self):
        x_min = 210.0  # coordinates for initialization; todo: make initialization
        x_max = 610.4
        y_min = 210.0
        y_max = 610.5
        skip_t = 1
        fig, axs = plt.subplots(1, 1)
        counter = -1
        clrs = ["r", "k", "b"]
        for kk in [2, 0, 3, 6]:
            counter = 0
            x_times = np.array(self.df["xp"][:, kk, ::skip_t])
            y_times = np.array(self.df["yp"][:, kk, ::skip_t])
            exit_ind = np.zeros(x_times.shape[0])
            for ii in range(0, x_times.shape[0]):
                id_ret = (
                    (x_times[ii, :] < x_min)
                    | (x_times[ii, :] > x_max)
                    | (y_times[ii, :] > y_max)
                    | (y_times[ii, :] < y_min)
                )
                if np.sum(id_ret) > 0:
                    exit_ind[ii] = np.where(id_ret)[0][0]

            print("n individuals leaving" + str(np.sum(exit_ind > 0)))

            exit_time = np.zeros(x_times.shape[1])
            for jj in range(0, x_times.shape[1]):
                id_ret = (
                    (x_times[:, jj] < x_min)
                    | (x_times[:, jj] > x_max)
                    | (y_times[:, jj] > y_max)
                    | (y_times[:, jj] < y_min)
                )
                exit_time[jj] = np.sum(id_ret)
                x_times[id_ret, :] = np.nan
                y_times[id_ret, :] = np.nan

            axs.plot(np.cumsum(exit_time))

        axs.set_xlabel("time step")
        axs.set_ylabel("n of individuals leaving")

        plt.legend(["w=0", "w=0.06", "w=0.12", "w=0.18"])
        plt.show()

        breakpoint()

    def plot_temp(self):
        fig, axs = plt.subplots()
        for kk in [2, 0, 3, 6]:
            t = np.array(self.df["w"][:, kk, :])
            t[t < -3000] = np.nan
            axs.plot(np.nanmean(t, 0))

        axs.set_xlabel("time step")
        axs.set_ylabel("w (m/s)")
        plt.legend(["w=0", "w=0.06", "w=0.12", "w=0.18"])
        plt.show()
        breakpoint()
        ff = 1
        fig, axs = plt.subplots(4, 1)
        counter = -1
        for kk in [2, 0, 3, 6]:
            counter += 1
            w_save = np.array(self.df["w"][:, kk, :])
            w_save[w_save < -3000] = np.nan
            axs[counter].hist(w_save.flatten(), bins=50)

        plt.show()

        fig, axs = plt.subplots(4, 1)
        counter = -1
        for kk in [2, 0, 3, 6]:
            counter += 1
            t = np.array(self.df["temp"][:, kk, :])
            t[t < -3000] = np.nan
            axs[counter].hist(t.flatten(), bins=50)

        # axs.set_xlabel('temp (C)')
        # axs.set_ylabel('frequency')
        # plt.legend(['w=0', 'w=0.06', 'w=0.12', 'w=0.18'])
        plt.show()
        breakpoint()

        fig, axs = plt.subplots()
        for kk in [2, 0, 3, 6]:
            t = np.array(self.df["temp"][:, kk, :])
            t[t < -3000] = np.nan
            axs.plot(np.nanmean(t, 0))

        axs.set_xlabel("time step")
        axs.set_ylabel("temp (C)")
        plt.legend(["w=0", "w=0.06", "w=0.12", "w=0.18"])
        plt.show()

        breakpoint()
        # np.array(self.df['yp'][:, kk, ::skip_t])

    def plot_depths(self, skip_t, kk):
        # x_times = np.array(self.df['xp'][:, kk, ::skip_t])
        # y_times = np.array(self.df['yp'][:, kk, ::skip_t])
        z_times = np.array(self.df["zp"][:, kk, ::skip_t])
        z_vals = np.arange(0, 400, 10)
        z_paths = np.zeros([z_vals.shape[0], z_times.shape[1]])
        id_z = np.digitize(z_times[:, :], z_vals)
        mx_idz = np.max(id_z)
        for i in range(0, id_z.shape[0]):
            for j in range(0, id_z.shape[1]):
                idx_p = id_z[i, j] < mx_idz
                if idx_p:
                    z_paths[id_z[i, j], j] += 1

        z_paths = (z_paths / (z_times.shape[0])) * 100
        fig, axs = plt.subplots()
        z_map = axs.pcolormesh(
            np.arange(0, z_paths.shape[1]),
            z_vals[:],
            z_paths,
            cmap=plt.get_cmap("hot"),
            vmax=30,
        )
        plt.gca().invert_yaxis()
        fig.colorbar(z_map)
        plt.show()
        breakpoint()
        # t_times = np.array(self.df['temp'][:, kk, ::skip_t])
        # w_times = np.array(self.df['w'][:, kk, ::skip_t])
        # t_times[t_times<-2000] = np.nan
        # w_times[w_times < -2000] = np.nan
        # plt.scatter(t_times, w_times)
        # [plt.plot(z_times[i,:],'b') for i in range(0, 1600)]
        # plt.show()
        breakpoint()

    def plot_depth_region(self, skip_t, kk):
        x_times = np.array(self.df["xp"][:, kk, ::skip_t])
        y_times = np.array(self.df["yp"][:, kk, ::skip_t])
        z_times = np.array(self.df["zp"][:, kk, ::skip_t])
        x_times[x_times > 1000] = np.nan
        y_times[y_times > 1000] = np.nan
        z_times[z_times > 1000] = np.nan
        dom_vals = np.zeros(self.depth.shape)
        dom_vals[:] = self.depth[:]
        dom_vals[~np.isnan(dom_vals)] = 0
        dom_vals2 = np.zeros(self.depth.shape)
        dom_vals2[:] = self.depth[:]
        dom_vals2[~np.isnan(dom_vals)] = 0

        for p in range(0, x_times.shape[0]):
            xp = x_times[p, :]
            yp = y_times[p, :]
            id_p = ~np.isnan(xp)
            dom_vals[yp[id_p].astype(int), xp[id_p].astype(int)] += 1
            dom_vals2[yp[id_p].astype(int), xp[id_p].astype(int)] += z_times[p, id_p]

        dom_vals2[dom_vals2 > 0] /= dom_vals[dom_vals2 > 0]
        print(str(x_times.shape[1]))
        # dom_vals[(dom_vals==0)] = np.nan
        fig, axs = plt.subplots(figsize=(12, 8), layout="constrained")
        cmap = plt.get_cmap("hot")
        cmap.set_bad("gray", 0.4)
        vmax = np.nanmax(dom_vals2)  # /2
        # vmax=0.0035
        dom_map = axs.pcolormesh(dom_vals2, cmap=cmap, alpha=1, vmax=vmax)
        axs.contour(
            self.depth,
            levels=self.bath_contours,
            colors="k",
            alpha=0.8,
            linewidths=1.5,
            zorder=2,
        )
        cbar = fig.colorbar(dom_map)
        cbar.ax.tick_params(labelsize=12)
        cbar.ax.set_ylabel("depth (m)", loc="center", size=12, weight="bold")
        self.save_plot(plt_name="depth_region")
        breakpoint()
        return

    def plot_anomalies(self, skip_t, kk):
        x_times = np.array(self.df["xp"][:, 2, ::skip_t])
        y_times = np.array(self.df["yp"][:, 2, ::skip_t])
        x_times[x_times > 1000] = np.nan
        y_times[y_times > 1000] = np.nan
        ref_dom_vals = np.zeros(self.depth.shape)
        ref_dom_vals[:] = self.depth[:]
        ref_dom_vals[~np.isnan(ref_dom_vals)] = 0
        for p in range(0, x_times.shape[0]):
            xp = x_times[p, :]
            yp = y_times[p, :]
            id_p = ~np.isnan(xp)
            ref_dom_vals[yp[id_p].astype(int), xp[id_p].astype(int)] += 1

        ref_dom_vals[~np.isnan(ref_dom_vals)] /= x_times.shape[0] * x_times.shape[1]

        x_times = np.array(self.df["xp"][:, kk, ::skip_t])
        y_times = np.array(self.df["yp"][:, kk, ::skip_t])
        x_times[x_times > 1000] = np.nan
        y_times[y_times > 1000] = np.nan
        dom_vals = np.zeros(self.depth.shape)
        dom_vals[:] = self.depth[:]
        dom_vals[~np.isnan(dom_vals)] = 0
        for p in range(0, x_times.shape[0]):
            xp = x_times[p, :]
            yp = y_times[p, :]
            id_p = ~np.isnan(xp)
            dom_vals[yp[id_p].astype(int), xp[id_p].astype(int)] += 1

        dom_vals[~np.isnan(dom_vals)] /= x_times.shape[0] * x_times.shape[1]

        anom_vals = dom_vals - ref_dom_vals
        # dom_vals[(dom_vals==0)] = np.nan
        fig, axs = plt.subplots(figsize=(12, 8), layout="constrained")
        cmap = plt.get_cmap("bwr")
        cmap.set_bad("gray", 0.4)
        vmax = np.nanmax(anom_vals) / 2
        # vmax=0.0035
        dom_map = axs.pcolormesh(anom_vals, cmap=cmap, alpha=1)
        axs.contour(
            self.depth,
            levels=self.bath_contours,
            colors="k",
            alpha=0.8,
            linewidths=1.5,
            zorder=2,
        )
        cbar = fig.colorbar(dom_map)
        cbar.ax.tick_params(labelsize=12)
        cbar.ax.set_ylabel("probability (%)", loc="center", size=12, weight="bold")
        self.save_plot(plt_name="anom_paths")
        return

    def plot_growth(self, skip_t, kk):
        x_times = np.array(self.df["xp"][:, kk, ::skip_t])
        y_times = np.array(self.df["yp"][:, kk, ::skip_t])
        length = np.array(self.df["lp"][:, kk, ::skip_t])
        x_times[x_times > 1000] = np.nan
        y_times[y_times > 1000] = np.nan
        length[length < -1000] = np.nan
        ldiff = (length[:, 1::] - length[:, 0:-1]) * 24
        growth = np.c_[np.zeros(ldiff.shape[0]), ldiff]
        dom_vals = np.zeros(self.depth.shape)
        dom_vals[:] = self.depth[:]
        dom_vals[~np.isnan(dom_vals)] = 0
        dom_vals_count = np.zeros(dom_vals.shape)
        dom_vals_count[:] = self.depth[:]
        dom_vals_count[~np.isnan(dom_vals)] = 0
        x_times[x_times >= self.i_max] = self.i_max - 1
        y_times[y_times >= self.j_max] = self.j_max - 1
        for p in range(0, x_times.shape[0]):
            xp = x_times[p, :]
            yp = y_times[p, :]
            id_p = (~np.isnan(xp)) & (~np.isnan(growth[p, :]))
            dom_vals[yp[id_p].astype(int), xp[id_p].astype(int)] += growth[p, id_p]
            dom_vals_count[yp[id_p].astype(int), xp[id_p].astype(int)] += 1

        dom_vals[(dom_vals > 0) & (dom_vals_count > 0)] /= dom_vals_count[
            (dom_vals > 0) & (dom_vals_count > 0)
        ]
        # print(str(x_times.shape[1]))
        # dom_vals[(dom_vals==0)] = np.nan
        fig, axs = plt.subplots(figsize=(12, 8), layout="constrained")
        cmap = cmr.nuclear
        cmap.set_bad("gray", 0.4)
        vmax = np.nanmax(dom_vals) / 1.1
        # vmax=0.0035
        dom_map = axs.pcolormesh(dom_vals, cmap=cmap, alpha=1, vmax=vmax)
        axs.contour(
            self.depth,
            levels=self.bath_contours,
            colors="k",
            alpha=0.8,
            linewidths=1.5,
            zorder=2,
        )
        cbar = fig.colorbar(dom_map)
        cbar.ax.tick_params(labelsize=12)
        cbar.ax.set_ylabel("growth (mm day**-1)", loc="center", size=12, weight="bold")
        self.save_plot(plt_name="growth")
        return

    def plot_dom_pathways(self, skip_t, kk):
        x_times = np.array(self.df["xp"][:, kk, ::skip_t])
        y_times = np.array(self.df["yp"][:, kk, ::skip_t])
        x_times[x_times > 1000] = np.nan
        y_times[y_times > 1000] = np.nan
        dom_vals = np.zeros(self.depth.shape)
        dom_vals[:] = self.depth[:]
        dom_vals[~np.isnan(dom_vals)] = 0
        x_times[x_times >= self.i_max] = self.i_max - 1
        y_times[y_times >= self.j_max] = self.j_max - 1
        for p in range(0, x_times.shape[0]):
            xp = x_times[p, :]
            yp = y_times[p, :]
            id_p = ~np.isnan(xp)

            dom_vals[yp[id_p].astype(int), xp[id_p].astype(int)] += 1

        dom_vals[~np.isnan(dom_vals)] /= x_times.shape[0] * x_times.shape[1]
        print(str(x_times.shape[1]))
        # dom_vals[(dom_vals==0)] = np.nan
        fig, axs = plt.subplots(figsize=(12, 8), layout="constrained")
        cmap = plt.get_cmap("hot")
        cmap.set_bad("gray", 0.4)
        # vmax = np.nanmax(dom_vals)/3
        vmax = 0.01
        dom_map = axs.pcolormesh(dom_vals, cmap=cmap, alpha=1, vmax=vmax)
        axs.contour(
            self.depth,
            levels=self.bath_contours,
            colors="k",
            alpha=0.8,
            linewidths=1.5,
            zorder=2,
        )
        axs.set_title(str(self.dates[0]), fontsize=14)
        cbar = fig.colorbar(dom_map)
        cbar.ax.tick_params(labelsize=12)
        cbar.ax.set_ylabel("probability (%)", loc="center", size=12, weight="bold")
        self.save_plot(plt_name="dom_paths")
        return

    def plot_trajectory_color(self, skip_n, skip_t, kk):
        fig, axs = plt.subplots(figsize=(12, 8), layout="constrained")
        # cmap = plt.get_cmap('Blues')
        cmap = cmr.ocean_r
        cmap.set_bad("gray", 1.0)
        d_map = axs.pcolormesh(self.depth, cmap=cmap, alpha=0.75, vmax=4000)
        axs.contour(
            self.depth,
            levels=self.bath_contours,
            colors="k",
            alpha=0.35,
            linewidths=1.5,
        )
        x_times = self.df["xp"][:, kk, ::skip_t]
        y_times = self.df["yp"][:, kk, ::skip_t]
        # time_stamps = num2date(self.df['time'][::skip_t], self.df['time'].unit)
        t_vals = np.arange(0, skip_t * (4 / 24) * x_times.shape[1], skip_t * (4 / 24))

        t = np.arange(x_times.shape[1])
        for i in range(0, self.df["xp"].shape[0], skip_n):
            x = x_times[i]
            y = y_times[i]
            axs.plot(x, y, "k", alpha=0.8, zorder=1, linewidth=1.5)
            axs.scatter(x, y, c=t, s=2.3, cmap=cmr.ember, zorder=2)
        scatter_map = axs.scatter(x, y, c=t_vals, s=2, cmap=cmr.ember, zorder=2)
        cbar1 = fig.colorbar(d_map)
        cbar1.ax.tick_params(labelsize=12)
        cbar1.ax.set_ylabel("depth (m)", loc="center", size=12, weight="bold")
        cbar = fig.colorbar(scatter_map)
        cbar.ax.tick_params(labelsize=12)
        cbar.ax.set_ylabel("time elapsed (days)", loc="center", size=12, weight="bold")
        self.save_plot(plt_name="trajectory_color")
        return

    def plot_currents(self):
        fig, axs = plt.subplots(figsize=(12, 8), layout="constrained")
        mag_uv0 = np.sqrt(np.square(self.u_east) + np.square(self.v_north))
        cmap = plt.get_cmap("viridis")
        cmap.set_bad(color="gray", alpha=0.6)
        d_map = axs.pcolormesh(mag_uv0, cmap=cmap, vmax=0.5)
        axs.contour(
            self.depth,
            levels=self.bath_contours,
            colors="k",
            alpha=0.35,
            linewidths=1.3,
        )
        sk = 12
        xx = np.arange(0, self.i_max, 1)
        yy = np.arange(0, self.j_max, 1)
        [ux, yx] = np.meshgrid(xx, yy)
        axs.quiver(
            ux[::sk, ::sk],
            yx[::sk, ::sk],
            self.u_east[::sk, ::sk],
            self.v_north[::sk, ::sk],
            edgecolors="k",
            alpha=0.7,
            linewidths=0.1,
        )
        # axs.scatter(self.x, self.y, s=0.4, c='r')
        fig.colorbar(d_map)
        self.save_plot(plt_name="current_map")
        return

    def animate_dom(self, kk, skip_t):
        import matplotlib.animation as animation

        fig, ax = plt.subplots(figsize=(12, 8), layout="constrained")
        # fig, ax = plt.subplots()
        cmap = plt.get_cmap("hot")
        cmap.set_bad("gray", 0.9)

        artists = []
        x_times = np.array(self.df["xp"][:, kk, ::skip_t])
        y_times = np.array(self.df["yp"][:, kk, ::skip_t])
        x_times[x_times > 1000] = np.nan
        y_times[y_times > 1000] = np.nan

        x_times[x_times >= self.i_max] = self.i_max - 1
        y_times[y_times >= self.j_max] = self.j_max - 1
        second_point = int(x_times.shape[1] / 5)
        stept = int(second_point / 5)
        it_a = 0
        it_b = second_point
        for j in range(0, int(x_times.shape[1] / stept)):
            dom_vals = np.zeros(self.depth.shape)
            dom_vals[:] = self.depth[:]
            dom_vals[~np.isnan(dom_vals)] = 0
            it_a += stept
            it_b += stept
            if it_b >= x_times.shape[1]:
                break
            print(str(it_a))
            print("to: " + str(it_b))
            print(str(np.shape(x_times[0, slice(it_a, it_b, 1)])))
            # if j == int(x_times.shape[1]/stept)-1:
            # breakpoint()
            for p in range(0, x_times.shape[0]):
                xp = x_times[p, slice(it_a, it_b, 1)]
                yp = y_times[p, slice(it_a, it_b, 1)]
                id_p = ~np.isnan(xp)
                dom_vals[yp[id_p].astype(int), xp[id_p].astype(int)] += 1
            dom_vals[~np.isnan(dom_vals)] /= np.nanmax(dom_vals)

            # vmax=0.0035
            # plt.contourf(df.columns, df.index, df, 800, colors='white')
            if j == 0:
                ax.contour(
                    self.depth,
                    levels=self.bath_contours,
                    colors="k",
                    alpha=0.35,
                    linewidths=1.5,
                )
                plt.title("datetime = " + str(self.dates[0]), fontsize=25)
            container = ax.contourf(
                dom_vals, cmap=cmap, alpha=1, levels=np.linspace(0, 0.03, 40)
            )

            print(str(self.dates[j * int(x_times.shape[1] / stept)]))

            # if j == 0:
            # vmax = np.nanmax(dom_vals) / 2
            # fig.colorbar(container)

            artists.append(container.collections)
        ani = animation.ArtistAnimation(
            fig=fig, artists=artists, interval=400, blit=True
        )
        ani.save(filename=self.save_folder + "ani_examp.gif", writer="Pillow")

        #

    def init_current_file(self):
        print("creating current file")
        samples_folder = "E:/fromIngrid/samplefiles_reruns/"
        samples_prefix = "samplesNSEW_"
        init_year = 2017
        init_month = 2
        filename = (
            samples_folder
            + samples_prefix
            + str(init_year)
            + "{:02d}".format(init_month)
            + ".nc"
        )
        df = nc.Dataset(filename)
        self.u_east = np.array(np.nanmean(df["u_east"][0:10, 0, :, :], 0))
        self.v_north = np.array(np.nanmean(df["v_north"][0:10, 0, :, :], 0))
        # self.u_east[self.u_east<-1000] = np.nan
        # self.v_north[self.v_north < -1000] = np.nan
        np.save(self.u_file, self.u_east)
        np.save(self.v_file, self.v_north)
        return

    def init_depth_file(self):
        print("creating depth file")
        samples_folder = "E:/fromIngrid/samplefiles_reruns/"
        samples_prefix = "samplesNSEW_"
        init_year = 2017
        init_month = 2
        filename = (
            samples_folder
            + samples_prefix
            + str(init_year)
            + "{:02d}".format(init_month)
            + ".nc"
        )
        df = nc.Dataset(filename)
        self.depth = np.array(df["depth"])
        self.depth[self.depth > 10000] = np.nan
        np.save(self.depth_file, self.depth)
        return

    def save_plot(self, plt_name):
        savefile = self.save_folder + plt_name + ".png"
        print("Saving file: " + savefile)
        plt.savefig(savefile, dpi=400)
        plt.close()
        return


def plot_physics(filename):
    x_min = 500
    x_max = 630
    x_slice = slice(x_min, x_max, 1)
    y_min = 475
    y_max = 625
    y_slice = slice(y_min, y_max, 1)

    df = nc.Dataset(filename)
    temp_mat = df["temperature"]

    # temp = temp_mat[0, 0, :, :]
    # temp[y_slice, x_slice] = np.nan
    # plt.pcolormesh(temp)
    # plt.show()

    time_list = num2date(df["time"], df["time"].units)
    depth_v = np.cumsum(df["LayerDepths"])

    fig, axs = plt.subplots()

    it_num = range(0, temp_mat.shape[0], 75)
    for i in it_num:
        print(str(i))
        temp_depth = temp_mat[i, :, :, :]
        # t_mean1 = np.nanmean(temp_depth[:, y_slice,  x_slice], 2)
        t_mean1 = np.nanmean(temp_depth[:, :, :], 2)
        t_mean2 = np.nanmean(t_mean1, 1) - 273.15
        axs.plot(t_mean2, depth_v, "k-")
        axs.set_ylabel("depth (m)")
        axs.set_xlabel("temperature (degC)")
    plt.title("month = " + str(time_list[0].month))
    plt.gca().invert_yaxis()
    plt.show()
    breakpoint()
