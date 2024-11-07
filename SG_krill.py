from netCDF4 import date2num
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt


class Krill:
    def __init__(self, nc_file):
        self.i_max = nc_file["xc"].shape[0]
        self.j_max = nc_file["yc"].shape[0]
        self.lay_depths = np.cumsum(np.array(nc_file["LayerDepths"][:]))
        self.depth = np.array(nc_file["depth"])
        self.depth[self.depth > 10000] = np.nan
        self.res = 800
        self.load_counter = -1
        self.bio_load_counter = -1
        self.init_temp = False
        self.kpar = 0.2  # mg m^-3
        self.gpar = 0.044  # mg mm^-2 day^-1 (assimilate ad libitum)
        self.lpar = 40  # mm
        self.lmax = 60  # mm
        self.r_ref = 0.0035  # day -1
        self.T_A = 8000  # K
        degC = 0  # C
        self.T_1 = degC + 273.15
        return

    def init_temp_response(self, readerSG):
        sqN_int = np.sqrt(readerSG.N).astype(int)
        t_scale = np.linspace(-1.5, 1.5, sqN_int)
        w_scale = np.linspace(1, np.sqrt(readerSG.N), sqN_int)
        [t_vals, w_vals] = np.meshgrid(t_scale, w_scale)
        self.t_min = np.reshape(t_vals, [readerSG.N])
        self.t_max = self.t_min + 3
        self.w_max = np.reshape(w_vals, [readerSG.N]) * 0.06
        readerSG.logger.info("min temps: " + str(self.t_min))
        readerSG.logger.info("max temps: " + str(self.t_max))
        readerSG.logger.info("max w: " + str(self.w_max))
        self.init_temp = True
        return

    def init_krill(self, reader_SG):
        self.x = np.zeros([reader_SG.n, reader_SG.N])
        self.y = np.zeros([reader_SG.n, reader_SG.N])
        self.z = np.zeros([reader_SG.n, reader_SG.N])
        self.l = np.zeros([reader_SG.n, reader_SG.N])
        length_list = [40, 40]
        depth_list = [50, 250]
        step_xy = reader_SG.n / np.sqrt(reader_SG.n)

        for kk in range(0, reader_SG.N, 1):
            kz = depth_list[kk]
            ll = length_list[kk]
            reader_SG.logger.info(
                "Initialised ensemble (#) "
                + str(kk)
                + " at depth (m) = "
                + str(kz)
                + " with length (mm) = "
                + str(ll)
            )
            counter = -1
            for ii in np.arange(0, step_xy, 1):
                x = reader_SG.x_min + (
                    (reader_SG.x_max - reader_SG.x_min) * ((ii + 1) / step_xy)
                )
                for jj in np.arange(0, step_xy, 1):
                    y = reader_SG.y_min + (
                        (reader_SG.y_max - reader_SG.y_min) * ((jj + 1) / step_xy)
                    )
                    counter = counter + 1
                    if self.depth[np.floor(y).astype(int), np.floor(x).astype(int)] > 0:
                        self.x[counter, kk] = x
                        self.y[counter, kk] = y
                        self.z[counter, kk] = kz
                        self.l[counter, kk] = ll
                    else:
                        self.x[counter, kk] = np.nan
                        self.y[counter, kk] = np.nan
                        self.z[counter, kk] = np.nan
                        self.l[counter, kk] = np.nan
        reader_SG.logger.info(
            "Initialised "
            + str(np.nansum(self.x > 0))
            + " krill with x and y values across "
            + str(reader_SG.N)
            + " ensemble members"
        )
        return

    def step_krill(self, reader_SG):
        time_index = reader_SG.time_index
        biotime_index = reader_SG.bio_index
        nc_file = reader_SG.nc_file
        dt = reader_SG.dt.seconds
        if biotime_index != self.bio_load_counter:
            self.chl = np.array(reader_SG.bio_ncfile["chl"][biotime_index, :, :, :])
            self.chl[self.chl > 10000] = np.nan
            self.bio_load_counter = biotime_index
            reader_SG.logger.info(
                "loading new biology variables for datetime: "
                + str(reader_SG.current_datetime)
            )
        if time_index != self.load_counter:
            self.u_east = np.array(nc_file["u_east"][time_index, :, :, :])
            self.v_north = np.array(nc_file["v_north"][time_index, :, :, :])
            self.temp = np.array(nc_file["temperature"][time_index, :, :, :])
            self.temp[self.temp < -3000] = np.nan
            self.v_north[self.v_north < -3000] = np.nan
            self.u_east[self.u_east < -3000] = np.nan
            self.load_counter = time_index
            reader_SG.logger.info(
                "loading new physics variables for datetime: "
                + str(reader_SG.current_datetime)
            )
        self.t_save = np.ones([self.x.shape[0], self.x.shape[1]]) * -32000
        self.w_save = np.ones([self.x.shape[0], self.x.shape[1]]) * -32000
        # self.g_save = np.ones([self.x.shape[0], self.x.shape[1]]) * -32000
        for kk in range(0, self.x.shape[1]):
            for ii in range(0, self.x.shape[0]):
                if (
                    (self.x[ii, kk] > 0)
                    & (self.x[ii, kk] < self.i_max)
                    & (self.y[ii, kk] < self.j_max)
                    & (self.y[ii, kk] > 0)
                ):
                    zv = self.z[ii, kk]
                    layer_ii = np.argmin((self.lay_depths - zv) ** 2)
                    u = self.u_east[
                        layer_ii,
                        np.floor(self.y[ii, kk]).astype(int),
                        np.floor(self.x[ii, kk]).astype(int),
                    ]
                    v = self.v_north[
                        layer_ii,
                        np.floor(self.y[ii, kk]).astype(int),
                        np.floor(self.x[ii, kk]).astype(int),
                    ]
                    t = (
                        self.temp[
                            layer_ii,
                            np.floor(self.y[ii, kk]).astype(int),
                            np.floor(self.x[ii, kk]).astype(int),
                        ]
                        - 273.15
                    )  # in celsius

                    # bio_index
                    if (np.isnan(t)) | (np.isnan(u)) | (np.isnan(v)):
                        self.x[ii, kk] = np.nan
                        self.y[ii, kk] = np.nan
                        self.z[ii, kk] = np.nan
                    elif (
                        (self.x[ii, kk] > self.i_max)
                        | (self.y[ii, kk] > self.j_max)
                        | (self.y[ii, kk] < 1)
                        | (self.x[ii, kk] < 1)
                    ):
                        self.x[ii, kk] = np.nan
                        self.y[ii, kk] = np.nan
                        self.z[ii, kk] = np.nan
                    elif np.isnan(
                        self.depth[
                            np.floor(self.y[ii, kk]).astype(int),
                            np.floor(self.x[ii, kk]).astype(int),
                        ]
                    ):
                        self.x[ii, kk] = np.nan
                        self.y[ii, kk] = np.nan
                        self.z[ii, kk] = np.nan
                    else:
                        x1 = int(self.x[ii, kk])
                        y1 = int(self.y[ii, kk])
                        z1 = int(zv)
                        l1 = self.l[ii, kk]

                        if reader_SG.feed_beh:
                            kgrowth = self.feeding(reader_SG, x1, y1, l1, layer_ii, t)
                            if np.isnan(kgrowth):
                                print("invalid growth value;")
                                self.l[ii, kk] += 0
                            else:
                                self.l[ii, kk] += (kgrowth / 86400) * dt

                        if reader_SG.dvm_beh:
                            w = 0.0
                        else:
                            w = 0.0
                            # self.dvm()

                        self.x[ii, kk] = self.x[ii, kk] + ((dt * u) / self.res)
                        self.y[ii, kk] = self.y[ii, kk] + ((dt * v) / self.res)
                        self.z[ii, kk] = self.z[ii, kk] + ((dt * w) / self.res)
                        self.t_save[ii, kk] = t
                        self.w_save[ii, kk] = w
                        # self.ing_save[ii, kk] = ing
                else:
                    self.x[ii, kk] = np.nan
                    self.y[ii, kk] = np.nan
                    self.z[ii, kk] = np.nan
        return

    def dvm(self, reader_SG, x1, y1, z1):
        lonid = reader_SG.light_lon_id[y1, x1]
        latid = reader_SG.light_lat_id[y1, x1]

    def feeding(self, reader_SG, x1, y1, l1, layer_ii, t):
        lonid = reader_SG.lon_id[y1, x1]
        latid = reader_SG.lat_id[y1, x1]
        dep_id = reader_SG.dep_id[layer_ii]
        d_min = np.nanmax([0, dep_id - 2])
        d_max = np.nanmin([self.chl.shape[0] - 2, dep_id + 2])
        d_slice = slice(int(d_min), int(d_max), 2)
        latid_min = np.nanmax([0, latid - 2])
        latid_max = np.nanmin([self.chl.shape[1], latid + 2])
        lat_slice = slice(int(latid_min), int(latid_max), 1)
        lonid_min = np.nanmax([0, lonid - 2])
        lonid_max = np.nanmin([self.chl.shape[2], lonid + 2])
        lon_slice = slice(int(lonid_min), int(lonid_max), 1)
        chl_val = np.nanmean(self.chl[d_slice, lat_slice, lon_slice])
        # print(str(chl_val))
        # get rates:
        if chl_val > 0:
            f = chl_val / (chl_val + self.kpar)
        else:
            f = 0
        # ing = self.gpar * f * self.lpar
        T = t + 273.15
        kgrowth = (
            (self.lmax - l1)
            * self.r_ref
            * np.exp((self.T_A / self.T_1) - (self.T_A / T))
            * f
        )
        # print('chl (mg m-3) = ' + str(chl_val))
        # print('f = ' + str(f))
        # print('T = ' + str(t))
        # print('ing (mg day-1) = ' + str(ing))
        # print('growth (mm day-1) = ' + str(kgrowth))
        return kgrowth

    def swim_temp(self, t, kk, grad_t):
        if (grad_t == 0) | (np.isnan(grad_t)):
            scale_w = np.random.choice([-1, 1])
            w = scale_w * self.w_max[kk]
        elif t < self.t_min[kk]:
            # scale_w = ((t - self.t_min) ** 2) / (2 ** 2)
            w = np.sign(grad_t) * self.w_max[kk]
        else:
            # scale_w = (((t - self.t_max) ** 2) / (2 ** 2))
            w = -1 * np.sign(grad_t) * self.w_max[kk]
        return w

    def save_step(self, save_counter, reader_SG):
        self.trajectory_file["xp"][:, :, save_counter] = self.x
        self.trajectory_file["yp"][:, :, save_counter] = self.y
        self.trajectory_file["zp"][:, :, save_counter] = self.z
        self.trajectory_file["w"][:, :, save_counter] = self.w_save
        self.trajectory_file["temp"][:, :, save_counter] = self.t_save
        if reader_SG.feed_beh:
            self.trajectory_file["lp"][:, :, save_counter] = self.l
        self.trajectory_file["time"][save_counter] = date2num(
            reader_SG.current_datetime, self.trajectory_file["time"].unit
        )
        return

    def init_netcdf(self, readerSG):
        trajectory_filename = readerSG.trajectory_folder + readerSG.save_file

        self.trajectory_file = nc.Dataset(trajectory_filename, mode="w")
        dimension_key_dict = {
            "trajectory": readerSG.n,
            "ensemble": readerSG.N,
            "obs": readerSG.save_number,
        }

        for dimension in dimension_key_dict:
            self.trajectory_file.createDimension(
                dimension, dimension_key_dict[dimension]
            )

        variable_key_dict = {
            "xp": {
                "datatype": "f4",
                "dimensions": ("trajectory", "ensemble", "obs"),
                "description": "x position of particle",
            },
            "yp": {
                "datatype": "f4",
                "dimensions": ("trajectory", "ensemble", "obs"),
                "description": "y position of particle",
            },
            "zp": {
                "datatype": "f4",
                "dimensions": ("trajectory", "ensemble", "obs"),
                "description": "z position of particle",
            },
            "lp": {
                "datatype": "f4",
                "dimensions": ("trajectory", "ensemble", "obs"),
                "description": "length based on instantaneous growth (mm day**-1) coupled T and Chl_a",
            },
            "temp": {
                "datatype": "f4",
                "dimensions": ("trajectory", "ensemble", "obs"),
                "description": "temperature of particle",
            },
            "w": {
                "datatype": "f4",
                "dimensions": ("trajectory", "ensemble", "obs"),
                "description": "speed of particle",
            },
            "time": {
                "datatype": "f4",
                "dimensions": ("obs",),
                "description": "datetime of particle",
            },
        }

        if readerSG.temp_beh:
            variable_key_dict.update(
                {
                    "t_min": {
                        "datatype": "f4",
                        "dimensions": ("ensemble",),
                        "description": "minimum temperature in comfort zone",
                    },
                    "t_max": {
                        "datatype": "f4",
                        "dimensions": ("ensemble",),
                        "description": "max temperature in comfort zone",
                    },
                    "w_max": {
                        "datatype": "f4",
                        "dimensions": ("ensemble",),
                        "description": "speed of discomfort response",
                    },
                }
            )

        # readerSG.logger.info('Dictionary structure for saving: ' + str(variable_key_dict))

        for variable in variable_key_dict:
            self.trajectory_file.createVariable(
                variable,
                variable_key_dict[variable]["datatype"],
                variable_key_dict[variable]["dimensions"],
            )
            self.trajectory_file[variable].description = variable_key_dict[variable][
                "description"
            ]

        if self.init_temp:
            self.trajectory_file["t_min"][:] = self.t_min
            self.trajectory_file["t_max"][:] = self.t_max
            self.trajectory_file["w_max"][:] = self.w_max

        # initialise time units
        time_unit_out = "seconds since 2014-04-01 00:00:00"
        self.trajectory_file["time"].setncattr("unit", time_unit_out)
        readerSG.logger.info("Initialising trajectory file: " + trajectory_filename)
        readerSG.logger.info(self.trajectory_file)
        return

    def plot_init(self, kk):
        fig, axs = plt.subplots(figsize=(12, 8), layout="constrained")
        cmap = plt.get_cmap("Blues")
        cmap.set_bad("gray", 0.8)
        d_map = axs.pcolormesh(self.depth, cmap=cmap, vmax=4500)
        fig.colorbar(d_map, label="depth (m)")
        axs.scatter(self.x[:, kk], self.y[:, kk], c="r", s=3)
        plt.show()
        return

    def plot_currents(self):
        fig, axs = plt.subplots(figsize=(12, 8), layout="constrained")
        mag_uv0 = np.sqrt(np.square(self.u_east) + np.square(self.v_north))
        d_map = axs.pcolormesh(mag_uv0, cmap=plt.get_cmap("viridis"), vmin=0)
        sk = 10
        xx = np.arange(0, self.i_max, 1)
        yy = np.arange(0, self.j_max, 1)
        [ux, yx] = np.meshgrid(xx, yy)
        axs.quiver(
            ux[::sk, ::sk],
            yx[::sk, ::sk],
            self.u_east[::sk, ::sk],
            self.v_north[::sk, ::sk],
            edgecolors="k",
            alpha=0.5,
            linewidths=0.05,
        )
        axs.scatter(self.x, self.y, s=0.4, c="r")
        fig.colorbar(d_map)
        plt.show()
        return
