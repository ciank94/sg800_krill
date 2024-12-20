import math
import matplotlib.pyplot as plt
import numpy as np
import pyproj
import netCDF4 as nc
import os
import datetime
from netCDF4 import num2date, date2num
from yaml import Loader, load as yml_load
import logging


class SGReader:
    def __init__(self, namelist_path, release_number):
        # initiate logger
        self.logger = logging.getLogger(__name__)
        self.logger.warning("========================")
        self.logger.warning("Beginning new simulation ")

        # load namelist:
        namelist_file = namelist_path + "namelist.yml"
        stream = open(namelist_file, "r")
        namelist = yml_load(stream, Loader)

        # switch list:
        switch_list = namelist["switches"]
        self.remote = switch_list["remote"]
        self.test = switch_list["test"]
        self.feed_beh = switch_list["feed_beh"]
        self.dvm_beh = switch_list["dvm_beh"]
        self.temp_beh = switch_list["temp_beh"]
        self.light_mapping = switch_list["light_mapping"]
        self.bio_mapping = switch_list["bio_mapping"]
        if (self.dvm_beh) & (self.temp_beh):
            self.logger.warning("dvm_beh overriding temp_beh")

        # file explorer:
        file_list = namelist["file_explorer"]
        self.samples_prefix = file_list["samples_prefix"]
        self.exp_tag = file_list["experiment_tag"]
        fsamples_prefix = []
        ftrajectory_prefix = []
        if self.remote:
            if file_list["remote_server"] == "idun":
                fsamples_prefix = file_list["remote"]["idun_dir_prefix"]
                ftrajectory_prefix = fsamples_prefix
            elif file_list["remote"]["server"] == "saga":
                fsamples_prefix = file_list["remote"]["saga_dir_prefix"]
                ftrajectory_prefix = fsamples_prefix
            else:
                os.error("not a valid server name")
        else:
            ftrajectory_prefix = file_list["local"]["local_dir_prefix"]
            if file_list["remote_server"] == "idun":
                fsamples_prefix = file_list["local"]["idun_dir_prefix"]
            elif file_list["remote_server"] == "saga":
                fsamples_prefix = file_list["local"]["saga_dir_prefix"]
            else:
                os.error("not a valid server name")

        samples_folder = fsamples_prefix + file_list["samples_folder_name"]
        trajectory_folder = (
            ftrajectory_prefix
            + file_list["sg_folder_name"]
            + file_list["trajectory_folder_name"]
        )
        self.samples_folder = samples_folder
        self.trajectory_folder = trajectory_folder

        self.logger.warning("samples folder = " + self.samples_folder)
        self.logger.warning("trajectory folder = " + self.trajectory_folder)
        self.logger.warning("========================")

        # time settings
        time_list = namelist["time_settings"]
        year1 = time_list["year"]
        month1 = time_list["month"]
        day1 = time_list["day"]
        day_inc = datetime.timedelta(days=release_number)  # IMPORTANT
        self.init_datetime = datetime.datetime(year1, month1, day1, 1, 0) + day_inc
        days1 = time_list["duration_days"]
        duration_days = datetime.timedelta(days=days1)  # simulation duration in days;
        self.duration = duration_days.days
        self.current_datetime = self.init_datetime
        self.file_month = self.init_datetime.month
        # load SG800 physics:
        self.filename = (
            self.samples_folder
            + self.samples_prefix
            + str(self.init_datetime.year)
            + "{:02d}".format(self.init_datetime.month)
            + ".nc"
        )
        self.nc_file = nc.Dataset(self.filename)
        self.imax = self.nc_file["xc"].shape[0]
        self.jmax = self.nc_file["yc"].shape[0]
        self.time = num2date(self.nc_file["time"], self.nc_file["time"].units)
        minutes = time_list["time_step_minutes"]
        self.dt = datetime.timedelta(hours=minutes / 60)
        save_step_hours = time_list["save_step_hours"]
        self.save_step = datetime.timedelta(hours=save_step_hours)  # save time step
        self.time_threshold = datetime.timedelta(seconds=0)
        self.simulation_steps = duration_days / self.dt
        self.save_number = int(duration_days / self.save_step)

        # ibm settings:
        ibm_list = namelist["ibm_settings"]
        self.n = ibm_list["n"]
        self.N = ibm_list["N"]
        self.x_min = ibm_list["x_min"]
        self.x_max = ibm_list["x_max"]
        self.y_min = ibm_list["y_min"]
        self.y_max = ibm_list["y_max"]

        # saving and logging:
        self.save_file_prefix = (
            "trajectory_"
            + str(self.init_datetime.year)
            + "{:02d}".format(self.init_datetime.month)
            + "{:02d}".format(self.init_datetime.day)
            + "_d"
            + str(self.duration)
            + "_"
            + self.exp_tag
        )
        self.save_file = self.save_file_prefix + ".nc"
        current_datenum = date2num(self.current_datetime, self.nc_file["time"].units)
        self.time_index = np.argmin((current_datenum - self.nc_file["time"][:]) ** 2)

        # set up logging to file
        log_filename = (
            "logger_"
            + str(self.init_datetime.year)
            + "{:02d}".format(self.init_datetime.month)
            + "{:02d}".format(self.init_datetime.day)
            + "_d"
            + str(self.duration)
            + "_"
            + self.exp_tag
        )
        logging.basicConfig(
            filename=self.trajectory_folder + log_filename + ".log",
            level=logging.INFO,
            filemode="w",
            format="[%(asctime)s] - %(message)s",
            datefmt="%H:%M:%S",
        )
        self.logger.info("========================")
        self.logger.info("=======filepaths========")
        self.logger.info("opening new SG file: " + self.filename)
        self.logger.warning("SG file: " + self.filename)
        self.logger.warning("trajectory file: " + self.save_file)

        # load biology file:
        self.read_bio()

        # load light file:
        self.read_light()

        # load initial values:
        self.log_init()
        return

    def read_light(self):
        self.light_prefix = "swrad_"
        self.light_filename = (
            self.samples_folder
            + self.light_prefix
            + str(self.current_datetime.year)
            + ".nc"
        )
        self.logger.info("opening light file " + self.light_filename)
        self.light_ncfile = nc.Dataset(self.light_filename)
        self.light_index = np.argmin(
            (
                date2num(self.current_datetime, self.light_ncfile["valid_time"].units)
                - self.light_ncfile["valid_time"][:]
            )
            ** 2
        )

        #
        # num2date(self.light_ncfile['valid_time'][self.light_index], self.light_ncfile['valid_time'].units)
        self.light_lon_mapping = (
            self.samples_folder
            + "lon_mapping_"
            + str(self.current_datetime.year)
            + ".npy"
        )
        self.light_lat_mapping = (
            self.samples_folder
            + "lat_mapping_"
            + str(self.current_datetime.year)
            + ".npy"
        )

        if self.light_mapping:
            self.logger.info("creating mapping from " + self.light_filename)
            lat2 = np.array(self.light_ncfile["latitude"][:])
            lon2 = np.array(self.light_ncfile["longitude"][:])
            [x_cmems, y_cmems] = np.meshgrid(lat2, lon2)
            x = np.arange(0, self.nc_file["xc"].shape[0])
            y = np.arange(0, self.nc_file["yc"].shape[0])
            [x_pos, y_pos] = np.meshgrid(x, y)
            [lat1, lon1] = geo2grid(x_pos.flatten(), y_pos.flatten(), "get_bl")
            lat_re = lat1.reshape([y.shape[0], x.shape[0]])
            lon_re = lon1.reshape([y.shape[0], x.shape[0]])
            lat_grid = np.zeros(lat_re.shape)
            lon_grid = np.zeros(lon_re.shape)
            for i in range(0, lat_re.shape[0]):
                for j in range(0, lat_re.shape[1]):
                    dist_v = haversine(y_cmems, x_cmems, lon_re[i, j], lat_re[i, j])
                    idx = np.argmin(dist_v)
                    [id_lon, id_lat] = np.unravel_index(idx, dist_v.shape)
                    if (id_lon < x_cmems.shape[0]) & (id_lat < x_cmems.shape[1]):
                        lon_grid[i, j] = (
                            id_lon  # value for where we should look in bio_file chl(id_lon, id_lat)
                        )
                        lat_grid[i, j] = id_lat

            np.save(self.light_lon_mapping, lon_grid)
            np.save(self.light_lat_mapping, lat_grid)
            self.logger.info("saved mapping from " + self.light_filename)
        self.logger.info("loading mapping from " + self.light_filename)
        self.light_lon_id = np.load(self.light_lon_mapping)
        self.light_lat_id = np.load(self.light_lat_mapping)

        # check section:
        # import matplotlib.pyplot as plt
        # import cartopy.crs as ccrs
        # light_array = np.array(self.light_ncfile['msnswrf'][self.light_index+20,:,:])
        # dp_vals = self.nc_file['depth'][:]
        # dp_vals[dp_vals>0] = 0
        # check_vals = np.zeros(self.light_lon_id.shape)
        # for ii in range(0, self.light_lon_id.shape[0]):
        #     for jj in range(0, self.light_lon_id.shape[1]):
        #         check_vals[ii,jj] += light_array[int(self.light_lat_id[ii,jj]), int(self.light_lon_id[ii,jj])]
        #         check_vals[ii,jj] += dp_vals[ii,jj]
        # dmap = plt.pcolormesh(check_vals, cmap=plt.get_cmap('hot'))
        # plt.colorbar(label='Wm**-2')
        # plt.show()
        #
        # lon = self.light_ncfile['longitude'][:]
        # lat = self.light_ncfile['latitude'][:]
        # fig, axs = plt.subplots(figsize=(8, 8), nrows=1, ncols=1, subplot_kw={'projection': ccrs.PlateCarree()})
        # fig1 = axs.pcolormesh(lon, lat, light_array, cmap=plt.get_cmap('hot'),
        #                       transform=ccrs.PlateCarree())
        # axs.coastlines()
        # fig.colorbar(fig1, label='short-wave rad' + ' units =' + 'Wm**-2')
        # axs.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
        # plt.show()
        # breakpoint()

        # check the overall light values;
        # light = self.light_ncfile['msnswrf']
        # l_mean = np.mean(light, 1)
        # l_mean2 = np.mean(l_mean, 1)
        # plt.plot(np.arange(0, l_mean2.shape[0])/24, l_mean2, 'r')
        # plt.xlabel('day of year')
        # plt.ylabel('W m **-2')
        # plt.show()
        # breakpoint()
        # print(str(np.mean(light[:])))
        # [print(str(np.mean(light[i, :, :]))) for i in range(0, light.shape[0])]
        # breakpoint()
        return

    def read_bio(self):
        self.bio_fileprefix = "CMEMS_SGBIO_Dfull_"
        self.bio_filename = (
            self.samples_folder
            + self.bio_fileprefix
            + str(self.current_datetime.year)
            + ".nc"
        )
        self.logger.info("opening biology file " + self.bio_filename)
        self.bio_ncfile = nc.Dataset(self.bio_filename)
        self.bio_index = np.argmin(
            (
                date2num(self.current_datetime, self.bio_ncfile["time"].units)
                - self.bio_ncfile["time"][:]
            )
            ** 2
        )

        # mapping between PISCES grid and SINMOD
        self.lon_mapping = (
            self.samples_folder
            + "lon_mapping_"
            + str(self.current_datetime.year)
            + ".npy"
        )
        self.lat_mapping = (
            self.samples_folder
            + "lat_mapping_"
            + str(self.current_datetime.year)
            + ".npy"
        )
        self.depth_mapping = (
            self.samples_folder
            + "depth_mapping_"
            + str(self.current_datetime.year)
            + ".npy"
        )

        if self.bio_mapping:
            self.logger.info("creating mapping from " + self.bio_filename)
            lat2 = np.array(self.bio_ncfile["latitude"][:])
            lon2 = np.array(self.bio_ncfile["longitude"][:])
            [x_cmems, y_cmems] = np.meshgrid(lat2, lon2)
            x = np.arange(0, self.nc_file["xc"].shape[0])
            y = np.arange(0, self.nc_file["yc"].shape[0])
            [x_pos, y_pos] = np.meshgrid(x, y)
            [lat1, lon1] = geo2grid(x_pos.flatten(), y_pos.flatten(), "get_bl")
            lat_re = lat1.reshape([y.shape[0], x.shape[0]])
            lon_re = lon1.reshape([y.shape[0], x.shape[0]])
            lat_grid = np.zeros(lat_re.shape)
            lon_grid = np.zeros(lon_re.shape)
            for i in range(0, lat_re.shape[0]):
                for j in range(0, lat_re.shape[1]):
                    dist_v = haversine(y_cmems, x_cmems, lon_re[i, j], lat_re[i, j])
                    idx = np.argmin(dist_v)
                    [id_lon, id_lat] = np.unravel_index(idx, dist_v.shape)
                    if (id_lon < x_cmems.shape[0]) & (id_lat < x_cmems.shape[1]):
                        lon_grid[i, j] = (
                            id_lon  # value for where we should look in bio_file chl(id_lon, id_lat)
                        )
                        lat_grid[i, j] = id_lat

            depth_bio = self.bio_ncfile["depth"][:]
            depth_layers = np.cumsum(self.nc_file["LayerDepths"][:])
            depth_grid = np.zeros(depth_layers.shape[0])
            for ii in range(0, depth_layers.shape[0]):
                idx_layer = np.argmin((depth_layers[ii] - depth_bio) ** 2)
                depth_grid[ii] = idx_layer

            np.save(self.depth_mapping, depth_grid)
            np.save(self.lon_mapping, lon_grid)
            np.save(self.lat_mapping, lat_grid)
            self.logger.info("saved mapping from " + self.bio_filename)
        self.logger.info("loading mapping from " + self.bio_filename)
        self.lon_id = np.load(self.lon_mapping)
        self.lat_id = np.load(self.lat_mapping)
        self.dep_id = np.load(self.depth_mapping)

        # check section:
        # import matplotlib.pyplot as plt
        # import cartopy.crs as ccrs
        # chl_array = np.array(self.bio_ncfile['chl'][self.bio_index+20,0, :,:])
        # chl_array[chl_array>3000] = np.nan
        # dp_vals = self.nc_file['depth'][:]
        # dp_vals[dp_vals>0] = 0
        # check_vals = np.zeros(self.lon_id.shape)
        # for ii in range(0, self.lon_id.shape[0]):
        #     for jj in range(0, self.lon_id.shape[1]):
        #         check_vals[ii,jj] += chl_array[int(self.lat_id[ii,jj]), int(self.lon_id[ii,jj])]
        #         check_vals[ii,jj] += dp_vals[ii,jj]
        # cmap = plt.get_cmap('Greens')
        # cmap.set_bad('gray', 0.8)
        # plt.pcolormesh(check_vals, cmap=cmap)
        # plt.colorbar(label='mg m**-3')
        # plt.show()

        # chlorophyll on original grid;
        # lon = self.bio_ncfile['longitude'][:]
        # lat = self.bio_ncfile['latitude'][:]
        # fig, axs = plt.subplots(figsize=(8, 8), nrows=1, ncols=1, subplot_kw={'projection': ccrs.PlateCarree()})
        # fig1 = axs.pcolormesh(lon, lat, chl_array, cmap=plt.get_cmap('Greens'),
        #                       transform=ccrs.PlateCarree())
        # axs.coastlines()
        # fig.colorbar(fig1, label='chl vals' + ' units =' + 'mg m**-3')
        # axs.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
        # plt.show()
        # breakpoint()

        # import matplotlib.pyplot as plt
        # fig, axs = plt.subplots()
        # T_A = 8000  # K
        # degC = 0  # C
        # T_1 = degC + 273.15
        # f=1
        # length_list = np.arange(1, 60, 0.1)
        # temp_list = [-1, 0, 1, 2, 3, 4, 5]
        # for t in temp_list:
        #    T = t + 273.15
        #    kgrowth = np.zeros(length_list.shape[0])
        #    counter = -1
        #    for l_v in length_list:
        #        counter += 1
        #        kgrowth[counter] = (60 - l_v) * 0.0035 * np.exp((T_A / T_1) - (T_A / T)) * f
        #    plt.plot(length_list, kgrowth)
        #    plt.legend(['-1', '0', '1', '2', '3', '4', '5'])
        # plt.xlabel('length (mm)')
        # plt.ylabel('instantaneous growth (mm day**-1)')
        # plt.show()

        # import matplotlib.pyplot as plt
        # Kval = 0.2
        # chl_vals = np.arange(0, 0.4, 0.01)
        # f = np.zeros(chl_vals.shape[0])
        # counter = -1
        # for cv in chl_vals:
        #    counter += 1
        #    f[counter] = cv/(cv + Kval)

        # plt.plot(chl_vals, f, 'k')
        # plt.xlabel('chl a (mg m**-3)')
        # plt.ylabel('f')
        # plt.show()
        return

    def log_init(self):
        self.logger.info("samples_filename = " + str(self.filename))
        self.logger.info("date_init = " + str(self.init_datetime))
        self.logger.info("==============================")
        self.logger.info("===========switches===========")
        self.logger.info("test = " + str(self.test))
        self.logger.info("dvm_beh = " + str(self.dvm_beh))
        self.logger.info("feed_beh = " + str(self.dvm_beh))
        self.logger.info("temp_beh = " + str(self.temp_beh))
        self.logger.info("=========ibm_settings=========")
        self.logger.info("n = " + str(self.n))
        self.logger.info("N = " + str(self.N))
        self.logger.info("x_min = " + str(self.x_min))
        self.logger.info("x_max = " + str(self.x_max))
        self.logger.info("y_min = " + str(self.y_min))
        self.logger.info("y_max = " + str(self.y_max))
        self.logger.info("imax =  " + str(self.imax))
        self.logger.info("jmax =  " + str(self.jmax))
        self.logger.info("=========time_settings========")
        self.logger.info("duration (days) = " + str(self.duration))
        self.logger.info("time step (hours) = " + str(self.dt))
        self.logger.info("save step (hours) = " + str(self.save_step))
        self.logger.info("simulation steps(#) = " + str(int(self.simulation_steps)))
        self.logger.info("save number(#) = " + str(int(self.save_number)))
        self.logger.info("==============================")
        return

    def update_time(self):
        self.current_datetime += self.dt
        if self.current_datetime.month != self.file_month:
            self.nc_file.close()
            self.update_file_month(
                samples_folder=self.samples_folder, samples_prefix=self.samples_prefix
            )
            new_file = True
        else:
            new_file = False

        bio_date = num2date(
            self.bio_ncfile["time"][self.bio_index], self.bio_ncfile["time"].units
        )

        if (self.current_datetime.day != bio_date.day) | new_file:
            self.bio_index = np.argmin(
                (
                    date2num(self.current_datetime, self.bio_ncfile["time"].units)
                    - self.bio_ncfile["time"][:]
                )
                ** 2
            )
            # bio_date = num2date(self.bio_ncfile['time'][self.bio_index], self.bio_ncfile['time'].units)

        if self.test:
            if (self.current_datetime.day != self.time[self.time_index].day) | new_file:
                current_datenum = date2num(
                    self.current_datetime, self.nc_file["time"].units
                )
                self.time_index = np.argmin(
                    (current_datenum - self.nc_file["time"][:]) ** 2
                )
        else:
            if (
                self.current_datetime.hour != self.time[self.time_index].hour
            ) | new_file:
                current_datenum = date2num(
                    self.current_datetime, self.nc_file["time"].units
                )
                self.time_index = np.argmin(
                    (current_datenum - self.nc_file["time"][:]) ** 2
                )
        return

    def update_file_month(self, samples_folder, samples_prefix):
        self.filename = (
            samples_folder
            + samples_prefix
            + str(self.current_datetime.year)
            + "{:02d}".format(self.current_datetime.month)
            + ".nc"
        )
        self.nc_file = nc.Dataset(self.filename)
        self.time = num2date(self.nc_file["time"], self.nc_file["time"].units)
        current_datenum = date2num(self.current_datetime, self.nc_file["time"].units)
        self.time_index = np.argmin((current_datenum - self.nc_file["time"][:]) ** 2)
        self.file_month = self.current_datetime.month
        self.logger.info("opening new SG file: " + self.filename)
        return


def geo2grid(lat, lon, case):
    case_types = ["get_xy", "get_bl"]
    if case not in case_types:
        raise ValueError(
            "Invalid case type in 3rd argument. Expected one of: %s" % case_types
        )

    a = 6378206.4  # Earth Radius
    fm = 294.97870  # Inverse flattening
    f = 1 / fm  # Flattening
    e = math.sqrt((2 * f) - (f**2))  # Eccentricity
    lon_0 = -45.0  # False origin longitude
    lat_0 = -44.0  # False origin latitude
    lat_1 = -40  # First parallel
    lat_2 = -68  # Second parallel
    x_0 = -175200.0 / 800  # Easting at false origin
    y_0 = 1484800.0 / 800  # Northing at false origin
    dx = 0.8
    imax = 825  # Weddell Sea domain

    rval = len(lat)
    xs = np.empty(rval)
    ys = np.empty(rval)
    for i in range(rval):
        FiF = lat_0 * math.pi / 180
        Fi1 = lat_1 * math.pi / 180
        Fi2 = lat_2 * math.pi / 180
        LamdaF = lon_0 * math.pi / 180

        if case == "get_xy":
            Fi = lat[i] * math.pi / 180
            Lamda = lon[i] * math.pi / 180

        EF = x_0 * dx * 1000
        NF = y_0 * dx * 1000

        if case == "get_bl":
            E = lat[i] * dx * 1000
            N = lon[i] * dx * 1000

        m1 = math.cos(Fi1) / math.sqrt(1 - e**2 * (math.sin(Fi1)) ** 2)
        m2 = math.cos(Fi2) / math.sqrt(1 - e**2 * (math.sin(Fi2)) ** 2)

        t1 = math.tan(math.pi / 4 - Fi1 / 2) / (
            ((1 - e * math.sin(Fi1)) / (1 + e * math.sin(Fi1))) ** (e / 2)
        )
        t2 = math.tan(math.pi / 4 - Fi2 / 2) / (
            ((1 - e * math.sin(Fi2)) / (1 + e * math.sin(Fi2))) ** (e / 2)
        )
        tF = math.tan(math.pi / 4 - FiF / 2) / (
            ((1 - e * math.sin(FiF)) / (1 + e * math.sin(FiF))) ** (e / 2)
        )

        if case == "get_xy":
            t = math.tan(math.pi / 4 - Fi / 2) / (
                ((1 - e * math.sin(Fi)) / (1 + e * math.sin(Fi))) ** (e / 2)
            )

        n = (math.log(m1) - math.log(m2)) / (math.log(t1) - math.log(t2))
        F = m1 / (n * t1**n)
        rF = a * F * tF**n

        if case == "get_xy":
            r = a * F * t**n

        if case == "get_bl":
            rm = np.sign(n) * np.sqrt((E - EF) ** 2 + (rF - (N - NF)) ** 2)
            tm = (rm / (a * F)) ** (1 / n)
            Tetam = np.arctan((E - EF) / (rF - (N - NF)))
            Fim = np.pi / 2 - 2 * np.arctan(Tetam)
            for j in range(1, 9):
                Fi = np.pi / 2 - 2 * np.arctan(
                    tm * ((1 - e * np.sin(Fim)) / (1 + e * np.sin(Fim))) ** (e / 2)
                )
                if np.abs(Fi - Fim) < 1e-7:
                    break
                else:
                    Fim = Fi
            Lamda = Tetam / n + LamdaF
            xs[i] = Fi * 180 / math.pi
            ys[i] = Lamda * 180 / math.pi

        if case == "get_xy":
            Teta = n * (Lamda - LamdaF)
            x = EF + r * math.sin(Teta)
            y = NF + rF - r * math.cos(Teta)
            xs[i] = x / (dx * 1000)
            ys[i] = y / (dx * 1000)

    return xs, ys


def pyproj_geo2grid(x, y, inverse):
    radius = 6378206.4
    x_0 = -175200 / 800
    y_0 = 1484800 / 800
    lon_0 = -45.0  # False origin longitude
    lat_0 = -44.0  # False origin latitude
    lat_1 = -40  # First parallel
    lat_2 = -68  # Second parallel
    proj4string = f"lcc +lat_0={lat_0} +lon_0={lon_0} +lat_1={lat_1} +lat_2={lat_2} +a={radius} +b={radius} +x_0={x_0} +y_0={y_0}"
    proj = pyproj.Proj(proj=proj4string)
    lon, lat = proj(x, y, inverse=inverse)
    return lon, lat


def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance in kilometers between two points
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    c = 2 * np.arcsin(np.sqrt(a))
    r = 6371  # Radius of earth in kilometers. Use 3956 for miles. Determines return value units.
    return c * r
