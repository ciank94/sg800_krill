import math
import numpy as np
import pyproj
import netCDF4 as nc
import os
import datetime
from netCDF4 import num2date, date2num
import logging

class SGReader:
    def __init__(self, trajectory_folder, samples_folder, samples_prefix, duration, date_init):
        self.filename = (samples_folder + samples_prefix + str(date_init.year) + "{:02d}".format(date_init.month) + '.nc')
        self.nc_file = nc.Dataset(self.filename)
        self.init_datetime = date_init
        self.duration = duration
        self.current_datetime = self.init_datetime
        self.file_month = self.init_datetime.month
        self.time = num2date(self.nc_file['time'], self.nc_file['time'].units)
        self.save_file_prefix = 'trajectory_' + str(date_init.year) + "{:02d}".format(
            date_init.month) + "{:02d}".format(date_init.day) + '_d' + str(duration)
        self.save_file = self.save_file_prefix + '.nc'
        current_datenum = date2num(self.current_datetime, self.nc_file['time'].units)
        self.time_index = np.argmin((current_datenum-self.nc_file['time'][:])**2)
        # set up logging to file
        self.logger = logging.getLogger(__name__)
        log_filename = 'logger_' + str(date_init.year) + "{:02d}".format(
            date_init.month) + "{:02d}".format(date_init.day) + '_d' + str(duration)
        logging.basicConfig(
            filename=trajectory_folder + log_filename + '.log',
            level=logging.INFO,
            filemode='w',
            #format='[%(asctime)s] {%(pathname)s:%(lineno)d} %(levelname)s - %(message)s',
            format='[%(asctime)s] - %(name)s - %(levelname)s - %(message)s',
            datefmt='%H:%M:%S'
        )
        self.logger.info('opening new SG file: ' + self.filename)

        # mapping between PISCES grid and SINMOD
        self.lon_mapping = samples_folder + 'lon_mapping_' + str(self.current_datetime.year) + '.npy'
        self.lat_mapping = samples_folder + 'lat_mapping_' + str(self.current_datetime.year) + '.npy'
        if not os.path.exists(self.lon_mapping):
            self.read_bio(samples_folder)

        return

    def read_bio(self, samples_folder):
        self.bio_fileprefix = 'CMEMS_SGBIO_Dfull_'
        self.bio_filename = samples_folder + self.bio_fileprefix + str(self.current_datetime.year) + '.nc'
        self.bio_ncfile = nc.Dataset(self.bio_filename)
        lat2 = np.array(self.bio_ncfile['latitude'][:])
        lon2 = np.array(self.bio_ncfile['longitude'][:])
        [x_cmems, y_cmems] = np.meshgrid(lon2, lat2)
        x = np.arange(0, self.nc_file['xc'].shape[0])
        y = np.arange(0, self.nc_file['yc'].shape[0])
        [y_pos, x_pos] = np.meshgrid(y, x)
        [lat1, lon1] = geo2grid(y_pos.flatten(), x_pos.flatten(), 'get_bl')
        self.logger.info('opening biology file ' + self.bio_filename)

        lat_re = lat1.reshape([x.shape[0], y.shape[0]])
        lon_re = lon1.reshape([x.shape[0], y.shape[0]])
        lat_grid = np.zeros(lat_re.shape)
        lon_grid = np.zeros(lon_re.shape)
        for i in range(0, lat_re.shape[0]):
            for j in range(0, lat_re.shape[1]):
                dist_v = haversine(x_cmems, y_cmems, lon_re[i, j], lat_re[i, j])
                idx = np.argmin(dist_v)
                [id_lon, id_lat] = np.unravel_index(idx, dist_v.shape)
                if (id_lon < x_cmems.shape[0]) & (id_lat < x_cmems.shape[1]):
                    lon_grid[i, j] = id_lon
                    lat_grid[i, j] = id_lat

        np.save(self.lon_mapping, lon_grid)
        np.save(self.lat_mapping, lat_grid)
        self.logger.info('saved mapping from ' + self.bio_filename)

        #todo: mapping between times and depth i.e. save a mapping from SINMOD time index->PISCES time index and
        # from SINMOD depth index -> PISCES depth index; then test that it works on server;






    def log_init(self, n, N, x_min, x_max, y_min, y_max, dt, save_step, simulation_steps, save_number):
        self.logger.info('samples_filename = ' + str(self.filename))
        self.logger.info('date_init = ' + str(self.init_datetime))
        self.logger.info('==============================')
        self.logger.info('=========init_params==========')
        self.logger.info('==============================')
        self.logger.info('1) IBM parameters')
        self.logger.info('n = ' + str(n))
        self.logger.info('N = ' + str(N))
        self.logger.info('x_min = ' + str(x_min))
        self.logger.info('x_max = ' + str(x_max))
        self.logger.info('y_min = ' + str(y_min))
        self.logger.info('y_max = ' + str(y_max))
        self.logger.info('==============================')
        self.logger.info('2) Time parameters')
        self.logger.info('duration (days) = ' + str(self.duration))
        self.logger.info('time step (hours) = ' + str(dt))
        self.logger.info('save step (hours) = ' + str(save_step))
        self.logger.info('simulation steps(#) = ' + str(int(simulation_steps)))
        self.logger.info('save number(#) = ' + str(int(save_number)))
        self.logger.info('==============================')
        return

    def update_time(self, dt, test, samples_folder, samples_prefix):
        self.current_datetime += dt
        if self.current_datetime.month != self.file_month:
            self.nc_file.close()
            self.update_file_month(samples_folder=samples_folder, samples_prefix=samples_prefix)
            new_file = True
        else:
            new_file = False

        if test:
            if (self.current_datetime.day != self.time[self.time_index].day) | new_file:
                current_datenum = date2num(self.current_datetime, self.nc_file['time'].units)
                self.time_index = np.argmin((current_datenum - self.nc_file['time'][:]) ** 2)
        else:
            if (self.current_datetime.hour != self.time[self.time_index].hour) | new_file:
                current_datenum = date2num(self.current_datetime, self.nc_file['time'].units)
                self.time_index = np.argmin((current_datenum - self.nc_file['time'][:]) ** 2)
        return

    def update_file_month(self, samples_folder, samples_prefix):
        self.filename = (samples_folder + samples_prefix + str(self.current_datetime.year) +
                    "{:02d}".format(self.current_datetime.month) + '.nc')
        self.nc_file = nc.Dataset(self.filename)
        self.time = num2date(self.nc_file['time'], self.nc_file['time'].units)
        current_datenum = date2num(self.current_datetime, self.nc_file['time'].units)
        self.time_index = np.argmin((current_datenum - self.nc_file['time'][:]) ** 2)
        self.file_month = self.current_datetime.month
        self.logger.info('opening new SG file: ' + self.filename)
        return

def geo2grid(lat, lon, case):
    case_types = ['get_xy', 'get_bl']
    if case not in case_types:
        raise ValueError("Invalid case type in 3rd argument. Expected one of: %s" % case_types)

    a = 6378206.4  # Earth Radius
    fm = 294.97870  # Inverse flattening
    f = 1 / fm  # Flattening
    e = math.sqrt((2 * f) - (f ** 2))  # Eccentricity
    lon_0 = -45.0  # False origin longitude
    lat_0 = -44.0  # False origin latitude
    lat_1 = -40  # First parallel
    lat_2 = -68  # Second parallel
    x_0 = -175200.0/800  # Easting at false origin
    y_0 = 1484800.0/800 # Northing at false origin
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

        if case == 'get_xy':
            Fi = lat[i] * math.pi / 180
            Lamda = lon[i] * math.pi / 180

        EF = x_0 * dx * 1000
        NF = y_0 * dx * 1000

        if case == 'get_bl':
            E = lat[i] * dx * 1000
            N = lon[i] * dx * 1000

        m1 = math.cos(Fi1) / math.sqrt(1 - e ** 2 * (math.sin(Fi1)) ** 2)
        m2 = math.cos(Fi2) / math.sqrt(1 - e ** 2 * (math.sin(Fi2)) ** 2)

        t1 = math.tan(math.pi / 4 - Fi1 / 2) / (((1 - e * math.sin(Fi1)) / (1 + e * math.sin(Fi1))) ** (e / 2))
        t2 = math.tan(math.pi / 4 - Fi2 / 2) / (((1 - e * math.sin(Fi2)) / (1 + e * math.sin(Fi2))) ** (e / 2))
        tF = math.tan(math.pi / 4 - FiF / 2) / (((1 - e * math.sin(FiF)) / (1 + e * math.sin(FiF))) ** (e / 2))

        if case == 'get_xy':
            t = math.tan(math.pi / 4 - Fi / 2) / (((1 - e * math.sin(Fi)) / (1 + e * math.sin(Fi))) ** (e / 2))

        n = (math.log(m1) - math.log(m2)) / (math.log(t1) - math.log(t2))
        F = m1 / (n * t1 ** n)
        rF = a * F * tF ** n

        if case == 'get_xy':
            r = a * F * t ** n

        if case == 'get_bl':
            rm = np.sign(n) * np.sqrt((E - EF) ** 2 + (rF - (N - NF)) ** 2)
            tm = (rm / (a * F)) ** (1 / n)
            Tetam = np.arctan((E - EF) / (rF - (N - NF)))
            Fim = np.pi / 2 - 2 * np.arctan(Tetam)
            for j in range(1, 9):

                Fi = np.pi / 2 - 2 * np.arctan(tm * ((1 - e * np.sin(Fim)) / (1 + e * np.sin(Fim))) ** (e / 2))
                if np.abs(Fi - Fim) < 1e-7:
                    break
                else:
                    Fim = Fi
            Lamda = Tetam / n + LamdaF
            xs[i] = Fi * 180 / math.pi
            ys[i] = Lamda * 180 / math.pi

        if case == 'get_xy':
            Teta = n * (Lamda - LamdaF)
            x = EF + r * math.sin(Teta)
            y = NF + rF - r * math.cos(Teta)
            xs[i] = x / (dx * 1000)
            ys[i] = y / (dx * 1000)

    return xs, ys

def pyproj_geo2grid(x, y, inverse):
    radius = 6378206.4
    x_0 = -175200/800
    y_0 = 1484800/800
    lon_0 = -45.0  # False origin longitude
    lat_0 = -44.0  # False origin latitude
    lat_1 = -40  # First parallel
    lat_2 = -68  # Second parallel
    proj4string = f'lcc +lat_0={lat_0} +lon_0={lon_0} +lat_1={lat_1} +lat_2={lat_2} +a={radius} +b={radius} +x_0={x_0} +y_0={y_0}'
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

