import math
import numpy as np
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