#-*- coding: utf-8 -*-


import numpy as np


def read_magnetic_moment(fpath):
    with open(fpath) as f:
        lines = f.readlines()

    mag_mom = []
    for line in lines:
        if all(x in line for x in ["atomic", "mx", "my", "mz"]):
            _, _, _, _, _, mx, my, mz = line.split()
            mag_mom.append([float(x) for x in [mx, my, mz]])

    return np.array(mag_mom)


def read_const_moment(fpath):
    with open(fpath) as f:
        lines = f.readlines()

    mag_mom = []
    for line in lines:
        if all(x in line for x in ["constrained", "moment"]):
            _, _, _, mx, my, mz = line.split()
            m = [float(x) for x in [mx, my, mz]]
            norm = np.linalg.norm(m)

            if np.isclose(norm, 0.0):
                mag_mom.append([0.0, 0.0, 0.0])
            else:
                mag_mom.append([mx/norm, my/norm, mz/norm])

    return np.array(mag_mom)


