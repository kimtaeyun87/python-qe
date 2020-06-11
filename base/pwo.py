#-*- coding: utf-8 -*-


import numpy as np


class PwoutIteration:
    def __init__(self, magnetic_moment):
        self.magnetic_moment = magnetic_moment

    @property
    def magnetic_moment(self):
        return self._magnetic_moment

    @magnetic_moment.setter
    def magnetic_moment(self, magnetic_moment):
        self._magnetic_moment = magnetic_moment


class Pwout:
    def __init__(self, fpath):
        self.iterations = _read_from(fpath)


def _read_from(fpath):
    """Read a pw output file.
    """
    with open(fpath) as f:
        lines = f.readlines()

    for line in lines:
        if 'iteration #' in line:
            _, _, num_iter, _, _, _, _, _ = line.split()
        if 'Tr[ns(na)]' in line:
            _, num_hub, _, _, _, _, _, _, _, _ = line.split()

    num_iter = int(num_iter)
    num_hub = int(num_hub)

    magnetic_moment = _read_magnetic_moment(lines)
    # first few iterations may not be converged sufficiently.
    magnetic_moment = magnetic_moment.reshape(
        (len(magnetic_moment)//num_hub, num_hub, 3))[-num_iter:]

    iterations = []
    for i in range(num_iter):
        iteration = PwoutIteration(magnetic_moment[i])
        iterations.append(iteration)

    return iterations


def _read_magnetic_moment(lines):
    mag_mom = []
    for line in lines:
        if 'atomic mx, my, mz' in line:
            _, _, _, _, _, mx, my, mz = line.split()
            mag_mom.append([float(x) for x in [mx, my, mz]])

    return np.array(mag_mom)


def _read_constrained_moment(lines):
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


