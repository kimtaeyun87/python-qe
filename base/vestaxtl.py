import numpy as np

from base.lattice import Lattice


def readVestaxtl(fpath):

    with open(fpath) as f:
        if f.readline().split()[0] != "TITLE":
            raise ValueError

        if f.readline().strip() != "CELL":
            raise ValueError

        cell_params = [float(x) for x in f.readline().split()]

        _ = int(f.readline().split()[-1])
        _ = f.readline().split()[-1]

        if f.readline().strip() != "ATOMS":
            raise ValueError

        line = f.readline()
        if not all([x in line for x in ["NAME", "X", "Y", "Z"]]):
            raise ValueError

        names = []
        positions = []

        while True:
            a = f.readline().split()
            if a[0] == "EOF": break

            names.append(a[0])
            positions.append([float(x) for x in a[1:]])

        positions = np.array(positions)
        positions = np.mod(positions, [1.0, 1.0, 1.0])

    unitcell = Lattice.get_unitcell_from_cell_params(cell_params)

    card = {}

    card["cell_parameters"] = {
        "option": "angstrom", "unitcell": unitcell}
    card["atomic_positions"] = {
        "option": "crystal", "position": positions, "name": names}

    return card