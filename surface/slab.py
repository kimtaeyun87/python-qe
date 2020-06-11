#-*- coding: utf-8 -*-


import copy
import numpy as np


atomsymtonum = {
    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10,
    "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar":18, "K": 19, "Ca": 20,
    "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30,
    "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36, "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40,
    "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
    "Sb": 51, "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57, "Ce": 58, "Pr": 59, "Nd": 60,
    "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65, "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70,
    "Lu": 71, "Hf": 72, "Ta": 72, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80,
    "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84, "At": 85, "Rn": 86, "Fr": 87, "Ra": 88, "Ac": 89, "Th": 90,
    "Pa": 91, "U": 92, "Np": 93, "Pu": 94
    }


def celldm2cell(celldm):
    A, B, C = celldm[:3]
    degtorad = np.pi/180.0
    alpha, beta, gamma = [x*degtorad for x in celldm[3:]]

    a = np.zeros((3, 3))

    a[0, :] = [A, 0.0, 0.0]
    a[1, :] = [B*np.cos(gamma), B*np.sin(gamma), 0.0]
    a[2, :] = [C*np.cos(beta),
               C*(np.cos(alpha) - np.cos(beta)*np.cos(gamma))/np.sin(gamma),
               C*np.sqrt(
                   1.0 + 2.0*np.cos(alpha)*np.cos(beta)*np.cos(gamma) -
                   np.cos(alpha)**2 - np.cos(beta)**2 - np.cos(gamma)**2
               )/np.sin(gamma)]
    return a


def cell2celldm(a):
    A, B, C = [np.linalg.norm(a[i, :]) for i in range(3)]

    cosbc = np.dot(a[1, :]/B, a[2, :]/C)
    cosac = np.dot(a[0, :]/A, a[2, :]/C)
    cosab = np.dot(a[0, :]/A, a[1, :]/B)

    radtodeg = 180.0/np.pi

    alpha = np.arccos(cosbc)*radtodeg
    beta  = np.arccos(cosac)*radtodeg
    gamma = np.arccos(cosab)*radtodeg

    return [A, B, C, alpha, beta, gamma]


def standardizeCell(a):
    celldm = cell2celldm(a)
    return celldm2cell(celldm)


def getReciprocal(a):
    a = np.array(a)
    V = np.linalg.det(a)
    if V < 0: raise ValueError("det(a) < 0")

    b = np.zeros((3, 3))
    b[:, 0] = np.cross(a[1, :], a[2, :])
    b[:, 1] = np.cross(a[2, :], a[0, :])
    b[:, 2] = np.cross(a[0, :], a[1, :])

    return b/V


def getGvector(a, g):
    a = np.array(a)
    g = np.array(g, dtype=int)
    b = getReciprocal(a)
    return np.dot(b, g)


def readVestaxtl(fpath):
    with open(fpath) as f:
        if f.readline().split()[0] != "TITLE":
            raise ValueError
        if f.readline().strip() != "CELL":
            raise ValueError
        celldm = [float(x) for x in f.readline().split()]
        symnum   = int(f.readline().split()[-1])
        symlabel = f.readline().split()[-1]
        if f.readline().strip() != "ATOMS":
            raise ValueError
        line = f.readline()
        if not all([x in line for x in ["NAME", "X", "Y", "Z"]]):
            raise ValueError
        atomlabel = []
        atompos   = []
        while True:
            namepos = f.readline().split()
            if namepos[0] == "EOF": break
            atomlabel.append(namepos[0])
            atompos.append([float(x) for x in namepos[1:]])
        atompos = np.array(atompos)
        atompos = np.mod(atompos, [1.0, 1.0, 1.0])

    cell = celldm2cell(celldm)

    card = {}
    card["cell_parameters"] = {"option": "angstrom", "cell": cell}
    card["atomic_positions"] = {"option": "crystal", "position": atompos, "label": atomlabel}

    return card


def writexsf(f, a, l, v, mode="crystal"):
    """ Write a xsf (xcrysden) file

    Args:

        f: File object

        a (3x3 array): Primitive cell
           [a1x, a1y, a1z]
           [a2x, a2y, a2z]
           [a3x, a3y, a3z]

        l (list): Atom label

        v (nx3 array): Atomic positions in Cartesian coordinate
           [v1x, v1y, v1z]
           [v2x, v2y, v2z]
                  :
                  :
           [vnx, vny, vnz]

        mode: "crystal", "slab", "polymer", "dot". Default to "crystal"

    See also:
        all coordinates in cartesian coordinate
    """
    f.write(f"{mode.upper():s}\n")
    f.write("PRIMVEC\n")

    for i in range(3):
        x, y, z = a[i, :]
        f.write(f"{x:16.8f}{y:16.8f}{z:16.8f}\n")

    f.write("CONVVEC\n")

    for i in range(3):
        x, y, z = a[i, :]
        f.write(f"{x:16.8f}{y:16.8f}{z:16.8f}\n")

    f.write("PRIMCOORD\n")
    f.write(f"{len(v):4d} 1\n")

    for i in range(len(v)):
        x, y, z = v[i, :]
        f.write(f"{atomsymtonum[l[i]]:4d}{x:16.8f}{y:16.8f}{z:16.8f}\n")


def cart2crys(a, v):
    """Convert Cartesian coordinates to crystal coordinates.

    Args:
        a (array): Unitcell.
        v (array): Atomic positions in Cartesian coordinates.

    Returns:
        w (array): Atomic positions in crystal coordinates.

    See also:
        b (array): Reciprocal unitcell.

        w[n, i] = sum_x v[n, i]*b[x, i]

        n: atom index
        i: crystal index
        x: cartesian index
    """
    a = np.array(a)
    v = np.array(v)
    b = getReciprocal(a)
    w = np.dot(v, b)

    return w


def normalize_crystal_coordinate(w):
    w = np.array(w)
    n, m = w.shape
    for i in range(n):
        for j in range(m):
            if np.isclose(w[i, j], 0.0):
                w[i, j] = 0.0
            elif np.isclose(w[i, j], 1.0):
                w[i, j] = 0.0
                #w[i, j] = 1.0
            else:
                w[i, j] = np.remainder(w[i, j], 1.0)
    return w


def crys2cart(a, w):
    """Convert crystal coordinates to Cartesian coordinates.

    Args:
        a (array): Unitcell.
        w (array): Atomic positions in crystal coordinate.

    Returns:
        v (array): Atomic positions in Cartesian coordinate.

    See also:
        v[n, x] = sum_i w[n, i]*a[i, x]

        n: atom index
        i: crystal index
        x: cartesian index
    """
    a = np.array(a)
    w = np.array(w)
    v = np.dot(w, a)

    return v


def getSurfaceCell(a, g):
    """Get new unitcell vectors that are compatible with the given miller index.

    Args:
        a (array): Unitcell.
        g (array): Miller index.

    Returns:
        asurf (array-like): A surface unitcell.
    """
    a = np.array(a)
    g = np.array(g).astype(int)
    g = g//np.gcd.reduce(g)

    asurf = np.zeros((3, 3))

    if np.linalg.det(a) < 0: raise ValueError("det(a) < 0")

    if g[0] == 0 and g[1] == 0 and g[2] == 0:
        raise ValueError

    elif g[0] == 0 and g[1] == 0:
        asurf[0, :] = a[0, :]
        asurf[1, :] = a[1, :]
        asurf[2, :] = a[2, :]

    elif g[1] == 0 and g[2] == 0:
        asurf[0, :] = a[1, :]
        asurf[1, :] = a[2, :]
        asurf[2, :] = a[0, :]

    elif g[2] == 0 and g[0] == 0:
        asurf[0, :] = a[2, :]
        asurf[1, :] = a[0, :]
        asurf[2, :] = a[1, :]

    elif g[0] == 0:
        asurf[0, :] = a[0, :]
        asurf[1, :] = (a[1, :]/g[1] - a[2, :]/g[2])*np.lcm(g[1], g[2])

        if np.abs(g[1])/np.linalg.norm(a[1, :]) > np.abs(g[2])/np.linalg.norm(a[2, :]):
            asurf[2, :] = a[1, :]
        else:
            asurf[2, :] = a[2, :]

    elif g[1] == 0:
        asurf[0, :] = a[1, :]
        asurf[1, :] = (a[2, :]/g[2] - a[0, :]/g[0])*np.lcm(g[2], g[0])

        if np.abs(g[2])/np.linalg.norm(a[2, :]) > np.abs(g[0])/np.linalg.norm(a[0, :]):
            asurf[2, :] = a[2, :]
        else:
            asurf[2, :] = a[0, :]

    elif g[2] == 0:
        asurf[0, :] = a[2, :]
        asurf[1, :] = (a[0, :]/g[0] - a[1, :]/g[1])*np.lcm(g[0], g[1])

        if np.abs(g[0])/np.linalg.norm(a[0, :]) > np.abs(g[1])/np.linalg.norm(a[1, :]):
            asurf[2, :] = a[0, :]
        else:
            asurf[2, :] = a[1, :]

    else:
        asurf[0, :] = (a[0, :]/g[0] - a[1, :]/g[1])*np.lcm(g[0], g[1])
        asurf[1, :] = (a[1, :]/g[1] - a[2, :]/g[2])*np.lcm(g[1], g[2])

        x = [np.abs(g[i])/np.linalg.norm(a[i, :]) for i in range(3)]
        asurf[2, :] = a[np.argmax(x), :]

    if np.linalg.det(asurf) < 0: raise ValueError("det(asurf) < 0")

    return asurf


def collectAtom(a, p, l, v, ns=1):
    """
    """
    a = np.array(a)
    p = np.array(p)
    l = copy.deepcopy(l)
    v = np.array(v)

    p1 = np.array(np.amin(p, axis=0), dtype=int) - ns
    p2 = np.array(np.amax(p, axis=0), dtype=int) + ns

    L = []
    V = []

    for i in range(p1[0], p2[0]):
        for j in range(p1[1], p2[1]):
            for k in range(p1[2], p2[2]):
                L += l
                V += (v + i*a[0, :] + j*a[1, :] + k*a[2, :]).tolist()

    V = np.array(V)
    return L, V



def reduceAtom(a, p, L, V):
    """
    """
    a = np.array(a)
    p = np.array(p)
    L = copy.deepcopy(L)
    V = np.array(V)

    p1 = np.amin(p, axis=0).astype(float)
    p2 = np.amax(p, axis=0).astype(float)

    l = []
    v = []

    W = cart2crys(a, V)

    for Li, Wi, Vi in list(zip(L, W, V)):
        Wi_ = np.zeros(3)
        for j in range(3):
            if np.isclose(p1[j], Wi[j]):
                Wi_[j] = p1[j]
            elif np.isclose(p2[j], Wi[j]):
                Wi_[j] = p2[j]
            else:
                Wi_[j] = Wi[j]

        if all(p1 <= Wi_) and all(Wi_ < p2):
            l.append(Li)
            v.append(crys2cart(a, Wi_))

    v = np.array(v)
    return l, v


def convertCell(a1, l1, v1, a2):
    """
    """
    a1 = np.array(a1)
    l1 = copy.deepcopy(l1)
    v1 = np.array(v1)
    a2 = np.array(a2)

    b1 = getReciprocal(a1)
    p21 = np.dot(a2, b1)

    L, V = collectAtom(a1, p21, l1, v1)

    p = [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]
    l2, v2 = reduceAtom(a2, p, L, V)

    return a2, l2, v2


def transformCoordinate(a1, l1, v1, a2):
    """
    """
    a1 = np.array(a1)
    l1 = copy.deepcopy(l1)
    v1 = np.array(v1)
    a2 = np.array(a2)

    l2 = l1
    r12 = np.dot(np.linalg.inv(a1), a2)
    v2 = np.dot(v1, r12)

    return a2, l2, v2


def changeBoundary(a, l, v, p):
    """
    """
    a = np.array(a)
    l = copy.deepcopy(l)
    v = np.array(v)
    p = np.array(p)

    L, V = collectAtom(a, p, l, v)
    l, v = reduceAtom(a, p, L, V)
    return a, l, v


def makeSlab(a, l, v, wmin, wmax, woff, dgap):
    """
    """
    a = np.array(a)
    l = copy.deepcopy(l)
    v = np.array(v)

    p = np.array([[0, 0, wmin], [1, 1, wmax]])
    t = np.array([0, 0, woff])

    a, l, v = changeBoundary(a, l, v, p+t)

    w = cart2crys(a, v)
    w -= t
    v = crys2cart(a, w)

    dsurf = a[2, 2]
    a[2, :]  = np.array([0.0, 0.0, (wmax - wmin)*dsurf + dgap])
    return a, l, v


if __name__ == "__main__":

    np.set_printoptions(formatter={
        "float": lambda x: "{:16.8f}".format(x),
        "int": lambda x: "{:4d}".format(x)})

    # 1. Build a 3QL slab of Bi2Se3
    import Bi2Se3

    aprim = Bi2Se3.rhombo(9.841, 24.27*np.pi/180.0)
    u1, u2 = (0.4, 0.21)

    wprim = np.array([
        [0.0, 0.0, 0.0],
        [ u2,  u2,  u2],
        [-u2, -u2, -u2],
        [ u1,  u1,  u1],
        [-u1, -u1, -u1]
    ], dtype=float)

    wprim = np.mod(wprim, [1.0, 1.0, 1.0])

    vprim = crys2cart(aprim, wprim)
    lprim = ["Se", "Se", "Se", "Bi", "Bi"]

    # Unitcell
    Xprim = (aprim, lprim, vprim)

    # Convert to surface cell
    gsurf = [1, 1, 1]
    Gsurf = getGvector(aprim, gsurf)
    asurf = getSurfaceCell(aprim, gsurf)
    Xsurf = convertCell(*Xprim, asurf)

    # Coordinate transformation
    Xsurf = transformCoordinate(*Xsurf, standardizeCell(asurf))

    # Make slab
    wmin, wmax, woff, gap = -1.5, 1.5, 0.0, 20.0
    Xslab = makeSlab(*Xsurf, wmin, wmax, woff, gap)

    with open("tmp.xsf", 'w') as f:
        writexsf(f, *Xslab, mode="slab")
    import subprocess
    subprocess.call(["/home/kimtaeyun/xcrysden/1.6.2/xcrysden", "--xsf", "tmp.xsf"])
    subprocess.call(["rm", "-f", "/tmp/tmp.xsf"])


    # 2. Build a BL NiPS3
    card = readVestaxtl("test/slab_test.xtl")
    a = card["cell_parameters"]["cell"]
    w = card["atomic_positions"]["position"]
    l = card["atomic_positions"]["label"]
    v = crys2cart(a, w)
    

    Xsurf = (a, l, v)
    wmin, wmax, woff, dgap = -1.0, 1.0, 0.5, 20.0
    Xslab = makeSlab(*Xsurf, wmin, wmax, woff, dgap)

    with open("tmp.xsf", 'w') as f:
        writexsf(f, *Xslab, mode="slab")
    import subprocess
    subprocess.call(["/home/kimtaeyun/xcrysden/1.6.2/xcrysden", "--xsf", "tmp.xsf"])
    subprocess.call(["rm", "-f", "/tmp/tmp.xsf"])

