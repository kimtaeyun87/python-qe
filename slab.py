#-*- coding: utf-8 -*-


import copy
import numpy as np


def getGvector(a, g):
    a = np.array(a)
    g = np.array(g, dtype=int)
    b = getReciprocal(a)
    return np.dot(b, g)


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

