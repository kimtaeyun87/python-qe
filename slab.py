#-*- coding: utf-8 -*-


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
    "Ti": 81, "Pb": 82, "Bi": 83, "Po": 84, "At": 85, "Rn": 86, "Fr": 87, "Ra": 88, "Ac": 89, "Th": 90,
    "Pa": 91, "U": 92, "Np": 93, "Pu": 94
    }


def celldmtocell(celldm):
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


def celltocelldm(s):
    A, B, C = [np.linalg.norm(s[i, :]) for i in range(3)]

    cosbc = np.dot(s[1, :]/B, s[2, :]/C)
    cosac = np.dot(s[0, :]/A, s[2, :]/C)
    cosab = np.dot(s[0, :]/A, s[1, :]/B)

    radtodeg = 180.0/np.pi

    alpha = np.arccos(cosbc)*radtodeg
    beta  = np.arccos(cosac)*radtodeg
    gamma = np.arccos(cosab)*radtodeg

    return [A, B, C, alpha, beta, gamma]


def getreciprocal(a):
    a = np.array(a)
    V = np.linalg.det(a)
    if V < 0: raise ValueError("det(a) < 0")

    b = np.zeros((3, 3))
    b[:, 0] = np.cross(a[1, :], a[2, :])
    b[:, 1] = np.cross(a[2, :], a[0, :])
    b[:, 2] = np.cross(a[0, :], a[1, :])

    return b/V


def readvestaxtl(fpath):
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

    cell = celldmtocell(celldm)
    out = {}
    out["cell_parameters"] = {"option": "angstrom", "cell": cell}
    out["atomic_positions"] = {"option": "crystal", "position": atompos, "label": atomlabel}

    return out


def writexsf(f, a, c, l, v, mode="crystal"):
    """
        write structure in xsf format
        a: primitive cell
        c: conventional cell
        l: atom labels
        v: atomic coordinates within primitive cell
        mode: "crystal", "slab", "polymer", "dot"
        all coordinates in cartesian coordinate
    """
    f.write(f"{mode.upper():s}\n")
    f.write("PRIMVEC\n")

    for i in range(3):
        x, y, z = a[i, :]
        f.write(f"{x:16.8f}{y:16.8f}{z:16.8f}\n")

    f.write("CONVVEC\n")

    for i in range(3):
        x, y, z = c[i, :]
        f.write(f"{x:16.8f}{y:16.8f}{z:16.8f}\n")

    f.write("PRIMCOORD\n")
    f.write(f"{len(v):4d} 1\n")

    for i in range(len(v)):
        x, y, z = v[i, :]
        f.write(f"{atomsymtonum[l[i]]:4d}{x:16.8f}{y:16.8f}{z:16.8f}\n")


def cartesiantocrystal(a, v):
    """
        a[i, x]: unitcell
        v[n, x]: atomic coordinates in cartesian basis

        b[x, i]: reciprocal unitcell

        n: atomic index
        i: crystal coordinate index
        x: cartesian coordinate index

        w[n, i]: atomic coordinates in crystal basis
        w[n, i] = sum_x v[n, i]*b[x, i]
    """
    a = np.array(a)
    v = np.array(v)
    b = getreciprocal(a)

    w = np.dot(v, b)

    return w


def crystaltocartesian(a, w):
    """
        a[i, x]: unitcell
        w[n, i]: atomic coordinates in crystal basis

        n: atomic index
        i: crystal coordinate index
        x: cartesian coordinate index

        v[n, x]: atomic coordinates in cartesian basis
        v[n, x] = sum_i w[n, i]*a[i, x]
    """
    a = np.array(a)
    w = np.array(w)

    v = np.dot(w, a)

    return v


def getsurface(a, g):
    """
        a: unitcell in cartesian coordinate
        g: miller index defining the surface
        s: supercell in cartesian coordinate
    """
    a = np.array(a)
    g = np.array(g, dtype=np.int32)
    g = g//np.gcd.reduce(g)
    s = np.zeros((3, 3))

    if np.linalg.det(a) < 0: raise ValueError("det(a) < 0")

    b = getreciprocal(a)

    G = np.dot(b, g)

    if g[0] == 0 and g[1] == 0 and g[2] == 0:
        raise ValueError

    elif g[0] == 0 and g[1] == 0:
        s[0, :] = a[0, :]
        s[1, :] = a[1, :]
        s[2, :] = a[2, :]

    elif g[1] == 0 and g[2] == 0:
        s[0, :] = a[1, :]
        s[1, :] = a[2, :]
        s[2, :] = a[0, :]

    elif g[2] == 0 and g[0] == 0:
        s[0, :] = a[2, :]
        s[1, :] = a[0, :]
        s[2, :] = a[1, :]

    elif g[0] == 0:
        s[0, :] = a[0, :]
        s[1, :] = (a[1, :]/g[1] - a[2, :]/g[2])*np.lcm(g[1], g[2])

        if np.abs(g[1])/np.linalg.norm(a[1, :]) > np.abs(g[2])/np.linalg.norm(a[2, :]):
            s[2, :] = a[1, :]
        else:
            s[2, :] = a[2, :]

    elif g[1] == 0:
        s[0, :] = a[1, :]
        s[1, :] = (a[2, :]/g[2] - a[0, :]/g[0])*np.lcm(g[2], g[0])

        if np.abs(g[2])/np.linalg.norm(a[2, :]) > np.abs(g[0])/np.linalg.norm(a[0, :]):
            s[2, :] = a[2, :]
        else:
            s[2, :] = a[0, :]

    elif g[2] == 0:
        s[0, :] = a[2, :]
        s[1, :] = (a[0, :]/g[0] - a[1, :]/g[1])*np.lcm(g[0], g[1])

        if np.abs(g[0])/np.linalg.norm(a[0, :]) > np.abs(g[1])/np.linalg.norm(a[1, :]):
            s[2, :] = a[0, :]
        else:
            s[2, :] = a[1, :]

    else:
        s[0, :] = (a[0, :]/g[0] - a[1, :]/g[1])*np.lcm(g[0], g[1])
        s[1, :] = (a[1, :]/g[1] - a[2, :]/g[2])*np.lcm(g[1], g[2])

        x = [np.abs(g[i])/np.linalg.norm(a[i, :]) for i in range(3)]
        s[2, :] = a[np.argmax(x), :]

    if np.linalg.det(s) < 0: raise ValueError("det(s) < 0")

    return G, s


def expandcell(a, p, l, w, ns=1):
    """
    expand unitcell to a much larger cell containing
    a[i, x]: unitcell
    p[j, i]: coefficients of new unitcell in terms of unitcell
    l[n]   : labels of atoms within a[i, x]
    w[n, x]: atomic coordinates within a[i, x]
    L[n]   : labels of atoms within larger cell
    W[n, x]: atomic coordinates within larger cell
    """

    """
    find all the atoms that are contained in the box formed
    by the following six vertices:
    1. min(p[:, 0])*a[0, :]
    2. min(p[:, 1])*a[1, :]
    3. min(p[:, 2])*a[2, :]
    4. max(p[:, 0])*a[0, :]
    5. max(p[:, 2])*a[1, :]
    6. max(p[:, 2])*a[2, :]
    """

    p1 = np.amin(p, axis=0).astype(int)
    p2 = np.amax(p, axis=0).astype(int)

    W = []
    L = []

    for i in range(p1[0]-ns, p2[0]+ns):
        for j in range(p1[1]-ns, p2[1]+ns):
            for k in range(p1[2]-ns, p2[2]+ns):
                W += (w + np.array([i, j, k], dtype=float)).tolist()
                L += l

    W = np.array(W)

    return L, W


def shrinkcell(s, p, L, W):
    w = []
    l = []

    p1 = np.amin(p, axis=0).astype(float)
    p2 = np.amax(p, axis=0).astype(float)

    for li, wi in list(zip(L, W)):
        mask = np.isclose(wi, p1)
        wi_clip = np.array([0.0 if mask[x] else wi[x] for x in range(3)])

        if all(p1 <= wi_clip) and all(wi_clip < p2):
            w.append(wi_clip)
            l.append(li)

    w = np.array(w)

    return l, w


if __name__ == "__main__":

    np.set_printoptions(formatter={
        "float": lambda x: "{:15.8f}".format(x),
        "int": lambda x: "{:4d}".format(x)})

    # building a slab of Bi2Se3 from scratch (bulk structure)

    # A0: primitive unitcell 
    # W0: atomic coordinates within A0 in crystal coordinates
    # V0: atomic coordinates within A0 in cartesian coordinates
    # L0: atom labels within A0
    import Bi2Se3

    A0 = Bi2Se3.rhombo(9.841, 24.27*np.pi/180.0)
    u1, u2 = (0.4, 0.21)

    W0 = np.zeros((5, 3))
    W0[0, :] = np.array([0.0, 0.0, 0.0])
    W0[1, :] = np.array([u2, u2, u2])
    W0[2, :] = -W0[1, :]
    W0[3, :] = np.array([u1, u1, u1])
    W0[4, :] = -W0[3, :]
    W0 = np.mod(np.array(W0), [1.0, 1.0, 1.0])

    V0 = crystaltocartesian(A0, W0)
    L0 = ["Se", "Se", "Se", "Bi", "Bi"]

    """
    G[x]: sum_i b[x, i]*h[i]
    s[i, x]: new unitcell for generating a surface perpendicular G[x]
    """
    # gs : miller index of the surface normal vector of a slab to be created
    # Gs : surface normal vector in cartesian coordinates
    # A0s: new unitcell for creation of a slab model
    gs = [1, 1, 1]
    Gs, A0s = getsurface(A0, gs)

    B0 = getreciprocal(A0)
    box = np.dot(A0s, B0)

    # W1: 
    # V1: 
    # L1: 
    L1, W1 = expandcell(A0, box, L0, W0, ns=1)
    V1 = crystaltocartesian(A0, W1)


    """
    remove all the atoms that are located outside the unitcell
    previously obtained
    """
    # W1s:
    # V1s:
    # L1s:
    W1s = cartesiantocrystal(A0s, V1)
    L1s = L1

    box = [[0.0, 0.0, 0.0],
           [1.0, 1.0, 1.0]]

    L0s, W0s = shrinkcell(A0s, box, L1s, W1s)
    V0s = crystaltocartesian(A0s, W0s)


    """
    transformation of the coordinate axes
    A0s[0, :] --> A0t[0, :] on xy plane, parallel to x
    A0s[1, :] --> A0t[1, :] on xy plane
    A0s[2, :] --> A0t[2, :] parallel to z (automatically)
    """
    A0t = celldmtocell(celltocelldm(A0s))

    """
    coordinate transformation matrix
    Rst[x, x']  = sum_i A0s^(-1)[x, i]*A0t[i, x']
    V0t[i, x'] = sum_x V0s[i, x]*Rst[x, x']
    """
    Rst = np.dot(np.linalg.inv(A0s), A0t)
    V0t = np.dot(V0s, Rst)
    W0t = cartesiantocrystal(A0t, V0t)
    L0t = L0s


    """
    set window parameters (along the z' direction)
    """
    wmin, wmax, wtrans, gap = -1.5, 1.5, 0.0, 10.0

    """
    find all the atoms within the slab of width (winmax - winmax)/|G|
    note that 1/|G| is the width between adjacent [h1, h2, h3] planes
    """
    box = [[0, 0, np.floor(wmin)],
           [1, 1, np.ceil(wmax)]]
    L1t, W1t = expandcell(A0t, box, L0t, W0t, ns=0)


    """
    remove all the atoms outside the window and apply (if given) a translation
    """
    box = [[0.0, 0.0, wmin],
           [1.0, 1.0, wmax]]
    L0t, W0t = shrinkcell(A0t, box, L1t, W1t)
    V0t = crystaltocartesian(A0t, W0t)

    A0t[2, :] *= wmax - wmin


    """
    finally, insert a vacuum layer in the z' direction and set s'[2, :]
    to be parallel to z' direction
    """
    A0t[2, :] = np.array([0.0, 0.0, (wmax - wmin)/np.linalg.norm(Gs)])
    A0t[2, :] += np.array([0.0, 0.0, gap])

    with open("/tmp/tmp.xsf", 'w') as f:
        writexsf(f, A0t, A0t, L0t, V0t, mode="slab")

    import subprocess
    subprocess.call(["/home/kimtaeyun/xcrysden/1.6.2/xcrysden", "--xsf", "/tmp/tmp.xsf"])
    subprocess.call(["rm", "-f", "/tmp/tmp.xsf"])

