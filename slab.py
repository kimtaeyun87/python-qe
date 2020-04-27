#-*- coding: utf-8 -*-


import numpy as np
#import pwscfin as pwi


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
    b = np.zeros((3, 3))

    V = np.linalg.det(a)
    b[:, 0] = np.cross(a[1, :], a[2, :])/V
    b[:, 1] = np.cross(a[2, :], a[0, :])/V
    b[:, 2] = np.cross(a[0, :], a[1, :])/V

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

    V = np.linalg.det(a)
    if V < 0: raise ValueError("det(a) < 0")

    b = np.zeros((3, 3))
    b[:, 0] = np.cross(a[1, :], a[2, :])/V
    b[:, 1] = np.cross(a[2, :], a[0, :])/V
    b[:, 2] = np.cross(a[0, :], a[1, :])/V

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


if __name__ == "__main__":

    np.set_printoptions(formatter={
        "float": lambda x: "{:15.8f}".format(x),
        "int": lambda x: "{:4d}".format(x)})

    #nml, card = pwi.readpwscfin("test/nosoc_0deg_n1/scf.nml")
    card = readvestaxtl("test/slab_test.xtl")


    a = card["cell_parameters"]["cell"]
    l = card["atomic_positions"]["label"]
    Omega = np.linalg.det(a)


    if card["atomic_positions"]["option"] == "crystal":
        w = card["atomic_positions"]["position"]
        v = crystaltocartesian(a, w)
    else:
        v = card["atomic_positions"]["position"]
        w = cartesiantocrystal(a, v)


    b = np.zeros((3, 3))
    b[:, 0] = np.cross(a[1, :], a[2, :])/Omega
    b[:, 1] = np.cross(a[2, :], a[0, :])/Omega
    b[:, 2] = np.cross(a[0, :], a[1, :])/Omega


    """
    G: sum_i b[:, i]*h[i]
    s: A unitcell whose first and second rows are on the [h1, h2, h3] plane
    """
    G, s = getsurface(a, [0, 0, 1])
    """
    p: 
    """
    p = np.dot(s, b).astype(int)

    """
    find all the atoms that are contained in the parallelopipe formed 
    by the six vertices:
    min(p[:, 0])*a[0, :], min(p[:, 1])*a[1, :], min(p[:, 2])*a[2, :],
    max(p[:, 0])*a[0, :], max(p[:, 2])*a[1, :], max(p[:, 2])*a[2, :].
    """
    W_ = []
    L_ = []
    for i in range(np.min(p[:, 0]), np.max(p[:, 0])):
        for j in range(np.min(p[:, 1]), np.max(p[:, 1])):
            for k in range(np.min(p[:, 2]), np.max(p[:, 2])):
                W_ += (w + np.array([i, j, k], dtype=float)).tolist()
                L_ += l
    V_ = crystaltocartesian(a, W_)


    """
    remove all the atoms that are located outside the unitcell
    previously obtained
    """
    W = []
    L = []
    for i, x in enumerate(cartesiantocrystal(s, V_)):
        mask = np.isclose(x, [0.0, 0.0, 0.0])
        y = [0.0 if mask[i] else x[i] for i in range(3)]
        if all(0 <= np.array(y)) and all(np.array(y) < 1.0):
            W.append(y)
            L.append(L_[i])
    V = crystaltocartesian(s, W)


    """
    transformation of the coordinate axes
    s[0, :] --> s'[0, :] parallel to x^, s'[0, :] perpendicular to z^
    s[1, :] --> s'[1, :] in x^, s'[1, :] perpendicular to z^
    """
    sp = celldmtocell(celltocelldm(s))
    """
    coordinate transformation matrix
    r[x, x'] = sum_i s^(-1)[x, i]*s'[i, x']
    """
    r =np.dot(np.linalg.inv(s), sp)
    """
    V'[i, x'] = sum_x V[i, x]*r[x, x']
    """
    Vp = np.dot(V, r)
    Wp = cartesiantocrystal(sp, Vp)
    Lp = L


    """
    set window parameters (along the z' direction)
    """
    winmin, winmax, wintrans, vacuum = -0.5, 0.5, 0.0, 20.0

    """
    find all the atoms within the slab of width (winmax - winmax)/|G|
    note that 1/|G| is the width between adjacent [h1, h2, h3] planes
    """
    Wp_ = []
    Lp_ = []
    for k in range(int(np.floor(winmin)), int(np.ceil(winmax))):
        Wp_ += (Wp + np.array([0, 0, k], dtype=float)).tolist()
        Lp_ += Lp
    Vp_ = crystaltocartesian(sp, Wp_)


    """
    remove all the atoms outside the window and apply (if given) a translation
    """
    Wp = []
    Lp = []
    for i, x in enumerate(Wp_):
        if np.isclose(x[2], winmin):
            y = np.array([x[0], x[1], winmin])
        else:
            y = np.array(x)
        if winmin <= y[2] and y[2] < winmax:
            Wp.append(y + np.array([0.0, 0.0, wintrans]))
            Lp.append(Lp_[i])
    Vp = crystaltocartesian(sp, Wp)
    sp[2, :] = np.array([0.0, 0.0, (winmax - winmin)/np.linalg.norm(G)])


    """
    finally, insert a vacuum layer in the z' direction
    """
    sp[2, :] += np.array([0.0, 0.0, vacuum])


    """
    write in xsf format for testing
    """
    #print("CRYSTAL")
    print("SLAB")
    print("PRIMVEC")
    for i in range(3):
        x, y, z = sp[i, :]
        print(f"{x:16.8f}{y:16.8f}{z:16.8f}")
    print("CONVVEC")
    for i in range(3):
        x, y, z = sp[i, :]
        print(f"{x:16.8f}{y:16.8f}{z:16.8f}")
    print("PRIMCOORD")
    print(f"{len(Lp):4d} 1")
    symtonum = {"Ni": 28, "P": 15, "S": 16}
    for i in range(len(Lp)):
        x, y, z = Vp[i, :]
        print(f"{symtonum[Lp[i]]:4d}{x:16.8f}{y:16.8f}{z:16.8f}")
