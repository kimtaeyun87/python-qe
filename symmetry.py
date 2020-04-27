#-*- coding: utf-8 -*-


import numpy as np
import slab
import pwscfout as pwo
import spglib


if __name__ == "__main__":

    card = slab.readvestaxtl("test/slab_test.xtl")

    a = card["cell_parameters"]["cell"]
    l = card["atomic_positions"]["label"]
    Omega = np.linalg.det(a)


    if card["atomic_positions"]["option"] == "crystal":
        w = card["atomic_positions"]["position"]
        v = slab.crystaltocartesian(a, w)
    else:
        v = card["atomic_positions"]["position"]
        w = slabmcartesiantocrystal(a, v)

    n = [1]*32 + [2]*32 + [3]*96
    s = [1.0]*8 + [-1.0]*16 + [1.0]*8 + [0.0]*128
    spgcell = (a, w, n, s)
    symdata = spglib.get_symmetry(spgcell)


    R = symdata["rotations"]
    t = symdata["translations"]


    i0 = []
    for n in range(32):
        for i in range(len(R)):
            Rw0t = np.dot(R[i], w[0]) + t[i]
            if np.linalg.norm(Rw0t - w[n]) < 1e-5:
                i0.append(i)
                break


    M = []
    for n in range(32):
        Mrow = []
        for m in range(32):
            wl = np.mod(np.dot(np.linalg.inv(R[i0[n]]), w[m] - t[i0[n]]), [1.0, 1.0, 1.0])
            for l in range(32):
                if np.linalg.norm(wl - w[l]) < 1e-5:
                    Mrow.append(l)
        M.append(Mrow)

    print(M)


    with open("test/nosoc_0deg_n1/scf.out") as f:
        m0 = pwo.readatomicmoment(f.readlines(), 32)

    with open("test/scf.out") as f:
        m = pwo.readatomicmoment(f.readlines(), 32)


    dm = np.array(m) - np.array(m0)
    dmz = dm[:, 2]


    dMz = np.zeros((32, 32))
    for m in range(32):
        for n in range(32):
            dMz[m, n] = dmz[M[m][n]]/0.01

    print(np.linalg.inv(dMz)[0]*13.6)
