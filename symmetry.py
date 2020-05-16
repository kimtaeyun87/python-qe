#-*- coding: utf-8 -*-


import numpy as np
import slab
import pwscfout as pwo
import spglib


# python implementation as reference
def get_map_IJ_to_I0(Nw, w, NR, R, t):
    """
    """
    map_i0 = []

    for i in range(Nw):
        x = []
        for n in range(NR):
            wi = np.dot(R[n], w[0]) + t[n]

            wi = np.mod(wi, [1.0, 1.0, 1.0])
            for k in range(3):
                if np.isclose(wi[k], 1.0, atol=1e-5):
                    wi[k] = 0.0

            y = np.linalg.norm(w[i] - wi)
            if np.isclose(y, 0.0, atol=1e-5):
                x.append(n)
        map_i0.append(x)

    map_IJ = np.empty((Nw, Nw), dtype=object)

    for I in range(Nw):
        for J in range(Nw):
            x = []
            for n in map_i0[I]:
                RI = R[n]
                tI = t[n]
                wi = np.dot(np.linalg.inv(RI), w[J] - tI)

                wi = np.mod(wi, [1.0, 1.0, 1.0])
                for k in range(3):
                    if np.isclose(wi[k], 1.0, atol=1e-5):
                        wi[k] = 0.0

                for i in range(Nw):
                    y = np.linalg.norm(w[i] - wi)
                    if np.isclose(y, 0.0, atol=1e-5):
                        x.append(i)
                        break
            map_IJ[I, J] = x

    return map_IJ


if __name__ == "__main__":

    np.set_printoptions(formatter={
        "float": lambda x: "{:16.8f}".format(x),
        "int": lambda x: "{:4d}".format(x)})

    import pwscfin as pwi
    nml, card = pwi.readpwscfin("test2/collinear_reference/scf.nml")
    # nml, card = pwi.readpwscfin("test/collinear_reference/scf.nml")
    # card = slab.readVestaxtl("test2/NiPS3_HT_BL.xtl")

    a = card["cell_parameters"]["cell"]
    l = card["atomic_positions"]["label"]
    if card["atomic_positions"]["option"] == "crystal":
        w = card["atomic_positions"]["position"]
    else:
        v = card["atomic_positions"]["position"]
        w = slab.cart2crys(a, v)
    w = slab.normalize_crystal_coordinate(w) 

    # n = [1]*32 + [2]*32 + [3]*96
    # s = [1.0]*8 + [-1.0]*16 + [1.0]*8 + [0.0]*128
    # Nhub = 32

    n = [1]*64 + [2]*64 + [3]*192
    s = [1.0]*16 + [-1.0]*32 + [1.0]*16 + [0.0]*256
    Nhub = 64

    spgcell = (a, w, n, s)
    symdata = spglib.get_symmetry(spgcell)


    R = symdata["rotations"]
    t = symdata["translations"]


    import ctypes as ct
    import numpy.ctypeslib as npct
    libsym = ct.CDLL('./libsymmetry.so', mode=ct.RTLD_GLOBAL)
    libsym.get_equiv_atom_map.argtypes = [
        ct.c_int,
        npct.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS'),
        ct.c_int,
        npct.ndpointer(dtype=np.int32, ndim=3, flags='C_CONTIGUOUS'),
        npct.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS')
    ]
    libsym.get_equiv_atom_map.restype = ct.POINTER(ct.POINTER(ct.c_int))

    equiv_atom_map = libsym.get_equiv_atom_map(Nhub, w, len(R), R, t)
    map_count = npct.as_array(equiv_atom_map[0], (Nhub,))
    max_count = max(map_count)
    map_IJ_to_I0 = npct.as_array(equiv_atom_map[1], (Nhub,Nhub,max_count))


    # map_IJ_to_I0 = get_map_IJ_to_I0(Nhub, w, len(R), R, t)


    m0 = pwo.readatomicmoment('test2/collinear_reference/scf.out', Nhub)
    m1 = pwo.readatomicmoment('test2/collinear_perturbation_magnitude_1/scf.out', Nhub)
    m2 = pwo.readatomicmoment('test2/collinear_perturbation_magnitude_2/scf.out', Nhub)
    m3 = pwo.readatomicmoment('test2/collinear_perturbation_magnitude_3/scf.out', Nhub)
    m4 = pwo.readatomicmoment('test2/collinear_perturbation_magnitude_4/scf.out', Nhub)
    m5 = pwo.readatomicmoment('test2/collinear_perturbation_magnitude_5/scf.out', Nhub)
    m6 = pwo.readatomicmoment('test2/collinear_perturbation_magnitude_6/scf.out', Nhub)
    m7 = pwo.readatomicmoment('test2/collinear_perturbation_magnitude_7/scf.out', Nhub)
    m8 = pwo.readatomicmoment('test2/collinear_perturbation_magnitude_8/scf.out', Nhub)

    p = [0.0, 0.01, 0.005, -0.005, -0.01]
    dx = np.linspace(-0.01, 0.01, 100)

    m_lx = np.zeros((5, Nhub*3))
    m_lx[0] = m0.flatten()
    m_lx[1] = m1.flatten()
    m_lx[2] = m3.flatten()
    m_lx[3] = m5.flatten()
    m_lx[4] = m7.flatten()

    m_lz = np.zeros((5, Nhub*3))
    m_lz[0] = m0.flatten()
    m_lz[1] = m2.flatten()
    m_lz[2] = m4.flatten()
    m_lz[3] = m6.flatten()
    m_lz[4] = m8.flatten()

    dm_dlx = np.reshape(np.polyfit(p, m_lx, deg=4)[-2, :], (Nhub, 3))
    dm_dlz = np.reshape(np.polyfit(p, m_lz, deg=4)[-2, :], (Nhub, 3))

    dm_dl = np.zeros((Nhub*2, Nhub*2))

    for I in range(Nhub):
        for J in range(Nhub):
            N = len(map_IJ_to_I0[I, J])

            xx = 0.0
            zx = 0.0

            xz = 0.0
            zz = 0.0

            for n in range(N):
                i = map_IJ_to_I0[I, J, n]

                xx += dm_dlx[i, 0]
                zx += dm_dlx[i, 2]

                xz += dm_dlz[i, 0]
                zz += dm_dlz[i, 2]

            dm_dl[I   , J   ] = xx/N
            dm_dl[I+Nhub, J   ] = zx/N
            dm_dl[I   , J+Nhub] = xz/N
            dm_dl[I+Nhub, J+Nhub] = zz/N

    dl_dm = np.linalg.inv(dm_dl)

    Lambda = 0.1
    d2F_dm2 = -dl_dm - Lambda*np.eye(Nhub*2)

    print("d2F/dm2 Matrix (eV)")
    for I in range(Nhub*2):
        for J in range(Nhub*2):
            print(f"({I:2d},{J:2d}): {d2F_dm2[I, J]*13.6:16.8f}")
        print('')

    JI0 = d2F_dm2[:Nhub,:Nhub][0]
    Ry2meV = 13.6*1e+3

    print("d2F/dmIdm0 (meV)")
    for I in range(0,Nhub):
        wI0 = np.mod(w[I] - w[0], [1.0, 1.0, 1.0])
        for i in range(3):
            if wI0[i] > 0.5:
                wI0[i] -= 1.0
            if np.isclose(wI0[i], 0.5):
                wI0[i] = -0.5
        vI0 = slab.crys2cart(a, wI0)
        print(f"( 1,{I+1:2d}):{JI0[I]*Ry2meV:16.8f}{np.linalg.norm(vI0):16.8f}")
    print('')

    print("d2F/dmi2 (meV)")
    print(f"{JI0[0]*Ry2meV:16.8f}")
    print('')

    sum_JI0_sgni = -np.sum(JI0[1:]*np.array(s[1:Nhub]))
    print("d2F/dmi2 - \u03A3_i {-J_i0 \u03C3_i} (meV)")
    print(f"{(JI0[0] - sum_JI0_sgni)*Ry2meV:16.8f}")
    print('')

    print("d2F/dmi2 (meV)")
    print(f"{d2F_dm2[Nhub,Nhub]*Ry2meV:16.8f}")
    print('')

    print('END')
