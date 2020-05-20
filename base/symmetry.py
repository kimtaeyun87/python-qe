#-*- coding: utf-8 -*-


import ctypes as ct
import numpy as np
import numpy.ctypeslib as npct


def get_map_IJ_to_I0(w, R, t, i0=0):

    Nw = len(w)
    NR = len(R)


    libsym = ct.CDLL('./libsymmetry.so', mode=ct.RTLD_GLOBAL)

    libsym.get_equiv_atom_map.argtypes = [
        ct.c_int,
        npct.ndpointer(dtype=float, ndim=2, flags='C_CONTIGUOUS'),
        ct.c_int,
        npct.ndpointer(dtype=np.int32, ndim=3, flags='C_CONTIGUOUS'),
        npct.ndpointer(dtype=float, ndim=2, flags='C_CONTIGUOUS')
    ]

    libsym.get_equiv_atom_map.restype = ct.POINTER(ct.POINTER(ct.c_int))

    equiv_atom_map = libsym.get_equiv_atom_map(Nw, w, NR, R, t, i0)

    map_count = npct.as_array(equiv_atom_map[0], (Nw,))
    max_count = max(map_count)

    map_IJ_to_I0 = npct.as_array(equiv_atom_map[1], (Nw,Nw,max_count))

    return map_IJ_to_I0


# python implementation as a reference
def _get_map_IJ_to_I0(Nw, w, NR, R, t, I0=0):
    """
    Nw: number of atoms.
    w: atomic positions in crystal coordinate.
    NR: number of symmetry operations.
    R: rotation part of symmetry operations.
    t: translation part of symmetry operations.
    I0: index of the reference atom.
    """
    map_i0 = []

    for i in range(Nw):
        x = []
        for n in range(NR):
            wi = np.dot(R[n], w[I0]) + t[n]
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

