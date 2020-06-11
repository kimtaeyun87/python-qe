#-*- coding: utf-8 -*-


import ctypes as ct
import numpy as np
import numpy.ctypeslib as npct
import os.path


def get_response_matrix(response_vector, equiv_atom_map):
    n = len(response_vector)
    response_matrix = np.zeros((n,n))

    for i in range(n):
        for j in range(n):
            num_map = len(equiv_atom_map[i, j])
            _tmp = 0.0
            for k in range(num_map):
                l = equiv_atom_map[i, j, k]
                _tmp += response_vector[l]
            response_matrix[i, j] = _tmp/num_map
    
    return response_matrix


def get_equiv_atom_map(crystal_coordinate, rotation, translation, index_reference_atom):
    """Wrapper for get_equiv_atom_map defined in equiv_map.c
    """
    libsym = ct.CDLL(os.path.dirname(__file__) + os.path.sep + 'libsymm.so', 
        mode=ct.RTLD_GLOBAL)

    libsym.get_equiv_atom_map.argtypes = [
        ct.c_int,
        npct.ndpointer(dtype=float, ndim=2, flags='C_CONTIGUOUS'),
        ct.c_int,
        npct.ndpointer(dtype=np.int32, ndim=3, flags='C_CONTIGUOUS'),
        npct.ndpointer(dtype=float, ndim=2, flags='C_CONTIGUOUS')
    ]

    libsym.get_equiv_atom_map.restype = ct.POINTER(ct.POINTER(ct.c_int))

    _equiv_atom_map = libsym.get_equiv_atom_map(
        len(crystal_coordinate), crystal_coordinate, len(rotation), rotation, 
        translation, index_reference_atom)

    map_count = npct.as_array(_equiv_atom_map[0], (len(crystal_coordinate),))
    max_count = max(map_count)

    equiv_atom_map = npct.as_array(
        _equiv_atom_map[1], (len(crystal_coordinate),len(crystal_coordinate),max_count))

    return equiv_atom_map


"""
# python implementation as a reference
def _get_equiv_atom_map(crystal_coordinates, rotations, translations, i0=0):
    # Nw: number of atoms.
    # w: atomic positions in crystal coordinate.
    # NR: number of symmetry operations.
    # R: rotation part of symmetry operations.
    # t: translation part of symmetry operations.
    # I0: index of the reference atom.
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
"""

