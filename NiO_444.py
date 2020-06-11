import numpy as np

def get_derivative_magnetic_moment(parameter_list, magnetic_moment_list):
    """
    Args:
        parameter_list: array-like parameters.
        magnetic_moment_list: list of magnetic moment arrays.
    Returns:
        derivative
    """
    natom, ndim = magnetic_moment_list[0].shape
    _tmp = np.zeros((len(p), natom*ndim))
    for i in range(len(parameter_list)):
        _tmp[i] = magnetic_moment_list[i].flatten()
    derivative = np.polyfit(parameter_list, _tmp, deg=len(p)-1)[-2, :]
    return np.reshape(derivative, (natom, ndim))


class Interpolator:
    def __init__(self, parameter, vector):
        self.parameter = parameter
        self.vector = vector

    def interpolate(self):
        assert(len(self.parameter) == len(self.vector))
        # check dimension of state vector
        n, m = self.vector.shape
        _tmp = np.zeros((len(self.parameter), n*m))


if __name__ == "__main__":
    import numpy as np
    import slab

    import base.fmt
    import base.crystal
    import base.pwi
    import base.pwo
    import symm.equivatom
    import os

    # base directory
    outdir = 'test_NiO/perturb_Ni1_2' + os.path.sep

    # load pwin data
    pwi = base.pwi.Pwin(outdir + 'scf.in')

    # hard-coded perturbation parameters
    p = [0.0, 0.01, -0.01, 0.005, -0.005]
    if len(p) < 2:
        raise ValueError("lenth of p must be >= 2")

    # scf calculations with transverse and longitudinal perturbations
    m_lx_list = ['scf.out', 'dm_dlx_1.out', 'dm_dlx_2.out', 'dm_dlx_3.out', 'dm_dlx_4.out']
    m_lz_list = ['scf.out', 'dm_dlz_1.out', 'dm_dlz_2.out', 'dm_dlz_3.out', 'dm_dlz_4.out']



    # derivative of magnetic moment with respect to perturbations
    dm_dlx = get_derivative_magnetic_moment(
        p, [base.pwo.Pwout(f'{outdir}/{x}').iterations[-1].magnetic_moment for x in m_lx_list])

    dm_dlz = get_derivative_magnetic_moment(
        p, [base.pwo.Pwout(f'{outdir}/{x}').iterations[-1].magnetic_moment for x in m_lz_list])

    # build parts of the response matrix
    equiv_atom_map = pwi.get_symmetrizing_map(index_reference_atom=0)
    dmx_dlx = symm.equivatom.get_response_matrix(dm_dlx[:, 0], equiv_atom_map)
    dmz_dlx = symm.equivatom.get_response_matrix(dm_dlx[:, 2], equiv_atom_map)
    dmx_dlz = symm.equivatom.get_response_matrix(dm_dlz[:, 0], equiv_atom_map)
    dmz_dlz = symm.equivatom.get_response_matrix(dm_dlz[:, 2], equiv_atom_map)

    # build the entire response matrix
    num_hub = pwi.get_num_hub()
    dm_dl = np.zeros((num_hub*2, num_hub*2))
    dm_dl[:num_hub, :num_hub] = dmx_dlx
    dm_dl[num_hub:, :num_hub] = dmz_dlx
    dm_dl[:num_hub, num_hub:] = dmx_dlz
    dm_dl[num_hub:, num_hub:] = dmz_dlz

    # build the second derivative matrix of the energy.
    dl_dm = np.linalg.inv(dm_dl)
    Lambda = 0.1
    d2F_dm2 = -dl_dm - Lambda*np.eye(num_hub*2)

    # unitcell and crystal coordinates for calculating distances
    a, w, _, m = pwi.get_spglib_cell()

    # starting magnetization for further information
    JI0 = d2F_dm2[:num_hub, :num_hub][0]
    Ry2meV = 13.6*1e+3

    #
    print("\nd2F/dmIdm0 (meV) distance (Å)")
    for I in range(0, num_hub):
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

    sum_JI0_sgni = -np.sum(JI0[1:]*np.array(m[1:num_hub]))
    print("d2F/dmi2 - \u03A3_i {-J_i0 \u03C3_i} (meV)")
    print(f"{(JI0[0] - sum_JI0_sgni)*Ry2meV:16.8f}")
    print('')

    print("d2F/dmi2 (meV)")
    print(f"{d2F_dm2[num_hub,num_hub]*Ry2meV:16.8f}")
    print('')

