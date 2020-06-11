# to calculate onsite exchange parameter
# we need at least two scf calculations
reference = 'test_NiO/perturb_Ni1_2/scf.out'
perturbed = 'test_NiO/perturb_Ni1_2/dm_dlz_1.out'

# hard-coded perturbation parameters
lambdas = [0.0, 0.01]

# get magnetic moments of the last iteration of each pwout file.
import base.pwo
m_ref = base.pwo.Pwout(reference).iterations[-1].magnetic_moment
m0 = base.pwo.Pwout(perturbed).iterations[0].magnetic_moment
m = base.pwo.Pwout(perturbed).iterations[-1].magnetic_moment

# calculate response vector delta_m
delta_m = (m - m_ref)/(lambdas[1] - lambdas[0])
delta_m0 = (m0 - m_ref)/(lambdas[1] - lambdas[0])

# find symmetry operations
import base.pwi
import symm.equivatom
import numpy as np

pwi = base.pwi.Pwin('test_NiO/perturb_Ni1_2/scf.in')
equiv_atom_map = pwi.get_symmetrizing_map(0)

# calculate response matrix delta_m_delta_l
delta_m_delta_l = symm.equivatom.get_response_matrix(delta_m[:, 2], equiv_atom_map)
delta_m0_delta_l = symm.equivatom.get_response_matrix(delta_m0[:, 2], equiv_atom_map)

delta_l_delta_m = np.linalg.inv(delta_m_delta_l)
delta_l_delta_m0 = np.linalg.inv(delta_m0_delta_l)

onsite_exchange = - delta_l_delta_m + delta_l_delta_m0
print('J = {:f} meV'.format(onsite_exchange[0, 0]*13.6*1000.0))

