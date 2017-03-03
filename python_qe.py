import numpy as np
from forcetheorem import get_eband
from crystal import get_Bi2X3_3QL_card

"""
pwi = read_pwi('data/scf.in')
print(get_pwi_atomic_positions(pwi))
set_pwi_value(pwi, 'CONTROL', 'prefix', "'par'")
write_pwi(pwi)
"""

dir = './data/'
eband_x, ef_x = get_eband(dir+'x.out')
eband_y, ef_y = get_eband(dir+'y.out')
eband_z, ef_z = get_eband(dir+'z.out')

print('Ex-Ez (meV) : %f' %((eband_x - eband_z)*1000))
print('Ey-Ez (meV) : %f' %((eband_y - eband_z)*1000))
print('Ex-Ey (meV) : %f' %((eband_x - eband_y)*1000))

#get_Bi2X3_3QL_card(9.841, 24.27*(np.pi/180.0), 0.4, 0.21, X='Se')