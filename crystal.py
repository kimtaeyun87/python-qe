import numpy as np

def rhombo(a, alpha):
    c = np.cos(alpha)
    tx = np.sqrt((1.0-c)*0.5)
    ty = np.sqrt((1.0-c)/6.0)
    tz = np.sqrt((1.0+2.0*c)/3.0)
    cell = np.array([[a*tx, -a*ty, a*tz], [0, 2.0*a*ty, a*tz], [-a*tx, -a*ty, a*tz]], dtype = 'float')
    return cell 

def rhombo_to_hexa(cell_r):
    T = np.array([[1, 0, -1], [-1, 1, 0], [1, 1, 1]], dtype = 'float')
    cell_h = np.dot(T, cell_r)
    return cell_h

def reciprocal(cell):
    recp = np.array([np.cross(cell[1],cell[2]), np.cross(cell[2],cell[0]), np.cross(cell[0],cell[1])], dtype = 'float')
    return recp/np.linalg.det(cell)

def get_supercell_card(atom_name, vec, cell, supercell):
    Nmax = 3; Nmin = -3
    d = np.dot(reciprocal(cell), vec)
    card = []
    for n1 in range(Nmin, Nmax):
        for n2 in range(Nmin, Nmax):
            for n3 in range(Nmin, Nmax):
                R = np.array([n1, n2, n3], dtype = 'float')
                vecR = np.dot(d+R, cell)
                dR = np.dot(reciprocal(supercell), vecR)
                eps = dR - np.around(dR)
                for i in range(3):
                    if (np.abs(eps[i]) < 1.0e-3): dR[i] -= eps[i]
                if (((dR[0]>=0) and (dR[0]<1)) and ((dR[1]>=0) and (dR[1]<1)) and ((dR[2]>=-0.5) and (dR[2]<0.5))):
                    vec_supercell = np.dot(dR, supercell)
                    card.append([atom_name, vec_supercell[0], vec_supercell[1], vec_supercell[2]])
    return card

def get_Bi2X3_3QL_card(a, alpha, u, v, X):
    card = []
    rhom = rhombo(a, alpha)
    hexa = rhombo_to_hexa(rhom)
    for atom_name, r in [(X, 0), (X, v), (X, -v), ('Bi', u), ('Bi', -u)]:
        d = np.array([r, r, r], dtype = 'float')
        vec = np.dot(d, rhom)
        card += get_supercell_card(atom_name, vec, rhom, hexa)
    card.sort(key = lambda ln: ln[3])
    print('ATOMIC_POSITIONS angstrom')
    for ln in card:
        print('%4s %f %f %f' %(ln[0], ln[1], ln[2], ln[3]))