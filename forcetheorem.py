def convert_ns_noncolin(nat, ldim, dir=None):
    if (dir == None):
        fn = './occup.txt'
    else:
        fn = dir+'/occup.txt'
    with open(fn, 'r') as f:
        fl = f.readlines()
    ns = []
    for ln in fl:
        for wd in ln.split():
            ns.append(complex(wd))
    ntot = len(ns)
    if (ntot != nat*2*ldim*ldim): return
    ns_nc = np.zeros(ntot*2, dtype='complex')
    for n in range(nat):
        for spin in range(2):
            for m2 in range(ldim):
                for m1 in range(ldim):
                    idx = m1 + m2*ldim + spin*ldim*ldim + n*2*ldim*ldim
                    idx_nc = m1 + m2*ldim + spin*ldim*ldim + n*4*ldim*ldim
                    ns_nc[idx_nc] = ns[idx]
    for i in range(ntot*2):
        print('(%E,%E)' %(ns_nc[i].real, ns_nc[i].imag))

def get_eband_tot(fn):
    with open(fn, 'r') as f:
        fl = f.readlines()
    for ln in fl:
        if (len(ln.split()) == 0): continue
        if (ln.split()[0] == 'eband_tot,'):
            eband_tot = float(ln.split()[4])
            eband_proj_tot = float(ln.split()[5])
    return eband_tot, eband_proj_tot

def get_eband(fn):
    with open(fn, 'r') as f:
        fl = f.readlines()
    for ln in fl:
        if (len(ln.split()) == 0): continue
        if (ln.split()[0] == 'eband,'):
            eband = float(ln.split()[4])
            Ef = float(ln.split()[5])
    return eband, Ef
