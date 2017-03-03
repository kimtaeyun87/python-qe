def read_pwi(fn):
    with open(fn, 'r') as f:
        fl = f.readlines()
    pwi = {}
    card = []
    for i, ln in enumerate(fl):
        nf = len(ln.split())
        if (nf == 0):
            continue
        else:
            if (ln.split()[0][0] == '!'): continue
        if (ln.split()[0][0] == '&'):
            section_name = ln.split()[0][1:]
            section_read = True
            card_read = False
            pwi[section_name] = {}
        elif (ln.split()[0][0] == '/'):
            section_read = False
            card_read = True
        else:
            if (section_read and card_read): print('line %d: section or card?' %i)
            if (section_read and (not card_read)):
                if (nf != 3): print('lind %d: number of fields not 3' %i)
                key, sgn, val = ln.split()
                if (sgn != '='): print('line %d: = might be omitted' %i)
                pwi[section_name][key] = val
            if ((not section_read) and card_read):
                card.append(ln.split())
            if ((not section_read) and (not card_read)):
                continue
    ntyp = int(pwi['SYSTEM']['ntyp'])
    nat = int(pwi['SYSTEM']['nat'])
    for i, ln in enumerate(card):
        if (ln[0] == 'ATOMIC_POSITIONS'):
            pwi['ATOMIC_POSITIONS'] = {}
            pwi['ATOMIC_POSITIONS']['option'] = ln[1]
            pwi['ATOMIC_POSITIONS']['value'] = card[i+1:i+1+nat]
        elif (ln[0] == 'ATOMIC_SPECIES'):
            pwi['ATOMIC_SPECIES'] = {}
            pwi['ATOMIC_SPECIES']['option'] = ''
            pwi['ATOMIC_SPECIES']['value'] = card[i+1:i+1+ntyp]
        elif (ln[0] == 'K_POINTS'):
            pwi['K_POINTS'] = {}
            pwi['K_POINTS']['option'] = ln[1]
            if (ln[1] == 'automatic'):
                pwi['K_POINTS']['value'] = card[i+1]
            elif (ln[1] == 'gamma'):
                continue
            else:
                nks = int(card[i+1])
                pwi['K_POINTS']['nks'] = card[i+1]
                pwi['K_POINTS']['value'] = card[i+2:i+2+nks]
    return pwi

def write_pwi(pwi):
    section_list = ['CONTROL', 'SYSTEM', 'ELECTRONS', 'IONS']
    card_list = ['ATOMIC_SPECIES', 'ATOMIC_POSITIONS', 'K_POINTS']
    for section_name in section_list:
        if (section_name in pwi):
            print('&'+section_name)
            for key in pwi[section_name]:
                print(key,'=',pwi[section_name][key])
            print('/')
        else:
            continue
    for card_name in card_list:
        if (card_name in pwi):
            print(card_name, pwi[card_name]['option'])
            if (card_name == 'ATOMIC_SPECIES'):
                for atyp in pwi[card_name]['value']:
                    print('{} {} {}'.format(atyp[0],atyp[1],atyp[2]))
            elif (card_name == 'ATOMIC_POSITIONS'):
                for atom in pwi[card_name]['value']:
                    print('{} {} {} {}'.format(atom[0],atom[1],atom[2],atom[3]))
            elif (card_name == 'K_POINTS'):
                if (pwi[card_name]['option'] == 'automatic'):
                    grid = pwi[card_name]['value']
                    print('{} {} {} {} {} {}'.format(grid[0],grid[1],grid[2],grid[3],grid[4],grid[5]))
                elif (pwi[card_name]['option'] == 'gamma'):
                    continue
                else:
                    print(pwi[card_name]['nks'])
                    for kpt in pwi[card_name]['value']:
                        print('{} {} {} {}'.format(kpt[0],kpt[1],kpt[2],kpt[3]))
        else:
            continue

def set_pwi_value(pwi, section_name, key, val):
    if ((section_name in pwi) and (key in pwi[section_name])):
        pwi[section_name][key] = val

def get_pwi_value(pwi, section_name, key):
    if ((section_name in pwi) and (key in pwi[section_name])):
        return pwi[section_name][key]

def get_pwi_atomic_positions(pwi):
    nat = int(get_pwi_value(pwi, 'SYSTEM', 'nat'))
    atom = get_pwi_value(pwi, 'ATOMIC_POSITIONS', 'value')
    position = np.zeros((nat, 3))
    for i in range(nat):
        position[i] = np.asarray(atom[i][1:])
    return position