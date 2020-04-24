#-*- coding: utf-8 -*-


import f90nml


def readcard(lines, i, n):
    try:
        tag, opt = lines[i].lower().split()
    except ValueError:
        if len(lines[i].split()) == 1:
            tag = lines[i].lower().rstrip()
            opt = ""
    data = lines[i+1:i+1+n]
    return tag, opt, data


def checkload(data):
    if isinstance(data, list):
        if len(data) == 0:
            return False
        else:
            return True
    else:
        raise TypeError


def readpwscfin(fpath):
    nml = f90nml.read(fpath)

    numatom = int(nml["system"]["nat"])
    numtype = int(nml["system"]["ntyp"])

    with open(fpath) as f:
        lines = f.readlines()

    atomtype   = []
    atommass   = []
    atompseudo = []

    atomlist   = []
    atompos    = []
    atomposopt = ""

    cell       = []
    cellopt    = ""

    extforce   = []

    kptgrid    = []
    kptlist    = []
    kptopt     = ""
    numkpt     = 0

    for i in range(len(lines)):
        cardhead = lines[i].lower().rstrip()

        if "atomic_species" in cardhead:
            tag, opt, data = readcard(lines, i, numtype)

            for j in range(numtype):
                atomtype.append(data[j].split()[0])
                atommass.append(float(data[j].split()[1]))
                atompseudo.append(data[j].split()[2])

        elif "atomic_positions" in cardhead:
            tag, opt, data = readcard(lines, i, numatom)

            atomposopt = opt

            for j in range(numatom):
                atomlist.append(data[j].split()[0])
                atompos.append([float(x) for x in data[j].split()[1:]])

        elif "cell_parameters" in cardhead:
            tag, opt, data = readcard(lines, i, 3)

            cellopt = opt

            for j in range(3):
                cell.append([float(x) for x in data[j].split()])

        elif "k_points" in cardhead:
            tag, opt, data = readcard(lines, i, 0)

            kptopt = opt

            if opt == "gamma":
                continue

            elif opt == "automatic":
                _, _, data = readcard(lines, i, 1)

                kptgrid += [int(x) for x in data[0].split()]

            else:
                _, _, data = readcard(lines, i, 1)

                numkpt = int(data[0].split()[0])

                _, _, data = readcard(lines, i, numkpt+1)

                for j in range(1, numkpt+1):
                    kptlist.append([float(x) for x in data[j].split()])

        elif "atomic_forces" in cardhead:
            tag, opt, data = readcard(lines, i, numatom)

            for j in range(numatom):
                extforce.append([float(x) for x in data[j].split()[1:]])

        else:
            continue

    # sanity check
    # control, system and electrons 
    # atomic_species, atomic_positions and k_points

    card = {}
    card["atomic_species"]   = {"type": atomtype, "mass": atommass, "pseudo": atompseudo}
    card["atomic_positions"] = {"opt": atomposopt, "list": atomlist, "pos": atompos}
    card["k_points"]         = {"opt": kptopt, "nkpt": numkpt, "grid": kptgrid, "list": kptlist}
    card["cell_parameters"]  = {"opt": cellopt, "cell": cell}
    card["atomic_forces"]    = {"list": atomlist, "extforce": extforce}

    return nml, card


def printpwscfin(nml, card):
    print(nml)

    headform = "{:s} {:s}"

    numatom = int(nml["system"]["nat"])
    numtype = int(nml["system"]["ntyp"])

    tag = "atomic_species"
    print(headform.format(tag, "").rstrip())
    for i in range(numtype):
        x      = card[tag]["type"][i]
        mass   = card[tag]["mass"][i]
        pseudo = card[tag]["pseudo"][i]
        print(f"    {x:4s}    {mass:8f}    {pseudo:80s}".rstrip())

    tag = "atomic_positions"
    opt = card[tag]["opt"]
    print(headform.format(tag, opt).rstrip())
    for i in range(numatom):
        x   = card[tag]["list"][i]
        pos = card[tag]["pos"][i]
        print(f"    {x:4s}" + "".join(f"    {y:8f}" for y in pos).rstrip())

    tag = "k_points"
    opt = card[tag]["opt"]
    print(headform.format(tag, opt).rstrip())
    if opt == "gamma":
        pass
    elif opt == "automatic":
        kptgrid = card[tag]["grid"]
        print("".join(f"{x:4d}" for x in kptgrid).rstrip())
    else:
        numkpt = card[tag]["nkpt"]
        print(f"{numkpt:4d}")
        for i in range(numkpt):
            kpt = card[tag]["list"][i]
            print("".join(f"    {x:8f}" for x in kpt).rstrip())


if __name__ == "__main__":
    # test
    name = "test/pwscfin_test.in"
    nml, card = readpwscfin(name)
    printpwscfin(nml, card)
