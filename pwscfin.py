#-*- coding: utf-8 -*-


import f90nml
import numpy as np


def readcard(lines, i, n):
    try:
        tag, opt = lines[i].lower().split()
    except ValueError:
        if len(lines[i].split()) == 1:
            tag = lines[i].lower().rstrip()
            opt = ""

    data = lines[i+1:i+1+n]

    return tag, opt, data


def readpwscfin(fpath):
    nml = f90nml.read(fpath)

    numatom = int(nml["system"]["nat"])
    numtype = int(nml["system"]["ntyp"])

    with open(fpath) as f:
        lines = f.readlines()

    typelabel  = []
    typemass   = []
    typepseudo = []

    atomlabel  = []
    atompos    = np.zeros((numatom, 3))
    atomposopt = ""

    cell    = np.zeros((3, 3))
    cellopt = ""

    atomforce   = np.zeros((numatom, 3))

    kptgrid    = []
    kptlist    = []
    kptopt     = ""
    numkpt     = 0

    for i in range(len(lines)):
        head = lines[i].lower().rstrip()

        if "atomic_species" in head:
            tag, opt, data = readcard(lines, i, numtype)

            for j in range(numtype):
                typelabel.append(data[j].split()[0])
                typemass.append(float(data[j].split()[1]))
                typepseudo.append(data[j].split()[2])

        elif "atomic_positions" in head:
            tag, opt, data = readcard(lines, i, numatom)
            atomposopt = opt

            for j in range(numatom):
                atomlabel.append(data[j].split()[0])
                atompos[j, :] = np.array([float(x) for x in data[j].split()[1:]])

        elif "cell_parameters" in head:
            tag, opt, data = readcard(lines, i, 3)
            cellopt = opt

            for j in range(3):
                cell[j, :] = np.array([float(x) for x in data[j].split()])

        elif "k_points" in head:
            tag, opt, data = readcard(lines, i, 0)
            kptopt = opt

            if opt == "gamma":
                continue

            elif opt == "automatic":
                _, _, data = readcard(lines, i, 1)

                kptgrid = [int(x) for x in data[0].split()]

            else:
                _, _, data = readcard(lines, i, 1)

                numkpt = int(data[0].split()[0])

                _, _, data = readcard(lines, i, numkpt+1)

                for j in range(1, numkpt+1):
                    kptlist.append(np.array([float(x) for x in data[j].split()]))

        elif "atomic_forces" in head:
            tag, opt, data = readcard(lines, i, numatom)

            for j in range(numatom):
                atomforce[j, :] = np.array([float(x) for x in data[j].split()[1:]])

        else:
            continue

    # sanity check
    # control, system and electrons 
    # atomic_species, atomic_positions and k_points

    card = {}
    card["atomic_species"]   = {"label": typelabel, "mass": typemass, "pseudo": typepseudo}
    card["atomic_positions"] = {"option": atomposopt, "label": atomlabel, "position": atompos}
    card["k_points"]         = {"option": kptopt, "number": numkpt, "grid": kptgrid, "list": kptlist}
    card["cell_parameters"]  = {"option": cellopt, "cell": cell}
    card["atomic_forces"]    = {"label": atomlabel, "force": atomforce}

    return nml, card


def printpwscfin(nml, card):
    print(nml)

    headform = "\n{:s} {:s}"

    numatom = int(nml["system"]["nat"])
    numtype = int(nml["system"]["ntyp"])

    tag = "atomic_species"
    print(headform.format(tag, "").rstrip())
    for i in range(numtype):
        x      = card[tag]["label"][i]
        mass   = card[tag]["mass"][i]
        pseudo = card[tag]["pseudo"][i]
        print(f"{x:4s}{mass:8.3f}{'':4s}{pseudo:80s}".rstrip())

    tag = "atomic_positions"
    opt = card[tag]["option"]
    print(headform.format(tag, opt).rstrip())
    for i in range(numatom):
        x   = card[tag]["label"][i]
        pos = card[tag]["position"][i]
        print(f"{x:4s}" + "".join(f"{y:16.8f}" for y in pos).rstrip())

    tag = "k_points"
    opt = card[tag]["option"]
    print(headform.format(tag, opt).rstrip())
    if opt == "gamma":
        pass
    elif opt == "automatic":
        kptgrid = card[tag]["grid"]
        print("".join(f"{x:4d}" for x in kptgrid).rstrip())
    else:
        numkpt = card[tag]["number"]
        print(f"{numkpt:4d}")
        for i in range(numkpt):
            kpt = card[tag]["list"][i]
            print("".join(f"{x:16.8f}" for x in kpt).rstrip())

    tag = "cell_parameters"
    opt = card[tag]["option"]
    print(headform.format(tag, opt).rstrip())
    for i in range(3):
        a = card[tag]["cell"][i]
        print(f"{a[0]:16.8f}{a[1]:16.8f}{a[2]:16.8f}".rstrip())


if __name__ == "__main__":
    name = "test/pwscfin_test.in"
    nml, card = readpwscfin(name)
    printpwscfin(nml, card)
