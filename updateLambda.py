import subprocess as subproc
import numpy as np


def grep(key, filename, args=""):
    cmd = "grep '{:s}' {:s} {:s}".format(key, filename, args) 
    return subproc.check_output(cmd, shell=True, universal_newlines=True)


def readAtomicForces(filename, nat):
    with open(str(filename), 'r') as f:
        flinesraw = f.readlines()
        for i, line in enumerate(flinesraw):
            if "ATOMIC_FORCES" in line:
                lambdaraw = flinesraw[i+1: i+1+nat]
                break
    atoms = []
    l = []
    for line in lambdaraw:
        atoms.append(line.split()[0])
        l.append(list(map(float, line.split()[1:4])))
    l = np.array(l)
    return atoms, l


def readAtomicMoments(filename, nathub):
    m = []
    for line in grep("atomic mx", filename, "| tail -{}".format(nathub)).splitlines():
        m.append(list(map(float, line.split()[5:8])))
    return np.array(m, dtype=float)


def readConstrainedMoments(filename, nat, nathub):
    constmomentsraw = grep("constrained moment", filename, args="| tail -{} | head -{}".format(nat, nathub))
    e = []
    for i, line in enumerate(constmomentsraw.splitlines()):
        e.append(list(map(float, line.split()[3:6])))
    e = np.array(e)
    for i in range(nathub):
        e[i] /= np.sqrt(np.sum(e[i]**2))
    return e


def calcLambdaParaPerp(l, e, nathub):
    lpara = []
    lperp = []
    for i in range(nathub):
        lpara.append(np.dot(l[i], e[i]))
        lperp.append(l[i] - lpara[i]*e[i])
    lpara = np.array(lpara)
    lperp = np.array(lperp)
    return lpara, lperp


def calcNextLambda(m, e, l, nathub):
    lpara, lperp = calcLambdaParaPerp(l, e, nathub)
    for i in range(nathub):
        mnorm = np.sqrt(np.sum(m[i]**2))
        mperp = m[i] - np.dot(m[i], e[i])*e[i]
        l[i] += -lpara[i]*mperp/mnorm
    return l


def printLambda(atoms, l, nat):
    print("\nATOMIC_FORCES")
    for i in range(nat):
        print("{:8s} {} {} {}".format(atoms[i], l[i, 0], l[i, 1], l[i, 2]))


if __name__ == "__main__":
    # example script
    scfout = "scf.out"
    scfin = "scf.nml"

    nathub = 32
    nat = int(grep("number of atoms/cell", scfout).split()[-1])

    m = readAtomicMoments(scfout, nathub)
    e = readConstrainedMoments(scfout, nat, nathub)
    atoms, l = readAtomicForces(scfin, nat)
    l = calcNextLambda(m, e, l, nathub)
    printLambda(atoms, l, nat)

