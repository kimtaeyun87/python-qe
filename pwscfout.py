#-*- coding: utf-8 -*-


import numpy as np
import subprocess


def grep(key, fpath, *args):
    cmd = ["grep", key, fpath] + [x for x in args]
    with subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True) as p:
        stdout = p.stdout.readlines()
    return stdout


def readatomicmoment(fpath, n=0):
    with open(fpath) as f:
        lines = f.readlines()

    m = []
    for i in range(len(lines)):
        if all(x in lines[i] for x in ["atomic", "mx", "my", "mz"]):
            _, _, _, _, _, mx, my, mz = lines[i].split()
            m.append([float(x) for x in [mx, my, mz]])
    return np.array(m[-n:])


def readconstrainedmoment(fpath, n=0):
    with open(fpath) as f:
        lines = f.readlines()

    m = []
    for i in range(len(lines)):
        if all(x in lines[i] for x in ["constrained", "moment"]):
            _, _, _, mx, my, mz = lines[i].split()
            mi = np.array([float(x) for x in [mx, my, mz]])
            if np.isclose(np.linalg.norm(mi), 0.0):
                m.append(np.zeros(3))
            else:
                m.append(mi/np.linalg.norm(mi))
    return np.array(m[-n:])


if __name__ == "__main__":
    name = "test/pwscfout_test.out"

    with open(name) as f:
        lines = f.readlines()

    magmom   = []
    etot     = []
    acc      = []
    hubconst = []
    magconst = []
    iternum  = []
    efermi   = 0.0
    totram   = (0.0, "MB")

    magmom = readatomicmoment(name)

    for i in range(len(lines)):
        #if all(x in lines[i] for x in ["atomic", "mx", "my", "mz"]):
        #    _, _, _, _, _, mx, my, mz = lines[i].split()
        #    magmom.append([float(x) for x in [mx, my, mz]])

        if "total energy" in lines[i]:
            try:
                _, _, _, e, _ = lines[i].replace("!", "").split()
                etot.append(float(e))
            except ValueError:
                if all(x in lines[i] for x in ["sum", "following", "terms"]):
                    pass

        if "accuracy" in lines[i]:
            _, _, _, _, x, _ = lines[i].split()
            acc.append(float(x))

        if all(x in lines[i] for x in ["Hub. E", "const"]):
            _, _, _, _, e, _ = lines[i].split()
            hubconst.append(float(e))

        if "constrained moment" in lines[i]:
            _, _, _, mx, my, mz = lines[i].split()
            magconst.append([float(x) for x in [mx, my, mz]])

        if "Fermi" in lines[i]:
            _, _, _, _, e, _ = lines[i].split()
            efermi = float(e)

        if all(x in lines[i] for x in ["Estimated", "total", "RAM"]):
            _, _, _, _, _, v, u = lines[i].split()
            totram = (float(v), u)
    
    print(magmom)
