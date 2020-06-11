from base.atom import _atom_symbol_to_number


def writexsf(fpath, a, l, v, mode="crystal"):
    """ Write a xsf (xcrysden) file
    Args:
        fpath: file
        a: unitcell
        l: name
        v: positions
        mode: "crystal", "slab", "polymer", "dot". Default to "crystal"
    """
    with open(fpath, 'w') as f:
        f.write(f"{mode.upper():s}\n")
        f.write("PRIMVEC\n")

        for i in range(3):
            x, y, z = a[i, :]
            f.write(f"{x:16.8f}{y:16.8f}{z:16.8f}\n")

        f.write("CONVVEC\n")

        for i in range(3):
            x, y, z = a[i, :]
            f.write(f"{x:16.8f}{y:16.8f}{z:16.8f}\n")

        f.write("PRIMCOORD\n")
        f.write(f"{len(v):4d} 1\n")

        for i in range(len(v)):
            x, y, z = v[i, :]
            f.write(f"{_atom_symbol_to_number[l[i]]:4d}{x:16.8f}{y:16.8f}{z:16.8f}\n")

