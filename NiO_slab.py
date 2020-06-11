if __name__ == '__main__':
    import numpy as np
    import slab

    np.set_printoptions(formatter={
        "float": lambda x: "{:16.8f}".format(x),
        "int": lambda x: "{:4d}".format(x)})

    # NiO conventional cell
    card = slab.readVestaxtl("./NiO.xtl")
    a = card["cell_parameters"]["cell"]
    w = card["atomic_positions"]["position"]
    l = card["atomic_positions"]["label"]
    v = slab.crys2cart(a, w)
    
    # Unitcell
    Xconv = (a, l, v)

    # Convert to surface cell
    gsurf = [1, 1, 1]
    Gsurf = slab.getGvector(a, gsurf)
    asurf = slab.getSurfaceCell(a, gsurf)
    Xsurf = slab.convertCell(*Xconv, asurf)
    asurf, lsurf, vsurf = Xsurf

    # Coordinate transformation
    Xsurf = slab.transformCoordinate(*Xsurf, slab.standardizeCell(asurf))
    asurf, lsurf, vsurf = Xsurf

    # Make slab
    asurf_super = asurf.copy()
    asurf_super[2] *= 2.0
    Xsurf_super = slab.convertCell(*Xsurf, asurf_super)
    _, lsurf_super, vsurf_super = Xsurf_super

    with open("tmp.xsf", 'w') as f:
        slab.writexsf(f, *Xsurf_super, mode='crystal')
    import subprocess
    subprocess.call(["/home/kimtaeyun/xcrysden/1.6.2/xcrysden", "--xsf", "tmp.xsf"])
    subprocess.call(["rm", "-f", "/tmp/tmp.xsf"])

    ii = np.argsort(lsurf_super)
    print(asurf_super)
    print(np.sort(lsurf_super))
    print(vsurf_super[ii])

