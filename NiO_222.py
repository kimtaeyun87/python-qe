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

a2, l2, v2 = slab.convertCell(a, l, v, a*2)
w2 = slab.cart2crys(a2, v2)

print(a2)
print(l2)
print(v2)
print(w2)
