
import numpy as np

from base.atom import Atom
from base.lattice import Lattice, get_reciprocal


class Crystal:
    def __init__(self, lattice, atomlist):
        self.lattice = lattice
        self.atomlist = atomlist

    @property
    def lattice(self):
        return self._lattice

    @lattice.setter
    def lattice(self, lattice):
        if not isinstance(lattice, Lattice):
            raise TypeError('lattice must be an instance of Lattice.')
        self._lattice = lattice

    @property
    def atomlist(self):
        return self._atomlist

    @atomlist.setter
    def atomlist(self, atomlist):
        if not isinstance(atomlist, list):
            raise TypeError('atomlist must be a \'list\'.')
        if not all(map(lambda atom: isinstance(atom, Atom), atomlist)):
            raise TypeError('atom in atomlist must be a \'instance of Atom\'.')
        self._atomlist = atomlist

    def get_atomlist_position(self):
        atomlist_position = []
        for atom in self.atomlist:
            atomlist_position.append(atom.position)
        return np.array(atomlist_position)
    
    def get_atomlist_name(self):
        atomlist_name = []
        for atom in self.atomlist:
            atomlist_name.append(atom.name)
        return atomlist_name

    def get_atomlist_as_tuplelist(self):
        atomlist_name = self.get_atomlist_name
        atomlist_position = self.get_atomlist_position
        return zip(atomlist_name, atomlist_position)

    def get_crystal_coordinate(self):
        unitcell = self.lattice.unitcell
        reciprocal = get_reciprocal(unitcell)
        position = self.get_atomlist_position()
        return np.dot(position, reciprocal)


def get_crystal_coordinate(position, reciprocal):
    if not isinstance(position, np.ndarray):
        raise TypeError('positions must be an ndarray.')
    if not isinstance(reciprocal, np.ndarray):
        raise TypeError('reciprocal must be an ndarray.')
    if reciprocal.shape != (3,3):
        raise TypeError('reciprocal must be an array of shape (3,3).')
    return np.dot(position, reciprocal)

