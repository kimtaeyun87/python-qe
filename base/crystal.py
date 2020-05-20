#-*- coding: utf-8 -*-

import numpy as np

from base.atom import Atom
from base.lattice import Lattice


class Crystal:
    def __init__(self, lattice, atoms):
        self.lattice = lattice
        self.atoms = atoms

    @property
    def lattice(self):
        return self._lattice

    @lattice.setter
    def lattice(self, lattice):
        if not isinstance(lattice, Lattice):
            raise TypeError('lattice must be an instance of Lattice.')
        self._lattice = lattice

    @property
    def atoms(self):
        return self._atoms

    @atoms.setter
    def atoms(self, atoms):
        if not isinstance(atoms, list):
            raise TypeError('atoms must be a \'list\' of instances of Atom.')
        if not all(map(lambda atom: isinstance(atom, Atom), atoms)):
            raise TypeError('atoms must be a list of \'instances of Atom\'.')
        self._atoms = atoms

    def get_atom_names(self):
        names = []
        for atom in self.atoms:
            names.append(atom.name)
        return np.array(names)

    def get_atom_positions(self):
        positions = []
        for atom in self.atoms:
            positions.append(atom.position)
        return np.array(positions)

    @staticmethod
    def get_crystal_coordinates(positions, reciprocal):
        if not isinstance(positions, np.ndarray):
            raise TypeError('positions must be an ndarray.')
        if not isinstance(reciprocal, np.ndarray):
            raise TypeError('reciprocal must be an ndarray.')
        if reciprocal.shape != (3,3):
            raise TypeError('reciprocal must be an array of shape (3,3).')
        return np.dot(positions, reciprocal)

