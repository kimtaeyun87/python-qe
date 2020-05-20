#-*- coding: utf-8 -*-


import re
import numpy as np


class Atom:
    """Atom defines a simple class for storing atom data.
    name: str
    position: np.ndarray (3,)
    """
    def __init__(self, name, position):
        self.name = name
        self.position = position

    @property
    def name(self, copy=False):
        return self._name

    @name.setter
    def name(self, name):
        if not isinstance(name, str):
            raise TypeError('name of an atom must be a string.')
        self._name = name

    @property
    def position(self):
        return self._position

    @position.setter
    def position(self, position):
        if not isinstance(position, np.ndarray):
            raise TypeError('position of an atom must be a ndarray.')
        if position.shape != (3,):
            raise ValueError('position of an atom must be an array of shape (3,).')
        self._position = position

    @staticmethod
    def guess_symbol_from_name(name):
        if not isinstance(name, str):
            raise TypeError('name of an atom must be a string.')

        symbol = re.split('(\d.*)', name)[0] # pylint: disable=anomalous-backslash-in-string
        if symbol in _atom_symbol_to_number.keys():
            return symbol
        else:
            raise KeyError(f'couldn\'t find symbol for \'{symbol}\' from the table.')

    @staticmethod
    def get_atom_number_from_symbol(symbol):
        if not isinstance(symbol, str):
            raise TypeError('symbol of an atom must be a string.')
        return _atom_symbol_to_number[symbol]

    @staticmethod
    def get_atom_mass_from_symbol(symbol):
        if not isinstance(symbol, str):
            raise TypeError('symbol of an atom must be a string.')
        return _atom_symbol_to_mass[symbol]


_atom_symbol_to_number = {
    'H' : 1,  'He': 2,  'Li': 3,  'Be': 4,  'B' : 5,  'C' : 6,  'N' : 7,  'O' : 8,
    'F' : 9,  'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P' : 15, 'S' : 16,
    'Cl': 17, 'Ar': 18, 'K' : 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V' : 23, 'Cr': 24,
    'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32,
    'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y' : 39, 'Zr': 40,
    'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48,
    'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I' : 53, 'Xe': 54, 'Cs': 55, 'Ba': 56,
    'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64,
    'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72,
    'Ta': 72, 'W' : 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80,
    'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88,
    'Ac': 89, 'Th': 90, 'Pa': 91, 'U' : 92, 'Np': 93, 'Pu': 94
    }

_atom_symbol_to_mass = {
    'H' : 1.008,  'He': 4.0026,  'Li': 6.94,  'Be': 9.0122,  
    'B' : 10.81,  'C' : 12.011,  'N' : 14.007,  'O' : 15.999,
    'F' : 18.998,  'Ne': 20.180, 'Na': 22.990, 'Mg': 24.305, 
    'Al': 13, 'Si': 14, 'P' : 15, 'S' : 16,
    'Cl': 17, 'Ar': 18, 'K' : 19, 'Ca': 20, 
    'Sc': 21, 'Ti': 22, 'V' : 23, 'Cr': 24,
    'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 
    'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32,
    'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 
    'Rb': 37, 'Sr': 38, 'Y' : 39, 'Zr': 40,
    'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 
    'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48,
    'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52,
    'I' : 53, 'Xe': 54, 'Cs': 55, 'Ba': 56,
    'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 
    'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64,
    'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 
    'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72,
    'Ta': 72, 'W' : 74, 'Re': 75, 'Os': 76, 
    'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80,
    'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 
    'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88,
    'Ac': 89, 'Th': 90, 'Pa': 91, 'U' : 92, 
    'Np': 93, 'Pu': 94
    }

