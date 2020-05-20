#-*- coding: utf-8 -*-


import numpy as np


class Lattice:
    """Lattice defines a class for storing lattice data.
    unitcell: np.ndarray (3,3)
    """
    def __init__(self, unitcell):
        self.unitcell = unitcell

    @property
    def unitcell(self):
        return self._unitcell

    @unitcell.setter
    def unitcell(self, unitcell):
        if not isinstance(unitcell, np.ndarray):
            raise TypeError('unitcell must be an ndarray.')
        if unitcell.shape != (3,3):
            raise TypeError('unitcell must be an anrray of shape (3,3).')
        if np.linalg.det(unitcell) < 0.0:
            raise ValueError('det(unicell) is negative.')
        if np.isclose(np.linalg.det(unitcell), 0.0):
            raise ValueError('det(unicell) is zero.')
        self._unitcell = unitcell

    @staticmethod
    def get_reciprocal(unitcell):
        if not isinstance(unitcell, np.ndarray):
            raise TypeError('unitcell must be an ndarray.')
        if unitcell.shape != (3,3):
            raise TypeError('unitcell must be an array of shape (3,3).')

        Volume = np.linalg.det(unitcell)
        if Volume < 0:
            raise ValueError("det(unitcell) is negative.")
        elif np.isclose(Volume, 0.0):
            raise ValueError("det(unitcell) is zero.")

        reciprocal = np.zeros((3, 3))
        reciprocal[:, 0] = np.cross(unitcell[1], unitcell[2])
        reciprocal[:, 1] = np.cross(unitcell[2], unitcell[0])
        reciprocal[:, 2] = np.cross(unitcell[0], unitcell[1])

        return reciprocal/Volume

    @staticmethod
    def get_unitcell_from_cell_params(cell_params):
        if not isinstance(cell_params, np.ndarray):
            raise TypeError('cell_params must be an ndarray.')
        if cell_params.shape != (6,):
            raise TypeError('cell_params must be an array of shape (6,).')

        A, B, C = cell_params[:3]
        sinα, sinβ, sinγ = np.sin(np.deg2rad(cell_params[3:])) # pylint: disable=unused-variable
        cosα, cosβ, cosγ = np.cos(np.deg2rad(cell_params[3:]))
        
        unitcell = np.zeros((3, 3))
        
        unitcell[0, :] = [A, 0.0, 0.0]
        unitcell[1, :] = [B*cosγ, B*sinγ, 0.0]
        unitcell[2, :] = [
            C*cosβ,
            C*(cosα - cosβ*cosγ)/sinγ,
            C*np.sqrt(1.0 + 2.0*cosα*cosβ*cosγ - cosα**2 - cosβ**2 - cosγ**2)/sinγ
        ]

        return unitcell

    @staticmethod
    def get_cell_params_from_unitcell(unitcell):
        if not isinstance(unitcell, np.ndarray):
            raise TypeError('unitcell must be an ndarray.')
        if unitcell.shape != (3, 3):
            raise TypeError("unitcell must be an array of shape (3,3).")

        A, B, C = [np.linalg.norm(unitcell[i]) for i in range(3)]
        
        cosα = np.dot(unitcell[1]/B, unitcell[2]/C)
        cosβ = np.dot(unitcell[0]/A, unitcell[2]/C)
        cosγ = np.dot(unitcell[0]/A, unitcell[1]/B)
        
        α = np.rad2deg(np.arccos(cosα))
        β = np.rad2deg(np.arccos(cosβ))
        γ = np.rad2deg(np.arccos(cosγ))
        
        return np.array([A, B, C, α, β, γ])

    @staticmethod
    def get_standard_unitcell(unitcell):
        cell_params = Lattice.get_cell_params_from_unitcell(unitcell)
        return Lattice.get_unitcell_from_cell_params(cell_params)


