#-*- coding: utf-8 -*-


import f90nml
from collections import namedtuple
import numpy as np

from base.atom import Atom
from base.lattice import Lattice
from base.crystal import Crystal


AtomicSpecies = namedtuple('AtomicSpecies', 'name, mass, pseudo')
AtomicPositions = namedtuple('AtomicPositions', 'option, name, position')
AtomicForces = namedtuple('AtomicForces', 'name, force')
CellParameters = namedtuple('CellParameters', 'option, unitcell')


class PwInput:
    def __init__(self, fpath):
        nml, card = _read_from(fpath)

        Control = namedtuple('Control', nml['control'].keys())
        System = namedtuple('System', nml['system'].keys())
        Electrons = namedtuple('Electrons', nml['electrons'].keys())

        # control and system are required
        self.control = Control(**nml['control'])
        self.system = System(**nml['system'])

        try:
            self.electrons = Electrons(**nml['electrons'])
        except:
            self.electrons = None

        # atomic_species and atomic_positions are required
        self.atomic_species = AtomicSpecies(**card['atomic_species'])
        self.atomic_positions = AtomicPositions(**card['atomic_positions'])

        try:
            self.atomic_forces = AtomicForces(**card['atomic_forces'])
        except KeyError:
            self.atomic_forces = None

        try:
            self.cell_parameters = CellParameters(**card['cell_parameters'])
        except KeyError:
            self.cell_parameters = None

    def get_crystal(self):
        unitcell = self.cell_parameters.unitcell
        names = self.atomic_positions.name
        positions = self.atomic_positions.position

        lattice = Lattice(np.array(unitcell))

        atoms = []
        for name, position in zip(names, positions):
            atoms.append(Atom(name, np.array(position)))
        
        return Crystal(lattice, atoms)

    def get_ishubbard(self):
        ntyp = self.system.ntyp
        names = self.atomic_species.name
        hubbard_u = self.system.hubbard_u

        if ntyp < len(hubbard_u):
            raise ValueError('ntyp must be greater than the number of hubbard atoms.')
        else:
            ishubbard = {}
            for i, name in enumerate(names):
                try:
                    u = hubbard_u[i]
                except IndexError:
                    u = None
                if not u:
                    ishubbard[name] = False
                else:
                    ishubbard[name] = True

        return ishubbard

    def get_nhubbard(self):
        ishubbard = self.get_ishubbard()
        names = self.atomic_positions.name

        nhubbard = 0
        for name in names:
            if ishubbard[name]:
                nhubbard += 1

        return nhubbard

    def get_starting_magnetization(self):
        ntyp = self.system.ntyp
        names = self.atomic_species.name
        _starting_magnetization = self.system.starting_magnetization

        if ntyp < len(_starting_magnetization):
            raise ValueError('ntyp must be greater than the number of hubbard atoms.')
        else:
            starting_magnetization = {}
            for i, name in enumerate(names):
                try:
                    m = _starting_magnetization[i]
                except IndexError:
                    m = 0.0

                starting_magnetization[name] = m

        return starting_magnetization

    def get_spglib_cell(self, magnetic=False):
        crystal = self.get_crystal()

        lattice = crystal.lattice
        unitcell = lattice.unitcell
        reciprocal = lattice.get_reciprocal(unitcell)

        positions = crystal.get_atom_positions()
        crystal_coordinates = crystal.get_crystal_coordinates(positions, reciprocal)

        starting_magnetization = self.get_starting_magnetization()

        numbers = []
        magnetization = []

        for atom in crystal.atoms:
            symbol = atom.guess_symbol_from_name(atom.name)
            numbers.append(atom.get_atom_number_from_symbol(symbol))
            magnetization.append(
                np.sign(starting_magnetization[atom.name]))

        return unitcell, crystal_coordinates, numbers, magnetization


def _read_from(fpath):
    nml = f90nml.read(fpath)

    with open(fpath) as f:
        lines = f.readlines()

    nat = nml['system']['nat']
    ntyp = nml['system']['ntyp']

    card = {}

    for i, line in enumerate(lines):

        if 'atomic_species' in line.lower():
            card['atomic_species'] = {}

            data = lines[i+1:i+1+ntyp]

            name = []
            mass = []
            pseudo = []

            for a in data:
                name.append(a.split()[0])
                mass.append(float(a.split()[1]))
                pseudo.append(a.split()[2])

            card['atomic_species']['name'] = name
            card['atomic_species']['mass'] = mass
            card['atomic_species']['pseudo'] = pseudo

        elif 'atomic_positions' in line.lower():
            card['atomic_positions'] = {}

            _, option = line.lower().split()
            card['atomic_positions']['option'] = option

            data = lines[i+1:i+1+nat]

            name = []
            position = []

            for a in data:
                name.append(a.split()[0])
                position.append([float(x) for x in a.split()[1:]])

            card['atomic_positions']['name'] = name
            card['atomic_positions']['position'] = position

        elif "cell_parameters" in line.lower():
            card['cell_parameters'] = {}

            _, option = line.lower().split()
            card['cell_parameters']['option'] = option

            data = lines[i+1:i+4]

            unitcell = []

            for a in data:
                unitcell.append([float(x) for x in a.split()])

            card['cell_parameters']['unitcell'] = unitcell

        elif 'k_points' in line.lower():
            card['k_points'] = {}

            _, option = line.lower().split()
            card['k_points']['option'] = option

            if option == 'gamma':
                pass
            elif option == 'automatic':
                data = lines[i+1]
                grid = [int(x) for x in data.split()]
                card['k_points']['grid'] = grid
            else:
                nkpt = int(lines[i+1].strip())
                data = lines[i+2:i+2+nkpt]

                kpts = []
                for a in data:
                    kpts.append([float(x) for x in a.split()])

                card['k_points']['nkpt'] = nkpt
                card['k_points']['list'] = kpts

        elif "atomic_forces" in line.lower():
            card['atomic_forces'] = {}

            data = lines[i+1:i+1+nat]

            name = []
            force = []

            for a in data:
                name.append(a.split()[0])
                force.append([float(x) for x in a.split()[1:]])

            card['atomic_forces']['name'] = name
            card['atomic_forces']['force'] = force

    return nml, card


def print_pw_input(nml, card):
    print(nml)

    headform = "\n{:s} {:s}"

    nat = int(nml['system']['nat'])
    ntyp = int(nml['system']['ntyp'])

    tag = 'atomic_species'
    print(headform.format(tag, "").rstrip())
    for i in range(ntyp):
        x      = card[tag]["name"][i]
        mass   = card[tag]["mass"][i]
        pseudo = card[tag]["pseudo"][i]
        print(f"{x:4s}{mass:8.3f}{'':4s}{pseudo:80s}".rstrip())

    tag = "atomic_positions"
    opt = card[tag]["option"]
    print(headform.format(tag, opt).rstrip())
    for i in range(nat):
        x   = card[tag]["name"][i]
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
        numkpt = card[tag]["nkpt"]
        print(f"{numkpt:4d}")
        for i in range(numkpt):
            kpt = card[tag]["list"][i]
            print("".join(f"{x:16.8f}" for x in kpt).rstrip())

    tag = "cell_parameters"
    opt = card[tag]["option"]
    print(headform.format(tag, opt).rstrip())
    for i in range(3):
        a = card[tag]["unitcell"][i]
        print(f"{a[0]:16.8f}{a[1]:16.8f}{a[2]:16.8f}".rstrip())

