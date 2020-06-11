import f90nml
from collections import namedtuple
import numpy as np
import spglib

from base.atom import Atom, guess_symbol_from_name, get_number_from_symbol
from base.lattice import Lattice, get_reciprocal
from base.crystal import Crystal, get_crystal_coordinate
from symm.equivatom import get_equiv_atom_map


AtomicSpecies = namedtuple('AtomicSpecies', 'name, mass, pseudo')
AtomicPositions = namedtuple('AtomicPositions', 'option, name, position')
AtomicForces = namedtuple('AtomicForces', 'name, force')
CellParameters = namedtuple('CellParameters', 'option, unitcell')


class Pwin:
    def __init__(self, fpath):
        nml, card = _read_from(fpath)

        Control = namedtuple('Control', nml['control'].keys())
        System = namedtuple('System', nml['system'].keys())
        Electrons = namedtuple('Electrons', nml['electrons'].keys())

        # control and system namelists are required
        self.control = Control(**nml['control'])
        self.system = System(**nml['system'])

        try:
            self.electrons = Electrons(**nml['electrons'])
        except:
            self.electrons = None

        # atomic_species and atomic_positions cards are required
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
        atomlist_name = self.atomic_positions.name
        atomlist_position = self.atomic_positions.position

        lattice = Lattice(np.array(unitcell))

        atomlist = []
        for name, position in zip(atomlist_name, atomlist_position):
            atomlist.append(Atom(name, np.array(position)))
        
        return Crystal(lattice, atomlist)

    def get_is_hub(self):
        num_atomtype = self.system.ntyp
        atomtype_name = self.atomic_species.name
        atomtype_hubu = self.system.hubbard_u

        if num_atomtype < len(atomtype_hubu):
            raise ValueError('ntyp must be greater than the number of hubbard atoms.')
        else:
            is_hub = {}
            for i, name in enumerate(atomtype_name):
                try:
                    u = atomtype_hubu[i]
                except IndexError:
                    u = None
                if not u:
                    is_hub[name] = False
                else:
                    is_hub[name] = True

        return is_hub

    def get_num_hub(self):
        is_hub = self.get_is_hub()
        atomlist_name = self.atomic_positions.name

        num_hub = 0
        for name in atomlist_name:
            if is_hub[name]:
                num_hub += 1

        return num_hub

    def get_starting_magnetization(self):
        num_atomtype = self.system.ntyp
        atomtype_name = self.atomic_species.name

        _starting_magnetization = self.system.starting_magnetization

        if num_atomtype < len(_starting_magnetization):
            raise ValueError('ntyp must be greater than the number of hubbard atoms.')
        else:
            starting_magnetization = {}
            for i, name in enumerate(atomtype_name):
                try:
                    m_i = _starting_magnetization[i]
                except IndexError:
                    m_i = 0.0

                starting_magnetization[name] = m_i

        return starting_magnetization

    def get_spglib_cell(self):
        crystal = self.get_crystal()

        unitcell = crystal.lattice.unitcell
        reciprocal = get_reciprocal(unitcell)

        position = crystal.get_atomlist_position()
        crystal_coordinate = get_crystal_coordinate(position, reciprocal)

        starting_magnetization = self.get_starting_magnetization()

        atomic_number = []
        magnetization = []

        for atom in crystal.atomlist:
            symbol = guess_symbol_from_name(atom.name)
            atomic_number.append(get_number_from_symbol(symbol))
            magnetization.append(
                np.sign(starting_magnetization[atom.name]))

        return unitcell, crystal_coordinate, atomic_number, magnetization

    def get_symmetry_data(self):
        spglib_cell = self.get_spglib_cell()
        return spglib.get_symmetry(spglib_cell)

    def get_symmetrizing_map(self, index_reference_atom):
        crystal = self.get_crystal()
        crystal_coordinate = crystal.get_crystal_coordinate()
        
        symmetry_data = self.get_symmetry_data()
        rotation = symmetry_data["rotations"]
        translation = symmetry_data["translations"]

        symmetrizing_map = get_equiv_atom_map(
            crystal_coordinate, rotation, translation, index_reference_atom)

        return symmetrizing_map


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


def _print_pw_input(nml, card):
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

