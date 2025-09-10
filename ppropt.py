import argparse
import json
from math import dist
from multiprocessing import Pool
from os import system, path
from time import time

import tqdm
from Bio import SeqUtils
from Bio.PDB import Select, PDBIO, PDBParser, Superimposer, NeighborSearch
from rdkit import Chem


def load_arguments():
    print("\nParsing arguments... ", end="")
    parser = argparse.ArgumentParser()
    parser.add_argument('--PDB_file',
                        type=str,
                        required=True,
                        help='PDB file with structure, which should be optimised.')
    parser.add_argument('--data_dir',
                        type=str,
                        required=True,
                        help='Directory for saving results.')
    parser.add_argument('--cpu',
                        type=int,
                        required=False,
                        default=1,
                        help='How many CPUs should be used for the calculation.')
    parser.add_argument('--run_full_xtb_optimisation',
                        action="store_true",
                        help='For testing the methodology. It also runs full xtb optimization with alpha carbons constrained.')
    parser.add_argument("--delete_auxiliary_files",
                        action="store_true",
                        help="Auxiliary calculation files can be large. With this argument, "
                             "the auxiliary files will be deleted during the calculation.")

    args = parser.parse_args()
    if not path.isfile(args.PDB_file):
        print(f"\nERROR! File {args.PDB_file} does not exist!\n")
        exit()
    if path.exists(args.data_dir):
        exit(f"\n\nError! Directory with name {args.data_dir} exists. "
             f"Remove existed directory or change --data_dir argument.")
    print("ok")
    return args


class AtomSelector(Select):
    """
    Support class for Biopython.
    After initialization, a set with all full ids of the atoms to be written into the substructure must be stored in self.full_ids.
    """
    def accept_atom(self, atom):
        return int(atom.full_id in self.full_ids)


def optimise_substructure(archive,
                          structure,
                          optimised_residue_index,
                          data_dir,
                          optimised_coordinates,
                          iteration):
    t1 = time() # todo smazat

    for atom, optimised_position in zip(structure.get_atoms(), optimised_coordinates):
        atom.coord = optimised_position


    io = PDBIO()
    io.set_structure(structure)
    substructure_data_dir = f"{data_dir}/sub_{optimised_residue_index+1}"

    kdtree = NeighborSearch(list(structure.get_atoms())) # optimalizovat todo
    optimised_residue = list(structure.get_residues())[optimised_residue_index]

    # define optimized atoms
    optimised_atoms = []
    for atom in optimised_residue:
        if atom.get_parent().id[1] == 1:  # in first residue optimise also -NH3 function group
            optimised_atoms.append(atom)
        else:
            if atom.name not in ["N", "H"]:
                optimised_atoms.append(atom)
    for atom in ["N", "H"]: # add atoms from following bonded residue to optimise whole peptide bond
        try:
            optimised_atoms.append(structure[0]["A"][optimised_residue.id[1]+1][atom])
        except KeyError: # because of last residue
            break

    atoms_in_7A = []
    atoms_in_12A = []
    for optimised_residue_atom in optimised_atoms:
        atoms_in_7A.extend(kdtree.search(center=optimised_residue_atom.coord,
                                         radius=7,
                                         level="A"))
        atoms_in_12A.extend(kdtree.search(center=optimised_residue_atom.coord,
                                         radius=12,
                                         level="A"))

    selector = AtomSelector()
    selector.full_ids = set([atom.full_id for atom in atoms_in_7A])
    io.save(file=f"{substructure_data_dir}/atoms_in_7A.pdb",
                 select=selector,
                 preserve_atom_numbering=True)
    selector.full_ids = set([atom.full_id for atom in atoms_in_12A])
    io.save(file=f"{substructure_data_dir}/atoms_in_12A.pdb",
                 select=selector,
                 preserve_atom_numbering=True)

    # load substructures by RDKit to determine bonds
    mol_min_radius = Chem.MolFromPDBFile(pdbFileName=f"{substructure_data_dir}/atoms_in_7A.pdb",
                                         removeHs=False,
                                         sanitize=False)
    mol_min_radius_conformer = mol_min_radius.GetConformer()
    mol_max_radius = Chem.MolFromPDBFile(pdbFileName=f"{substructure_data_dir}/atoms_in_12A.pdb",
                                         removeHs=False,
                                         sanitize=False)
    mol_max_radius_conformer = mol_max_radius.GetConformer()

    # dictionaries allow quick and precise matching of atoms from mol_min_radius and mol_max_radius
    mol_min_radius_coord_dict = {}
    for i, mol_min_radius_atom in enumerate(mol_min_radius.GetAtoms()):
        coord = mol_min_radius_conformer.GetAtomPosition(i)
        mol_min_radius_coord_dict[(coord.x, coord.y, coord.z)] = mol_min_radius_atom
    mol_max_radius_coord_dict = {}
    for i, mol_max_radius_atom in enumerate(mol_max_radius.GetAtoms()):
        coord = mol_max_radius_conformer.GetAtomPosition(i)
        mol_max_radius_coord_dict[(coord.x, coord.y, coord.z)] = mol_max_radius_atom

    # find atoms from mol_min_radius with broken bonds
    atoms_with_broken_bonds = []
    for mol_min_radius_atom in mol_min_radius.GetAtoms():
        coord = mol_min_radius_conformer.GetAtomPosition(mol_min_radius_atom.GetIdx())
        mol_max_radius_atom = mol_max_radius_coord_dict[(coord.x, coord.y, coord.z)]
        if len(mol_min_radius_atom.GetNeighbors()) != len(mol_max_radius_atom.GetNeighbors()):
            atoms_with_broken_bonds.append(mol_max_radius_atom)

    # create a substructure that will have only C-C bonds broken
    carbons_with_broken_bonds_coord = []  # hydrogens will be added only to these carbons
    substructure_coord_dict = mol_min_radius_coord_dict
    while atoms_with_broken_bonds:
        atom_with_broken_bonds = atoms_with_broken_bonds.pop(0)
        bonded_atoms = atom_with_broken_bonds.GetNeighbors()
        for bonded_atom in bonded_atoms:
            coord = mol_max_radius_conformer.GetAtomPosition(bonded_atom.GetIdx())
            if (coord.x, coord.y, coord.z) in substructure_coord_dict:
                continue
            else:
                if atom_with_broken_bonds.GetSymbol() == "C" and bonded_atom.GetSymbol() == "C":
                    carbons_with_broken_bonds_coord.append(
                        mol_max_radius_conformer.GetAtomPosition(atom_with_broken_bonds.GetIdx()))
                    continue
                else:
                    atoms_with_broken_bonds.append(bonded_atom)
                    substructure_coord_dict[(coord.x, coord.y, coord.z)] = bonded_atom

    # create substructure in Biopython library
    substructure_atoms = [kdtree.search(center=coord,
                                        radius=0.1,
                                        level="A")[0] for coord in substructure_coord_dict.keys()]
    selector.full_ids = set([atom.full_id for atom in substructure_atoms])
    io.save(file=f"{substructure_data_dir}/substructure_{iteration}.pdb",
            select=selector,
            preserve_atom_numbering=True)
    substructure = PDBParser(QUIET=True).get_structure(id="structure",
                                                       file=f"{substructure_data_dir}/substructure_{iteration}.pdb")
    substructure_atoms = list(substructure.get_atoms())


    # definitions of which atoms should be constrained during optimization
    # optimized atoms are not constrained during optimisation and written into the overall structure
    # flexible atoms are not constrained during optimisation and are not written into the overall structure
    # constrained atoms are constrained during optimisation and are not written into the overall structure
    constrained_atoms_indices = []
    optimised_atoms_indices = []
    rigid_atoms = [] # constrained atoms without atoms close to carbons with broken bonds
    rigid_atoms_indices = []
    for i, atom in enumerate(substructure.get_atoms(),
                             start=1):
        if atom.name == "CA":
            atom.mode = "constrained"
        elif atom in optimised_atoms:
            atom.mode = "optimised"
        elif any(dist(atom.coord, ra.coord) < 4 for ra in optimised_atoms):
            atom.mode = "flexible"
        else:
            atom.mode = "constrained"

        if atom.mode == "optimised":
            optimised_atoms_indices.append(i)
        elif atom.mode == "constrained":
            constrained_atoms_indices.append(i)
            if min([dist(atom.coord, x) for x in carbons_with_broken_bonds_coord]) > 2:
                rigid_atoms.append(atom)
                rigid_atoms_indices.append(i)

    # prepare xtb settings file
    xtb_settings_template = f"""$constrain
    atoms: xxx
    force constant=10.0
    $end
    $opt
    maxcycle={len(optimised_atoms_indices)+iteration}
    microcycle={len(optimised_atoms_indices)+iteration + 10}
    $end
    """
    substructure_settings = xtb_settings_template.replace("xxx", ", ".join([str(i) for i in constrained_atoms_indices]))
    with open(f"{substructure_data_dir}/xtb_settings_{iteration}.inp", "w") as xtb_settings_file:
        xtb_settings_file.write(substructure_settings)

    # optimise substructure by xtb
    run_xtb = (f"cd {substructure_data_dir} ;"
               f"ulimit -s unlimited ;"
               f"export OMP_STACKSIZE=1G ; "
               f"export OMP_NUM_THREADS=1,1 ;"
               f"export OMP_MAX_ACTIVE_LEVELS=1 ;"
               f"export MKL_NUM_THREADS=1 ;"
               f"xtb substructure_{iteration}.pdb --gfnff --input xtb_settings_{iteration}.inp --opt vtight --alpb water --verbose > xtb_output_{iteration}.txt 2> xtb_error_output_{iteration}.txt   ; rm gfnff*")
    t2 = time() # todo smazat
    system(run_xtb)
    t3 = time() - t2 # todo smazat

    from pathlib import Path
    file_path = Path(f"{substructure_data_dir}/xtbopt.pdb")
    if not file_path.exists(): # TODO !!!!!!!!
        print(substructure_data_dir, iteration)
        exit()
    system(f"cd {substructure_data_dir} ; mv xtbopt.log xtbopt_{iteration}.log ; mv xtbopt.pdb xtbopt_{iteration}.pdb")


    # superimpose optimised and original substructures
    optimised_substructure = PDBParser(QUIET=True).get_structure("substructure",
                                                                 f"{substructure_data_dir}/xtbopt_{iteration}.pdb")
    optimised_substructure_atoms = list(optimised_substructure.get_atoms())
    optimised_rigid_atoms = [optimised_substructure_atoms[rigid_atom_index - 1] for rigid_atom_index in rigid_atoms_indices]
    sup = Superimposer()
    sup.set_atoms(rigid_atoms, optimised_rigid_atoms)
    sup.apply(optimised_substructure.get_atoms())

    # write coordinates of optimised atoms
    optimised_coordinates = []
    for optimised_atom_index in optimised_atoms_indices:
        optimised_atom_coord = optimised_substructure_atoms[optimised_atom_index - 1].coord
        original_atom_index = substructure_atoms[optimised_atom_index - 1].serial_number - 1
        optimised_coordinates.append((original_atom_index, optimised_atom_coord))

    # check convergence
    converged = False
    if len(archive) > 1:
        max_diffs = []
        for x in range(1, 3):
            diffs = [dist(a,b) for a,b in zip([x[1] for x in optimised_coordinates], archive[-x])]
            max_diffs.append(max(diffs))
        if any([x<0.01 for x in max_diffs]):
            converged = True

    t4 = time() - t1 # todo smazat

    # print(t4, t3)

    return optimised_coordinates, converged, optimised_residue_index



class PRO:
    def __init__(self,
                 data_dir: str,
                 PDB_file: str,
                 cpu: int,
                 delete_auxiliary_files: bool):
        self.data_dir = data_dir
        self.PDB_file = PDB_file
        self.cpu = cpu
        self.delete_auxiliary_files = delete_auxiliary_files

    def optimise(self):
        print(f"Loading of structure from {self.PDB_file}... ", end="")
        self._load_molecule()
        print("ok")

        print("Optimisation...")
        self.optimised_coordinates = [atom.coord for atom in self.structure.get_atoms()]
        with Pool(self.cpu) as pool:
            #
            # for iteration in tqdm.tqdm(range(1, 51),
            #                            desc="Structure optimisation",
            #                            unit="iter.",
            #                            smoothing=0,
            #                            delay=0.1,
            #                            mininterval=0.4,
            #                            maxinterval=0.4):
            converged = [False for _ in self.structure.get_residues()]
            archive = [[] for _ in self.structure.get_residues()]
            for iteration in range(1, 50):
                iteration_results = pool.starmap(optimise_substructure, [(archive[ori], self.structure, ori, self.data_dir, self.optimised_coordinates, iteration) for ori, substructure_converged in enumerate(converged) if not substructure_converged])
                for optimised_coordinates, sc, optimised_residue_index in iteration_results:
                    for optimised_atom_index, optimised_atom_coordinates in optimised_coordinates:
                        self.optimised_coordinates[optimised_atom_index] = optimised_atom_coordinates
                    archive[optimised_residue_index].append([x[1] for x in optimised_coordinates])
                    converged[optimised_residue_index] = sc
                for atom, coord in zip(self.structure.get_atoms(), self.optimised_coordinates):
                    atom.coord = coord
                self.io.save(f"{self.data_dir}/optimised_PDB/{path.basename(self.PDB_file[:-4])}_optimised_{iteration}.pdb")
                if all(converged):
                    print(iteration) # todo hlaška
                    break


        self.io.save(f"{self.data_dir}/optimised_PDB/{path.basename(self.PDB_file[:-4])}_optimised_final.pdb")
        return






        print("Storage of the optimised structure... ", end="")
        logs = sorted([json.loads(open(f).read()) for f in glob(f"{self.data_dir}/sub_*/residue.log")],
                      key=lambda x: x['residue index'])
        atom_counter = 0
        for optimised_residue, log in zip(self.residues, logs):
            d = 0
            for optimised_atom in optimised_residue.get_atoms():
                d += dist(optimised_atom.coord, self.original_atoms_positions[atom_counter])**2
                atom_counter += 1
            residual_rmsd = (d / len(list(optimised_residue.get_atoms())))**(1/2)
            log["residual_rmsd"] = residual_rmsd
            if residual_rmsd > 1:
                log["category"] = "Highly optimised residue"
        with open(f"{self.data_dir}/residues.logs", "w") as residues_logs:
            residues_logs.write(json.dumps(logs, indent=2))
        self.io.save(f"{self.data_dir}/optimised_PDB/{path.basename(self.PDB_file[:-4])}_optimised.pdb")
        if self.delete_auxiliary_files:
            system(f"for au_file in {self.data_dir}/sub_* ; do rm -fr $au_file ; done &")
        print("ok\n\n")

    def _load_molecule(self):
        system(f"mkdir {self.data_dir};"
               f"mkdir {self.data_dir}/inputed_PDB;"
               f"mkdir {self.data_dir}/optimised_PDB;"
               f"cp {self.PDB_file} {self.data_dir}/inputed_PDB")
        try:
            structure = PDBParser(QUIET=True).get_structure("structure", self.PDB_file)
            io = PDBIO()
            io.set_structure(structure)
            self.io = io
            self.structure = io.structure
        except KeyError:
            exit(f"\nERROR! PDB file {self.PDB_file} does not contain any structure.\n")
        self.residues = list(self.structure.get_residues())
        self.atoms = list(self.structure.get_atoms())
        self.original_atoms_positions = [atom.coord for atom in self.structure.get_atoms()]

        # creation of data directories for substructure
        for residue_number in range(1, len(list(self.structure.get_residues()))+1):
            substructure_data_dir = f"{self.data_dir}/sub_{residue_number}"
            system(f"mkdir {substructure_data_dir}")





def run_full_xtb_optimisation(PRO):
    alpha_carbons_indices = []
    structure = PDBParser(QUIET=True).get_structure("structure", PRO.PDB_file)
    for i, atom in enumerate(structure.get_atoms(), start=1):
        if atom.name == "CA":
            alpha_carbons_indices.append(str(i))


    system(f"mkdir {PRO.data_dir}/full_xtb_optimisation ")
    system(f"mkdir {PRO.data_dir}/full_xtb_optimisation/original ")
    with open(f"{PRO.data_dir}/full_xtb_optimisation/original/xtb_settings.inp", "w") as xtb_settings_file:
        xtb_settings_file.write(f"$constrain\n    force constant=10.0\n    atoms: {",".join(alpha_carbons_indices)}\n$end""")

    system(f"""cd {PRO.data_dir}/full_xtb_optimisation/original;
               export OMP_NUM_THREADS=1,1 ;
               export MKL_NUM_THREADS=1 ;
               export OMP_MAX_ACTIVE_LEVELS=1 ;
               export OMP_STACKSIZE=200G ;
               ulimit -s unlimited ;
               xtb ../../inputed_PDB/{path.basename(PRO.PDB_file)} --opt --alpb water --verbose --gfnff --input xtb_settings.inp --verbose > xtb_output.txt 2> xtb_error_output.txt""")



    system(f"mkdir {PRO.data_dir}/full_xtb_optimisation/proptimus ")
    with open(f"{PRO.data_dir}/full_xtb_optimisation/proptimus/xtb_settings.inp", "w") as xtb_settings_file:
        xtb_settings_file.write(f"$constrain\n    force constant=10.0\n    atoms: {",".join(alpha_carbons_indices)}\n$end""")

    system(f"""cd {PRO.data_dir}/full_xtb_optimisation/proptimus ;
               export OMP_NUM_THREADS=1,1 ;
               export MKL_NUM_THREADS=1 ;
               export OMP_MAX_ACTIVE_LEVELS=1 ;
               export OMP_STACKSIZE=200G ;
               ulimit -s unlimited ;
               xtb ../../optimised_PDB/{path.basename(PRO.PDB_file[:-4])}_optimised_final.pdb --opt --alpb water --verbose --gfnff --input xtb_settings.inp --verbose > xtb_output.txt 2> xtb_error_output.txt""")



    s1 = PDBParser(QUIET=True).get_structure(id="structure", file=PRO.PDB_file)
    s2 = PDBParser(QUIET=True).get_structure(id="structure", file=f"{PRO.data_dir}/full_xtb_optimisation/original/xtbopt.pdb")
    sup = Superimposer()
    sup.set_atoms([a for a in s1.get_atoms() if a.name == "CA"], [a for a in s2.get_atoms() if a.name == "CA"])
    sup.apply(s2.get_atoms())
    d = []
    for a1, a2 in zip(s1.get_atoms(), s2.get_atoms()):
        d.append(a1 - a2)
    print(f"\n\n\noriginal/xtb_original difference: {sum(d)/len(d)}")

    s1 = PDBParser(QUIET=True).get_structure(id="structure", file=f"{PRO.data_dir}/optimised_PDB/{path.basename(PRO.PDB_file[:-4])}_optimised_final.pdb")
    s2 = PDBParser(QUIET=True).get_structure(id="structure", file=f"{PRO.data_dir}/full_xtb_optimisation/proptimus/xtbopt.pdb")
    sup = Superimposer()
    sup.set_atoms([a for a in s1.get_atoms() if a.name == "CA"], [a for a in s2.get_atoms() if a.name == "CA"])
    sup.apply(s2.get_atoms())
    d = []
    for a1, a2 in zip(s1.get_atoms(), s2.get_atoms()):
        d.append(a1 - a2)
    print(f"proptimus/xtb_proptimus difference: {sum(d) / len(d)}")



    s1 = PDBParser(QUIET=True).get_structure(id="structure", file=f"{PRO.data_dir}/full_xtb_optimisation/original/xtbopt.pdb")
    s2 = PDBParser(QUIET=True).get_structure(id="structure", file=f"{PRO.data_dir}/full_xtb_optimisation/proptimus/xtbopt.pdb")
    sup = Superimposer()
    sup.set_atoms([a for a in s1.get_atoms() if a.name == "CA"], [a for a in s2.get_atoms() if a.name == "CA"])
    sup.apply(s2.get_atoms())
    d = []
    for a1, a2 in zip(s1.get_atoms(), s2.get_atoms()):
        d.append(a1 - a2)
    print(f"xtb_original/xtb_proptimus difference: {sum(d) / len(d)}")




if __name__ == '__main__':
    args = load_arguments()
    t = time()
    Proptimus = PRO(args.data_dir, args.PDB_file, args.cpu, args.delete_auxiliary_files)
    Proptimus.optimise()
    proptimus_time = time() - t
    print(proptimus_time)

    if args.run_full_xtb_optimisation:
        print("Running full xtb optimisation...", end="")
        t = time()
        run_full_xtb_optimisation(Proptimus)
        full_optimisation_time = time() - t




# dopsat rm auxiliary files

# co se děje když to nezkonverguje?

# na závěr doplnit pet iterace

# omezit počet iterací na 50

"""
P0DL07     0.028  2
B4U768	   0.101  50 
Q980V5	   0.037  59
P83999	   0.026  36
P62103	   0.036  20
P06123	   0.037  50
A0RLG5     0.053  172
A3DEN2     0.039  94
B2UFK5     0.028  175
B9K6M6     0.045  216
O69250     0.025  103
P00239     0.041  119
P00248     0.048  134
Q7MQC4     0.022  184
"""


#jak spustit
#testovací data
#co to má vyplivnout
# michal mikuš