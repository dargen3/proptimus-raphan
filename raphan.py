#!/usr/bin/env python3
import argparse
from math import dist
from multiprocessing import Pool
from os import system, path
from pathlib import Path
from time import time
import json

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
    parser.add_argument('--constrained_alpha_carbons_optimisations',
                        action="store_true",
                        help='For testing the methodology. '
                             'The original structure and the PROPTIMUS RAPHANgfnff optimized structure will be optimized with constrained alpha carbons.'
                             'Short comparison will be stored in <data_dir>/comparison.json. '
                             'Please note that optimization with constrained alpha carbons is computationally expensive for larger protein structures.')
    parser.add_argument("--delete_auxiliary_files",
                        action="store_true",
                        help="Auxiliary calculation files can be large. With this argument, "
                             "the auxiliary files will be deleted during the calculation."
                             "Do not use in combination with the argument --constrained_alpha_carbons_optimisations!")

    args = parser.parse_args()
    if not path.isfile(args.PDB_file):
        print(f"\nERROR! File {args.PDB_file} does not exist!\n")
        exit()
    if path.exists(args.data_dir):
        exit(f"\nError! Directory with name {args.data_dir} exists. "
             f"Remove existed directory or change --data_dir argument.")
    print("ok")
    return args


class AtomSelector(Select):
    """
    Support class for Biopython.
    After initialization, a set with all full ids of the atoms and set with all full ids of the residues to be
    written into the substructure must be stored in self.full_ids and self.res_full_ids.
    """
    def accept_atom(self, atom):
        return int(atom.full_id in self.full_ids)

    def accept_residue(self, residue):
        return int(residue.full_id in self.res_full_ids)


class Substructure_data:
    def __init__(self,
                 atoms_30A,
                 data_dir,
                 optimised_residue,
                 optimised_residue_index,
                 structure,
                 io):
        self.atoms_30A = atoms_30A
        self.data_dir = data_dir
        self.optimised_residue_index = optimised_residue_index
        self.optimised_residue = optimised_residue
        self.archive = []
        self.converged = False
        self.io = io

        self.optimised_atoms = set()
        for atom in optimised_residue:
            if atom.get_parent().id[1] == 1:  # in first residue optimise also -NH3 function group
                self.optimised_atoms.add(atom)
            else:
                if atom.name not in ["N", "H"]:
                    self.optimised_atoms.add(atom)
        for atom in ["N", "H"]:  # add atoms from following bonded residue to optimise whole peptide bond
            try:
                self.optimised_atoms.add(structure[0]["A"][optimised_residue.id[1] + 1][atom])
            except KeyError:  # because of last residue
                break

        self.final_optimised_atoms = set([atom for atom in optimised_residue])


def optimise_substructure(substructure_data,
                          iteration,
                          phase):

    if phase == "optimisation":
        optimised_atoms = substructure_data.optimised_atoms
        minimum_radius = 6
    elif phase == "final refinement":
        optimised_atoms = substructure_data.final_optimised_atoms
        minimum_radius = 8

    # find effective neighbourhood
    kdtree = NeighborSearch(substructure_data.atoms_30A)
    atoms_in_minimum_radius = []
    atoms_in_12A = []
    for optimised_residue_atom in optimised_atoms:

        atoms_in_minimum_radius.extend(kdtree.search(center=optimised_residue_atom.coord,
                                         radius=minimum_radius,
                                         level="A"))
        atoms_in_12A.extend(kdtree.search(center=optimised_residue_atom.coord,
                                         radius=12,
                                         level="A"))

        

    # create pdb files to can load them with RDKit
    selector = AtomSelector()
    selector.full_ids = set([atom.full_id for atom in atoms_in_minimum_radius])
    selector.res_full_ids = set([atom.get_parent().full_id for atom in atoms_in_minimum_radius])
    substructure_data.io.save(file=f"{substructure_data.data_dir}/atoms_in_minimum_radius.pdb",
                              select=selector,
                              preserve_atom_numbering=True)
    selector.full_ids = set([atom.full_id for atom in atoms_in_12A])
    selector.res_full_ids = set([atom.get_parent().full_id for atom in atoms_in_12A])
    substructure_data.io.save(file=f"{substructure_data.data_dir}/atoms_in_12A.pdb",
                              select=selector,
                              preserve_atom_numbering=True)

    # load pdb files with RDKit
    mol_min_radius = Chem.MolFromPDBFile(pdbFileName=f"{substructure_data.data_dir}/atoms_in_minimum_radius.pdb",
                                         removeHs=False,
                                         sanitize=False)
    mol_min_radius_conformer = mol_min_radius.GetConformer()
    mol_max_radius = Chem.MolFromPDBFile(pdbFileName=f"{substructure_data.data_dir}/atoms_in_12A.pdb",
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
    carbons_with_broken_bonds_coord = []
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
                    carbons_with_broken_bonds_coord.append(mol_max_radius_conformer.GetAtomPosition(atom_with_broken_bonds.GetIdx()))
                    continue
                else:
                    atoms_with_broken_bonds.append(bonded_atom)
                    substructure_coord_dict[(coord.x, coord.y, coord.z)] = bonded_atom

    # create substructure in Biopython library
    substructure_atoms = [kdtree.search(center=coord,
                                        radius=0.1,
                                        level="A")[0] for coord in substructure_coord_dict.keys()]
    selector.full_ids = set([atom.full_id for atom in substructure_atoms])
    selector.res_full_ids = set([atom.get_parent().full_id for atom in substructure_atoms])
    substructure_data.io.save(file=f"{substructure_data.data_dir}/substructure_{iteration}.pdb",
                              select=selector,
                              preserve_atom_numbering=True)
    substructure = PDBParser(QUIET=True).get_structure(id="structure",
                                                       file=f"{substructure_data.data_dir}/substructure_{iteration}.pdb")
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
            if any(dist(atom.coord, x) < 2 for x in carbons_with_broken_bonds_coord):
                continue
            else:
                rigid_atoms.append(atom)
                rigid_atoms_indices.append(i)

    # prepare xtb settings file
    xtb_settings_template = f"""$constrain
    atoms: xxx
    force constant=10.0
    $end
    $opt
    maxcycle={len(optimised_atoms_indices)+iteration}
    microcycle={len(optimised_atoms_indices)+iteration+1}
    $end
    """
    substructure_settings = xtb_settings_template.replace("xxx", ", ".join([str(i) for i in constrained_atoms_indices]))
    with open(f"{substructure_data.data_dir}/xtb_settings_{iteration}.inp", "w") as xtb_settings_file:
        xtb_settings_file.write(substructure_settings)

    # optimise substructure by xtb
    run_xtb = (f"cd {substructure_data.data_dir} ;"
               f"ulimit -s unlimited ;"
               f"export OMP_STACKSIZE=1G ; "
               f"export OMP_NUM_THREADS=1,1 ;"
               f"export OMP_MAX_ACTIVE_LEVELS=1 ;"
               f"export MKL_NUM_THREADS=1 ;"
               f" xtb substructure_{iteration}.pdb --gfnff --input xtb_settings_{iteration}.inp --opt tight --alpb water --verbose > xtb_output_{iteration}.txt 2> xtb_error_output_{iteration}.txt ; rm gfnff*")
    system(run_xtb)

    # check xtb convergence
    if not Path(f"{substructure_data.data_dir}/xtbopt.pdb").exists():
        return None, False, None
    system(f"cd {substructure_data.data_dir} ; mv xtbopt.log xtbopt_{iteration}.log ; mv xtbopt.pdb xtbopt_{iteration}.pdb")

    # superimpose original and optimised structures
    optimised_substructure = PDBParser(QUIET=True).get_structure("substructure",
                                                                 f"{substructure_data.data_dir}/xtbopt_{iteration}.pdb")
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

    central_residue_coordinates = []
    for atom in optimised_substructure_atoms:
        if atom.get_parent().id[1] == substructure_data.optimised_residue_index:
            central_residue_coordinates.append(atom.coord)

    # check raphan convergence
    raphan_converged = False
    if len(substructure_data.archive) > 1:
        max_diffs = []
        for x in range(1, 3):
            diffs = [dist(a,b) for a,b in zip([x[1] for x in optimised_coordinates], substructure_data.archive[-x])]
            max_diffs.append(max(diffs))
        if any([x<0.01 for x in max_diffs]):
            raphan_converged = True
    return optimised_coordinates, raphan_converged, substructure_data


class Raphan:
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
        self._load_molecule()

        # optimise
        self.optimised_coordinates = [atom.coord for atom in self.structure.get_atoms()]
        with Pool(self.cpu) as pool:
            bar = tqdm.tqdm(total=100,
                            desc="Structure optimisation",
                            unit=" iteration")
            for iteration in range(1, 50):
                bar.update(1)
                iteration_results = pool.starmap(optimise_substructure, [(substructure, iteration, "optimisation") for substructure in self.substructures_data if not substructure.converged])
                for optimised_coordinates, convergence, substructure_data in iteration_results:
                    if optimised_coordinates is None: # xtb did not converge
                        continue
                    for optimised_atom_index, optimised_atom_coordinates in optimised_coordinates:
                        self.optimised_coordinates[optimised_atom_index] = optimised_atom_coordinates
                    self.substructures_data[substructure_data.optimised_residue_index-1].archive.append([x[1] for x in optimised_coordinates])
                    self.substructures_data[substructure_data.optimised_residue_index-1].converged = convergence
                for atom, coord in zip(self.structure.get_atoms(), self.optimised_coordinates):
                    atom.coord = coord
                self.io.save(f"{self.data_dir}/optimised_PDB/{path.basename(self.PDB_file[:-4])}_optimised_{iteration}.pdb")
                if all([substructure_data.converged for substructure_data in self.substructures_data]):
                    break

            for substructure_data in self.substructures_data:
                substructure_data.converged = False

            for iteration in range(iteration+1, iteration+51):
                bar.update(1)
                iteration_results = pool.starmap(optimise_substructure, [(substructure, iteration, "final refinement") for substructure in self.substructures_data if not substructure.converged])
                for optimised_coordinates, convergence, substructure_data in iteration_results:
                    if optimised_coordinates is None: # xtb did not converge
                        continue
                    for optimised_atom_index, optimised_atom_coordinates in optimised_coordinates:
                        self.optimised_coordinates[optimised_atom_index] = optimised_atom_coordinates
                    self.substructures_data[substructure_data.optimised_residue_index-1].archive.append([x[1] for x in optimised_coordinates])
                    self.substructures_data[substructure_data.optimised_residue_index-1].converged = convergence
                for atom, coord in zip(self.structure.get_atoms(), self.optimised_coordinates):
                    atom.coord = coord
                self.io.save(f"{self.data_dir}/optimised_PDB/{path.basename(self.PDB_file[:-4])}_optimised_{iteration}.pdb")
                if all([substructure_data.converged for substructure_data in self.substructures_data]):
                    bar.update(100 - iteration)
                    bar.refresh()
                    break
            bar.close()
            self.iterations = iteration

            # control of unconverged residues
            unconverged_substructures = [str(substructure_data.optimised_residue_index) for substructure_data in self.substructures_data if not substructure_data.converged]
            if unconverged_substructures:
                print(f"WARNING! OPTIMISATION FOR RESIDUES WITH INDICE(S) {', '.join(unconverged_substructures)} DID NOT CONVERGE!")

        print(f"Saving optimised structure to {self.data_dir}/optimised_PDB/{path.basename(self.PDB_file[:-4])}_optimised.pdb... ", end="")
        self.io.save(f"{self.data_dir}/optimised_PDB/{path.basename(self.PDB_file[:-4])}_optimised.pdb")
        print("ok")

        if self.delete_auxiliary_files:
            print("Deleting auxiliary files...", end="")
            system(f"cd {self.data_dir};"
                   f"mv optimised_PDB/{path.basename(self.PDB_file[:-4])}_optimised.pdb .;"
                   f"rm -r sub_* optimised_PDB input_PDB")
            print(" ok")


    def _load_molecule(self):
        print(f"Loading of structure from {self.PDB_file}... ", end="")

        # open PDB file by Biopython
        try:
            structure = PDBParser(QUIET=True).get_structure("structure", self.PDB_file)
            io = PDBIO()
            io.set_structure(structure)
            self.io = io
            self.structure = io.structure
        except KeyError:
            exit(f"\nERROR! PDB file {self.PDB_file} does not contain any structure or file is corrupted.\n")

        # prepared substructures data for optimisation
        self.substructures_data = []
        kdtree = NeighborSearch(list(self.structure.get_atoms()))
        for residue_index, residue in enumerate(self.structure.get_residues(), start=1):
            atoms_30A = kdtree.search(center=residue.center_of_mass(geometric=True),
                                      radius=30, # radius of AMK (6A) + outer substructure radius (12A) + maximum shift of atom (10A) + extra (2A)
                                      level="A")
            self.substructures_data.append(Substructure_data(atoms_30A=atoms_30A,
                                                             data_dir=f"{self.data_dir}/sub_{residue_index}",
                                                             optimised_residue=residue,
                                                             optimised_residue_index=residue_index,
                                                             structure=self.structure,
                                                             io=self.io))

        # creation of data directories
        system(f"mkdir {self.data_dir};"
               f"mkdir {self.data_dir}/input_PDB;"
               f"mkdir {self.data_dir}/optimised_PDB;"
               f"cp {self.PDB_file} {self.data_dir}/input_PDB")
        for residue_number in range(1, len(list(self.structure.get_residues()))+1):
            system(f"mkdir {self.data_dir}/sub_{residue_number}")
        print("ok")


def run_constrained_alpha_optimisations(raphan):
    print("Running constrained alpha optimisation...", end="")

    # find alpha_carbons_indices to constrain them
    alpha_carbons_indices = []
    structure = PDBParser(QUIET=True).get_structure("structure", raphan.PDB_file)
    for i, atom in enumerate(structure.get_atoms(), start=1):
        if atom.name == "CA":
            alpha_carbons_indices.append(str(i))

    # optimise original structure
    system(f"mkdir {raphan.data_dir}/constrained_alpha_carbons_optimisations ")
    system(f"mkdir {raphan.data_dir}/constrained_alpha_carbons_optimisations/original ")
    with open(f"{raphan.data_dir}/constrained_alpha_carbons_optimisations/original/xtb_settings.inp", "w") as xtb_settings_file:
        xtb_settings_file.write(f"$constrain\n    force constant=10.0\n    atoms: {','.join(alpha_carbons_indices)}\n$end")
    t = time()
    system(f"""cd {raphan.data_dir}/constrained_alpha_carbons_optimisations/original;
               export OMP_NUM_THREADS=1,1 ;
               export MKL_NUM_THREADS=1 ;
               export OMP_MAX_ACTIVE_LEVELS=1 ;
               export OMP_STACKSIZE=5G ;
               ulimit -s unlimited ;
               xtb ../../input_PDB/{path.basename(raphan.PDB_file)} --opt --alpb water --verbose --gfnff --input xtb_settings.inp --verbose > xtb_output.txt 2> xtb_error_output.txt""")
    GFNFFca_time = time() - t

    # optimise structure already optimised with raphan
    system(f"mkdir {raphan.data_dir}/constrained_alpha_carbons_optimisations/raphan ")
    with open(f"{raphan.data_dir}/constrained_alpha_carbons_optimisations/raphan/xtb_settings.inp", "w") as xtb_settings_file:
        xtb_settings_file.write(f"$constrain\n    force constant=10.0\n    atoms: {','.join(alpha_carbons_indices)}\n$end")
    system(f"""cd {raphan.data_dir}/constrained_alpha_carbons_optimisations/raphan ;
               export OMP_NUM_THREADS=1,1 ;
               export MKL_NUM_THREADS=1 ;
               export OMP_MAX_ACTIVE_LEVELS=1 ;
               export OMP_STACKSIZE=5G ;
               ulimit -s unlimited ;
               xtb ../../optimised_PDB/{path.basename(raphan.PDB_file[:-4])}_optimised.pdb --opt --alpb water --verbose --gfnff --input xtb_settings.inp --verbose > xtb_output.txt 2> xtb_error_output.txt""")

    # compare original structure with structure optimised by GFNFFca
    try:
        s1 = PDBParser(QUIET=True).get_structure(id="structure", file=raphan.PDB_file)
        s2 = PDBParser(QUIET=True).get_structure(id="structure", file=f"{raphan.data_dir}/constrained_alpha_carbons_optimisations/original/xtbopt.pdb")
        sup = Superimposer()
        sup.set_atoms([a for a in s1.get_atoms() if a.name == "CA"], [a for a in s2.get_atoms() if a.name == "CA"])
        sup.apply(s2.get_atoms())
        differences = []
        for a1, a2 in zip(s1.get_atoms(), s2.get_atoms()):
            differences.append(a1 - a2)
        original_GFNFFca_difference = round(float(sum(differences)/len(differences)), 4)
    except FileNotFoundError:
        original_GFNFFca_difference = None

    # compare structure optimised by PROPTIMUS RAPHANgfnff and structure optimised by PROPTIMUS RAPHANgfnff + GFNFFca
    try:
        s1 = PDBParser(QUIET=True).get_structure(id="structure", file=f"{raphan.data_dir}/optimised_PDB/{path.basename(raphan.PDB_file[:-4])}_optimised.pdb")
        s2 = PDBParser(QUIET=True).get_structure(id="structure", file=f"{raphan.data_dir}/constrained_alpha_carbons_optimisations/raphan/xtbopt.pdb")
        sup = Superimposer()
        sup.set_atoms([a for a in s1.get_atoms() if a.name == "CA"], [a for a in s2.get_atoms() if a.name == "CA"])
        sup.apply(s2.get_atoms())
        differences = []
        for a1, a2 in zip(s1.get_atoms(), s2.get_atoms()):
            differences.append(a1 - a2)
        PROPTIMUS_RAPHANgfnff__PROPTIMUS_RAPHANgfnff_GFNFFca_difference = round(float(sum(differences) / len(differences)), 4)
    except:
        PROPTIMUS_RAPHANgfnff__PROPTIMUS_RAPHANgfnff_GFNFFca_difference = None


    # compare structure optimised by GFNFFca and structure optimised by PROPTIMUS RAPHANgfnff + GFNFFca
    try:
        s1 = PDBParser(QUIET=True).get_structure(id="structure", file=f"{raphan.data_dir}/constrained_alpha_carbons_optimisations/original/xtbopt.pdb")
        s2 = PDBParser(QUIET=True).get_structure(id="structure", file=f"{raphan.data_dir}/constrained_alpha_carbons_optimisations/raphan/xtbopt.pdb")
        sup = Superimposer()
        sup.set_atoms([a for a in s1.get_atoms() if a.name == "CA"], [a for a in s2.get_atoms() if a.name == "CA"])
        sup.apply(s2.get_atoms())
        differences = []
        for a1, a2 in zip(s1.get_atoms(), s2.get_atoms()):
            differences.append(a1 - a2)
        GFNFFca__PROPTIMUS_RAPHANgfnff_GFNFFca_difference = round(float(sum(differences) / len(differences)), 4)
    except:
        GFNFFca__PROPTIMUS_RAPHANgfnff_GFNFFca_difference = None

    # write results
    results = {"PROPTIMUS RAPHANgfnff time": round(raphan.calculation_time, 4),
               "PROPTIMUS RAPHANgfnff iterations": raphan.iterations,
               "GFN-FFca time": round(GFNFFca_time, 4),
               "original / GFN-FFca MAD": original_GFNFFca_difference,
               "PROPTIMUS RAPHANgfnff / PROPTIMUS RAPHANgfnff + GFN-FFca MAD": PROPTIMUS_RAPHANgfnff__PROPTIMUS_RAPHANgfnff_GFNFFca_difference,
               "GFN-FFca / PROPTIMUS RAPHANgfnff + GFN-FFca": GFNFFca__PROPTIMUS_RAPHANgfnff_GFNFFca_difference}
    with open(f"{raphan.data_dir}/comparison.json", 'w') as data_json:
        json.dump(results, data_json, indent=4)
    print(" ok")


if __name__ == '__main__':
    args = load_arguments()
    t = time()
    raphan = Raphan(args.data_dir, args.PDB_file, args.cpu, args.delete_auxiliary_files)
    raphan.optimise()
    raphan.calculation_time = time() - t

    if args.constrained_alpha_carbons_optimisations:
        run_constrained_alpha_optimisations(raphan)
    print("\n")
