# PROPTIMUS RAPHAN

## Description

PROPTIMUS RAPHAN is a rapid alternative to optimisation with constrained alpha carbons. This approach divides a protein structure into overlapping substructures, allowing each to be optimised independently. As a result, the computation time is linear with respect to the size of the structure. Our approach can achieve results comparable to the overall optimisation of the structure with constrained alpha carbons in significantly less time. PROPTIMUS RAPHAN employs an almost quantum-mechanically accurate force field, [GFN-FF](https://onlinelibrary.wiley.com/doi/10.1002/anie.202004239). This force field is generic, physics-based, and suitable for large molecular systems. The details about the methodology are described on the [wiki](https://github.com/sb-ncbr/proptimus_raphan/wiki).

## Getting Started

**1. Download and run the [Miniconda](https://www.anaconda.com/docs/getting-started/miniconda/main) installer**

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

**2. Clone GitHub repository:**

```
git clone https://github.com/sb-ncbr/proptimus_raphan
```

**3. Set up virtual environment**

* Create a virtual environment named `proptimus_raphan_env`:

```
conda env create -f proptimus_raphan/environment.yml
```

* Activate the virtual environment:

```
conda activate proptimus_raphan_env
```

## Executing the calculation

### Required arguments

`--PDB_file PDB_FILE`          PDB file with structure, which should be optimised.

`--data_dir DATA_DIR`          Directory for saving results.

### Optional arguments

`--cpu CPU`                                       How many CPUs should be used for the calculation.

`--delete_auxiliary_files`                        Auxiliary calculation files can be large. With this argument, the auxiliary files will be deleted during the calculation. Do not use in combination with the argument `‑‑constrained_alpha_carbons_optimisations`!

`--constrained_alpha_carbons_optimisations`       To test the methodology. The original structure and the PROPTIMUS RAPHAN–optimised structure will be optimised with constrained alpha carbons. Short comparison will be stored in <data_dir>/comparison.json. Please note that optimisation with constrained alpha carbons is computationally expensive for larger protein structures.

### Example of executing the calculation:

```bash
cd proptimus_raphan
python raphan.py --PDB_file examples/P0DL07.pdb --data_dir P0DL07_test
```
For testing purposes, you can compare the optimised structure with the reference:

```bash
diff examples/P0DL07_optimised.pdb P0DL07_test/optimised_PDB/P0DL07_optimised.pdb
```
Please note that the resulting structures may vary slightly. This is due to numerical instabilities caused by running on different hardware.

## License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/sb-ncbr/proptimus_raphan/blob/main/LICENSE) file for details.
