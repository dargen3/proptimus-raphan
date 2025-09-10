# Per Residuum Optimisation - Rapid Alternative to Protein Structure Optimisation with Constrained Alpha Carbons (PROPTIMuS-RAPHAN)

PROPTIMUS-RAPHAN je python program implementující algoritmus Per Residuum Optimisation pro rychlou optimalizaci proteinových sturktur s omezenemýni alpha uhlíky. 

## How to run

```bash
$  python ppropt.py --cpu 16 --PDB_file examples/P0DL07.pdb --data_dir P0DL07_test --run_full_xtb_optimisation
```

Terminál by měl vypsat něco jako:

original/xtb_original difference: 0.5124319791793823
proptimus/xtb_proptimus difference: 0.028605064377188683
xtb_original/xtb_proptimus difference: 0.09751877188682556

a struktury examples/P0DL07_optimised.pdb a P0DL07_test/optimised_structure/P0DL07_optimised_final.pdb by měly být totožné.

