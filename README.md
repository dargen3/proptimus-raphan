# PROPTIMus-RAPHAN

Per Residuum Optimisation - Rapid Alternative to Protein Structure Optimisation with Constrained Alpha Carbons (PROPTIMus-RAPHAN) je python program implementující algoritmus Per Residuum Optimisation pro rychlou optimalizaci proteinových sturktur s omezenemýni alpha uhlíky. 

## How to run

```bash
$  python raphan.py --PDB_file examples/P0DL07.pdb --data_dir P0DL07_test --run_full_xtb_optimisation
```
Struktury examples/P0DL07_optimised.pdb a P0DL07_test/optimised_structure/P0DL07_optimised_final.pdb by měly být totožné.

V adresáři P0DL07 by měl být soubor data.json, který by měl vypadat nějak takto :
```
{
    "raphan time": 2.480051279067993,       # může se mírně lišit
    "xtb(original) time": 1.620025396347046,   # může se mírně lišit
    "xtb(raphan) time": 0.7602341175079346,  # může se mírně lišit
    "original/xtb(original) MAD": 0.5124319791793823,     # mělo by být stejné
    "proptimus/xtb(raphan) MAD": 0.028605064377188683,    # mělo by být stejné
    "xtb(raphan)/xtb(original)": 0.09751877188682556      # mělo by být stejné
}
```
