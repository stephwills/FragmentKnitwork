
### Fragalysis formatted files
This assumes you have a `FRAGALYSIS_DATA_DIR` which contains target directories that have been downloaded from Fragalysis. This adopts the old style of formatting that was used in Fragalysis, whereby fragments are named as follows: TARGET-FRAGMENT_CHAIN (e.g. A71EV2A-x0556_0A).

In this format, the `FRAGALYSIS_DATA_DIR` will contain a target dir (the target should be named consistently throughout the files), which contains an aligned dir, which then contains individual directories for each fragment, within which we can find files for the ligand and protein.
The following are required -  a mol file for the fragment, an apo desolv file for the protein, an apo desolv file with hydrogens added for the protein, and a file containing the SMILES for the fragment.
An example structure is shown below.
```
├── TARGET
│   ├── aligned
│   │   ├── TARGET-x0001_0A
│   │   │   ├──TARGET-x0001_0A.mol
│   │   │   ├──TARGET-x0001_0A.apo-desolv.pdb
│   │   │   ├──TARGET-x0001_0A.apo-desolv-Hs.pdb
│   │   │   ├──TARGET-x0001_0A_smiles.txt
│   │   ├── TARGET-x0002_0A
│   │   │   ├──TARGET-x0002_0A.mol
│   │   │   ├──TARGET-x0002_0A.apo-desolv.pdb
│   │   │   ├──TARGET-x0002_0A.apo-desolv-Hs.pdb
│   │   │   ├──TARGET-x0002_0A_smiles.txt
```
