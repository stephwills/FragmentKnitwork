
### Fragalysis formatted files
This assumes you have a `FRAGALYSIS_DATA_DIR` which contains target directories that have been downloaded from Fragalysis. This adopts the old style of formatting that was used in Fragalysis, whereby fragments are named as follows: TARGET-FRAGMENT_CHAIN (e.g. A71EV2A-x0556_0A).


With the fragalysis update the file format has changed, so indicate what version you are using in utils/Config.py.

### v1
In this format, the `FRAGALYSIS_DATA_DIR` will contain a target dir (the target should be named consistently throughout the files), which contains an aligned dir, which then contains individual directories for each fragment, within which we can find files for the ligand and protein.
The following are required -  a mol file for the fragment, an apo desolv file for the protein, an apo desolv file with hydrogens added for the protein, and a file containing the SMILES for the fragment.
An example structure is shown below.
```
├── TARGET
│   ├── aligned
│   │   ├── TARGET-x0001_0A
│   │   │   ├──TARGET-x0001_0A.mol
│   │   │   ├──TARGET-x0001_0A_apo-desolv.pdb
│   │   │   ├──TARGET-x0001_0A_apo-desolv-Hs.pdb
│   │   │   ├──TARGET-x0001_0A_smiles.txt
│   │   ├── TARGET-x0002_0A
│   │   │   ├──TARGET-x0002_0A.mol
│   │   │   ├──TARGET-x0002_0A_apo-desolv.pdb
│   │   │   ├──TARGET-x0002_0A_apo-desolv-Hs.pdb
│   │   │   ├──TARGET-x0002_0A_smiles.txt
```
### v2
The fragments and files are named differently in version 2. Ax00001A represents the fragment name, which is used throughout all file names.
The aligned dir is now called aligned_files.
```
├── TARGET
│   ├── aligned_files
│   │   ├── Ax0001a
│   │   │   ├──Ax0001a_ligand.mol
│   │   │   ├──Ax0001a_apo-desolv.pdb
│   │   │   ├──Ax0001a_apo-desolv-Hs.pdb
│   │   │   ├──Ax0001a_ligand.smi
│   │   ├── TARGET-x0002_0A
│   │   │   ├──Ax0001a.mol
│   │   │   ├──Ax0001a_apo-desolv.pdb
│   │   │   ├──Ax0001a_apo-desolv-Hs.pdb
│   │   │   ├──Ax0001a_ligand.smi
```
