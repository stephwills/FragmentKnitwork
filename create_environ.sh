conda create --name fragnet
conda install -c conda-forge rdkit -y
conda install -c https://username:password@west.rosettacommons.org/pyrosetta/conda/release/ pyrosetta -y
pip install tqdm joblib fragmenstein MDanalysis prolif
git clone https://github.com/schrodinger/pymol-open-source.git && cd pymol-open-source && python setup.py install && cd .. && rm -r pymol-open-source