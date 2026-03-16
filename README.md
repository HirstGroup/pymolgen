PyMolGen
=================================================================================================================

Installation
-------------

1. Install [Conda](https://conda.io/projects/conda/en/latest/index.html)
2. Clone this Git repository
3. Open a shell, and go to the repository and create the Conda environment:
   
        $ conda env create -f env_pymolgen.yml

4. Optionally install [openeye][https://docs.eyesopen.com/toolkits/python/quickstart-python/install.html] python toolkits and [Lilly's MedChem Rules][https://github.com/IanAWatson/Lilly-Medchem-Rules/] to filter molecules.

5. Activate the environment:
   
        $ conda activate pymolgen

6. Use the package.

Examples
---------

Example runs, please remove the -v or --verbose option for long runs

To run a molecular generation with fragment_builder.py:

python ~/pymolgen/fragment_builder.py --fragments_sdf ~/pymolgen/datasets/database1000/fragments1000.sdf --fragments_txt ~/pymolgen/datasets/database1000/fragments1000.txt --frequencies_txt ~/pymolgen/datasets/database1000/frequencies1000.txt  --parent_file ~/pymolgen/datasets/database1000/phenylisoxazole.sdf --parent_fragment_file_list ~/pymolgen/datasets/database1000/benzene.sdf ~/pymolgen/datasets/database1000/benzene.sdf --parent_mapping_1 16 0 17 0  --remove_hydrogens 21 22 --remove_hydrogens_parent_fragment 11 11 --fragments_used_file fragments_used.txt --intermediates --mw_check -o output.inchi --cap --unique --candidate_file candidates.txt -n 100 --verbose

To generate a fragment database from an SDF file of molecules:

python ~/pymolgen/fragment_mol.py -i ~/pymolgen/datasets/database1000/database1000.sdf -o 1000 -v
