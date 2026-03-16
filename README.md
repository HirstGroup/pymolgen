PyMolGen
=================================================================================================================

Installation
-------------

1. Install [Conda](https://conda.io/projects/conda/en/latest/index.html)
2. Clone this Git repository
3. Open a shell, and go to the repository and create the Conda environment:
   
        $ conda env create -f env_pymolgen.yml

4. Optionally install [openeye](https://docs.eyesopen.com/toolkits/python/quickstart-python/install.html) python toolkits and [Lilly's MedChem Rules](https://github.com/IanAWatson/Lilly-Medchem-Rules/) to filter molecules.

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

Tests
------

The PyMolGen package uses pytest for testing. Main tests can be found in the `test/` directory. All core molecular generation functionality of the package can be tested with the basic `env_pymolgen.yml` enviroment. Tests involving molecule filtering require the openeye and Lilly's MedChem Rules installed. Tests involving the ML scoring model of generated molecules described in the original PyMolGen publication require the installation of that ML code. 

Structure of the package and main entry points
-----------------------------------------------

The main entry points are located in the root directory. They provide the functionality for the derivation of the fragment database and the generation of new molecules based on that database.

* `fragment_mol.py`: derive the fragment combination rules from a database
* `fragment_mol_combine.py`: allows for combining the fragment combination rules if performed in parts in an original database (for performance purposes)
* `fragment_builder.py`: generate new molecules from a fragment database, using the all-atom description of fragments
* `fragment_molecule_builder.py`: generate new molecules from a fragment database, using the fragment graph description of fragments, allowing for faster performance and the use of build probabilities for random and systematic generation

`datasets/` contains the fragment database derived from ChEMBL, `analysis/` and `utils/` provide useful scripts for the analysis and manipulation of fragment databases, with tests located in their corresponding `test/` subdirectories. Main tests are located in `test/`. 