from pymolgen.fragment_molecule import *

f = FragmentMolecule()

f.add_fragment(10, [0,1]) # 0
f.add_fragment(20, [2,2]) # 1
f.add_fragment(20, [2,2]) # 2
f.add_fragment(30, [4])   # 3
f.add_fragment(40, [5])   # 4

f.add_bond(0, 1, 0, 2)
f.add_bond(1, 2, 2, 2)
f.add_bond(2, 3, 2, 4)
f.add_bond(0, 4, 1, 5)

print(f.list_free_valence_points())
