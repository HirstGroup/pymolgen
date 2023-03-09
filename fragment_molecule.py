import argparse

from pymolgen.fragment_graph import *
from pymolgen.fragment_builder import get_fragment_database, get_frag_mapping, get_bond_frequencies, update_bond_frequencies
from pymolgen.fragment_mol import get_canonical_mapping
from pymolgen.molecule_formats import *

class FragmentMolecule:

    def __init__(self):
        self._graph = FragmentGraph()

    def add_fragment(self, frag_id: int, attachment_point_list) -> int:
        node_id = len(self._graph.fragments)
        self._graph.add_fragment(node_id, attachment_point_list)
        self._graph.add_node_attribute(node_id, "frag_id", frag_id)
        return node_id

    def add_bond(self, fragment_from: int, fragment_to: int, attach_from: int, attach_to: int):
        self._graph.add_bond(fragment_from, fragment_to, attach_from, attach_to)

    def list_bonds(self):
        return self._graph.bonds

    def list_free_valence_points(self):
        return self._graph.free_valence_points

    def get_frag_id(self, node_id):
        return self._graph.fragments[node_id].get_attribute('frag_id')

    def list_frag_id(self):
        frag_id_list = []
        for node_id in range(len(self._graph.fragments)):
            frag_id_list.append(self._graph.fragments[node_id].get_attribute('frag_id'))
        return frag_id_list

    def get_canonical_mapping(self, frag_id, fragment_database):
        return get_canonical_mapping(self._graph.fragments[frag_id].get_molecule(fragment_database).graph)

    def __str__(self):
        out = ''
        for i in range(len(self._graph.fragments)):
            if len(out) > 0:
                out += '-'
            out += str(self._graph.fragments[i].get_attribute('frag_id'))
        return out

def convert_fragment_molecule_to_mol(FragmentMolecule, fragment_database):

    mol = Molecule()

    frag_mol_list = []

    for i in FragmentMolecule.list_frag_id():
        frag_mol_list.append(fragment_database[i])      

    new_frag_bond_list = []

    frag_len_list = [len(i.graph.nodes) for i in frag_mol_list]

    added_frag_len_list = [0]

    for i in range(1,len(frag_len_list)):
        added_frag_len_list.append(sum(frag_len_list[:i]))

    for bond in FragmentMolecule.list_bonds():
        i = bond[0]
        j = bond[1]
        k = bond[2]
        l = bond[3]

        k += added_frag_len_list[i]
        l += added_frag_len_list[j]

        new_frag_bond_list.append((i,j,k,l))

    graphs = [x.graph for x in frag_mol_list]

    mol.graph = copy.deepcopy(networkx.disjoint_union_all(graphs))

    for bond in new_frag_bond_list:
        k = bond[2]
        l = bond[3]
        mol.graph.add_edge(k, l, order=1)        

    return mol