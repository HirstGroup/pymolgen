import argparse
import networkx

from networkx.algorithms import isomorphism

from pymolgen.fragment_graph import *
from pymolgen.fragment_builder import get_fragment_database, get_frag_mapping, get_bond_frequencies, update_bond_frequencies
from pymolgen.fragment_mol import get_canonical_mapping
from pymolgen.molecule_formats import *

class FragmentMolecule:

    def __init__(self, build_probability=None):
        self._graph = FragmentGraph(build_probability)

    def __hash__(self):

        # do not cap before converting to networkx since free_valence_points will be stored as attribute
        #f = self.cap()

        g = self._graph.convert_to_networkx()

        return int(networkx.weisfeiler_lehman_graph_hash(g, node_attr='frag_id', edge_attr='atoms'), 16)

    def __str__(self):
        out = ''
        for i in range(len(self._graph.fragments)):
            if len(out) > 0:
                out += '-'
            out += str(self._graph.fragments[i].get_attribute('frag_id'))
        out += ':'
        for i in range(len(self._graph.bonds)):
            if i > 0:
                out += ';'
            out += f'{self._graph.bonds[i]}'
        out += f':{self.get_build_probability()}'
        return out

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        g1 = self.convert_to_networkx()
        g2 = other.convert_to_networkx()

        gm = isomorphism.GraphMatcher(g1, g2, node_match=lambda n1,n2:n1['frag_id']==n2['frag_id'], edge_match= lambda e1,e2: e1['atoms'] == e2['atoms'])

        return gm.is_isomorphic()

    def __lt__(self, other):
        if abs(self.get_build_probability() - other.get_build_probability()) < 1e-10: 
            return self.__hash__() < other.__hash__()
        else:
            return self.get_build_probability() < other.get_build_probability()

    def add_fragment(self, frag_id: int, attachment_point_list, canonical_mapping=None) -> int:
        node_id = len(self._graph.fragments)
        self._graph.add_fragment(node_id, attachment_point_list, canonical_mapping)
        self._graph.add_node_attribute(node_id, "frag_id", frag_id)
        return node_id

    def add_bond(self, fragment_from: int, fragment_to: int, attach_from: int, attach_to: int, attachment_probability: float = None):
        self._graph.add_bond(fragment_from, fragment_to, attach_from, attach_to, attachment_probability)

    def list_bonds(self):
        return self._graph.bonds

    def list_free_valence_points(self):
        return self._graph.free_valence_points

    def get_total_free_valence(self):
        total_free_valence = 0
        for i in self.list_free_valence_points():
            total_free_valence += len(i)
        return total_free_valence

    def get_frag_id(self, node_id):
        return self._graph.fragments[node_id].get_attribute('frag_id')

    def list_frag_id(self):
        frag_id_list = []
        for node_id in range(len(self._graph.fragments)):
            frag_id_list.append(self._graph.fragments[node_id].get_attribute('frag_id'))
        return frag_id_list

    def get_canonical_mapping(self, frag_id):
        return self._graph.fragments[frag_id].get_canonical_mapping()

    def get_build_probability(self):
        return self._graph.build_probability

    def cap(self):

        f = copy.deepcopy(self)

        free_valence_list = f.list_free_valence_points()

        for i in range(len(free_valence_list)):
            for j in free_valence_list[i]:
                id = f.add_fragment(-1, [0], {0:0})
                f.add_bond(i, id, j, 0)

        return f


    def convert_to_networkx(self):

        f = self.cap()

        return f._graph.convert_to_networkx()

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

    mol.graph = networkx.disjoint_union_all(graphs)

    for bond in new_frag_bond_list:
        k = bond[2]
        l = bond[3]
        mol.graph.add_edge(k, l, order=1)        

    return mol