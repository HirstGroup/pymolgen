import copy
import networkx

from typing import Tuple, Dict, List

from pymolgen.fragment_mol import get_canonical_mapping
from pymolgen.molecule import Molecule


class FragmentGraphNode:

    def __init__(self, attachment_points: List[int]):
        self._attachment_points = attachment_points
        self._attributes = dict()
        self._molecule = None

    @property
    def attachment_points(self):
        return list(self._attachment_points)

    def get_molecule(self, fragment_database):
        if self._molecule is None:
            self._molecule = fragment_database[self._attributes['frag_id']]
        return self._molecule

    def set_attribute(self, key: str, val):
        self._attributes[key] = val

    def get_attribute(self, key: str):
        return self._attributes[key]


class FragmentGraph:

    def __init__(self):
        self._fragments: Dict[int, FragmentGraphNode] = dict()
        self._bonds: List[Tuple(int, int, int, int)] = []
        self._attachment_point_list = []
        self._free_valence_points = []

    @property
    def fragments(self):
        return dict(self._fragments)

    @property
    def bonds(self):
        return list(self._bonds)

    @property
    def attachment_point_list(self):
        return list(self._attachment_point_list)

    @property
    def free_valence_points(self):
        return list(self._free_valence_points)

    def add_fragment(self, id: int, attachment_points: List[int]):
        if len(attachment_points) < 1:
            raise ValueError("A fragment must have at least 1 attachment point")
        self._fragments[id] = FragmentGraphNode(attachment_points)
        self._attachment_point_list.append(attachment_points)
        self._free_valence_points.append(attachment_points)

    def add_bond(self, fragment_from: int, fragment_to: int, attach_from: int, attach_to: int):

        if fragment_from > fragment_to:
            # Ensure bonds are always stored 
            # in acending-fragment order
            tmp = fragment_to
            fragment_to = fragment_from
            fragment_from = tmp
            tmp = attach_to
            attach_to = attach_from
            attach_from = tmp

        # Check acending order and that fragments don't bond to themselves
        assert fragment_from < fragment_to

        # Check bond is between existing fragments
        assert 0 <= fragment_from < len(self._fragments)
        assert 0 <= fragment_to < len(self._fragments)

        # Check attachment points are valid
        assert attach_from in self._fragments[fragment_from].attachment_points
        assert attach_to in self._fragments[fragment_to].attachment_points
        assert attach_from in self._attachment_point_list[fragment_from]
        assert attach_to in self._attachment_point_list[fragment_to]

        # Check that the attachment points are free
        assert attach_from in self._free_valence_points[fragment_from]
        assert attach_to in self._free_valence_points[fragment_to]

        # Make bon
        self._bonds.append((fragment_from, fragment_to, attach_from, attach_to))
        self._free_valence_points[fragment_from].remove(attach_from)
        self._free_valence_points[fragment_to].remove(attach_to)

    def add_node_attribute(self, node_id, atribute_name, atribute_value):
        self.fragments[node_id].set_attribute(atribute_name, atribute_value)

def convert_fragment_graph_to_mol(FragmentGraph, fragment_database):

    mol = Molecule()

    frag_mol_list = []

    for i in FragmentGraph.fragments:
        frag_mol_list.append(fragment_database[i])      

    new_frag_bond_list = []

    frag_len_list = [len(i.graph.nodes) for i in FragmentGraph.bonds]

    added_frag_len_list = [0]

    for i in range(1,len(frag_len_list)):
        added_frag_len_list.append(sum(frag_len_list[:i]))

    for bond in FragmentGraph.bonds:
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