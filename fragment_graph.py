import copy
import networkx

from typing import Tuple, Dict, List

from pymolgen.fragment_mol import get_canonical_mapping
from pymolgen.molecule import Molecule

from functools import partial
print = partial(print, flush=True)

class FragmentGraphNode:

    def __init__(self, attachment_points: List[int]):
        self._attachment_points = list(attachment_points)
        self._attributes = dict()
        self._molecule = None
        self._canonical_mapping = None

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

    def set_canonical_mapping(self, fragment_database):
        if self._canonical_mapping is None:
            self._canonical_mapping = get_canonical_mapping(fragment_database[self._attributes['frag_id']].graph)
        return self._canonical_mapping

    def manual_canonical_mapping(self, dict):
        self._canonical_mapping = dict
        return self._canonical_mapping

    def get_canonical_mapping(self):
        return self._canonical_mapping

class FragmentGraph:

    def __init__(self, build_probability=None):
        self._fragments: Dict[int, FragmentGraphNode] = dict()
        self._bonds: List[Tuple(int, int, int, int)] = []
        self._attachment_point_list = []
        self._free_valence_points = []
        if build_probability is not None:
            self._build_probability = build_probability
        else:
            self._build_probability = 1.0

    def __len__(self):
        return len(self._fragments)

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

    @property
    def build_probability(self):
        return self._build_probability

    def add_fragment(self, id: int, attachment_points: List[int], canonical_mapping=None):
        self._fragments[id] = FragmentGraphNode(attachment_points)
        self._fragments[id].manual_canonical_mapping(canonical_mapping)
        self._attachment_point_list.append(attachment_points)
        self._free_valence_points.append(attachment_points)

    def add_bond(self, fragment_from: int, fragment_to: int, attach_from: int, attach_to: int, attachment_probability: float = None):

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

        # Make bond
        self._bonds.append((fragment_from, fragment_to, attach_from, attach_to))
        self._free_valence_points[fragment_from].remove(attach_from)
        self._free_valence_points[fragment_to].remove(attach_to)
        if attachment_probability is not None:
            self._build_probability *= attachment_probability

    def add_node_attribute(self, node_id, atribute_name, atribute_value):
        self.fragments[node_id].set_attribute(atribute_name, atribute_value)


    def convert_to_networkx(self):
        g = networkx.Graph()

        canonicalise = True

        for i in range(len(self.fragments)):

            g.add_node(i, frag_id=self.fragments[i].get_attribute('frag_id'))
            if self.fragments[i].get_attribute('frag_id') != -1 and self.fragments[i].get_canonical_mapping() is None:
                canonicalise = False

        for bond in self.bonds:
            i = bond[0]
            j = bond[1]
            k = bond[2]
            l = bond[3]

            left_id, right_id = self.fragments[i].get_attribute('frag_id'), self.fragments[j].get_attribute('frag_id')

            if canonicalise is True:
                # convert to canonical atoms if atom does not have valence greater than 1
                # this is done to convert aromatic rings to canonical atoms but not saturated rings
                if left_id != -1 and self.fragments[i].attachment_points.count(k) == 1:
                    k = self.fragments[i].get_canonical_mapping()[k]
                if right_id != -1 and self.fragments[j].attachment_points.count(l) == 1:
                    l = self.fragments[j].get_canonical_mapping()[l]
            else:
                print('Canonicalise is False')

            if left_id < right_id:
                atoms_attr = f'{left_id}:{k}, {right_id}:{l}'
            elif left_id == right_id:
                if k < l:
                    atoms_attr = f'{left_id}:{k}, {right_id}:{l}'
                else:
                    atoms_attr = f'{left_id}:{l}, {right_id}:{k}'
            else:
                atoms_attr = f'{right_id}:{l}, {left_id}:{k}'

            g.add_edge(i, j, atoms=atoms_attr)


        bond_atoms = networkx.get_edge_attributes(g, 'atoms')

        return g


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

def convert_fragment_database_to_graph(fragment_database):

    print('Converting fragment database to graph ...')

    f = FragmentGraph()
    for x in range(len(fragment_database)):
        if x % 100 == 0:
            print('%s' %x, end = ' ')
        f.add_fragment(x, fragment_database[x].free_valence_list)
        f.fragments[x].set_attribute('frag_id', x)
        f.fragments[x].set_canonical_mapping(fragment_database)

    print('Converting fragment database to graph FINISHED')
    return f