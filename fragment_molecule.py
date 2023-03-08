import argparse

from pymolgen.fragment_graph import *
from pymolgen.fragment_builder import get_fragment_database, get_frag_mapping, get_bond_frequencies, update_bond_frequencies
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

    def list_free_valence_points(self):
        return self._graph.free_valence_points

    def __str__(self):
        out = ''
        for i in range(len(self._graph.fragments)):
            if len(out) > 0:
                out += '-'
            out += str(self._graph.fragments[i].get_attribute('frag_id'))
        return out
