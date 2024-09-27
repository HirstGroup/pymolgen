from pymolgen.fragment_builder import get_fragment_database
from pymolgen.fragment_graph import convert_fragment_database_to_graph

class FragmentDatabaseEquivalent:

	def __init__(self, fragment_database_graph):

		self.fragments = dict()
		self.equivalent_fragments = dict()
		self.hydrogenated_fragments = dict()

		for frag_id, fragment in fragment_database_graph.fragments.items():

			self.fragments[frag_id] = fragment

	def generate_equivalent_fragments(self):

		for frag_id, fragment in self.fragments:

			



