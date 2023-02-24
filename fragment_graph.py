from typing import Tuple, Dict, List

class FragmentGraphNode:

	def __init__(self, attachment_points: List[int]):
		self._attachment_points = attachment_points

	@property
	def attachment_points(self):
		return list(self._attachment_points)

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

		# Make bond
		self._bonds.append((fragment_from, fragment_to, attach_from, attach_to))
		self._free_valence_points[fragment_from].remove(attach_from)
		self._free_valence_points[fragment_to].remove(attach_to)
