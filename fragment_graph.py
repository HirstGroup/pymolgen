class FragmentGraph:

	def __init__(self):

	    self.fragments = [-1] #list of fragment indexes from database, parent fragment shown as -1
	    self.frag_edge_list = [] #list of bonds between fragments
	    self.frag_free_valence_list = [[]] #list of free valence points of each fragment		

	def add_fragment(self, fragment):

		self.fragments.append(fragment)

		pass

