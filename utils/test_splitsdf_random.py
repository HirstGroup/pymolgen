for j in range(10):

	inchi_list = []

	for i in range(10):

		n = j * 10 + i

		with open('pepe%s.inchi' %n) as infile:
			for line in infile:
				inchi_list.append(line.strip('\n'))

	print(len(inchi_list))

	print(len(set(inchi_list)))