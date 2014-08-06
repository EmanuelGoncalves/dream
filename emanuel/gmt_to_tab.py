file_name = 'go_terms/c2.all.v4.0.symbols.gmt'

with open(file_name, 'Ur') as f:
	file = open(file_name + '.tab', 'w')
	file.write('term' + '\t' + 'gene' + '\n')

	for line in f:
		values = line.split('\t')

		for i in range(2, len(values)):
			file.write(values[0] + '\t' + values[i].strip() + '\n')

	file.close()