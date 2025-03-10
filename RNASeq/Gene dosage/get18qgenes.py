genes = []

with open("mart_export-18q22-tergenes.txt", "r") as infile:
	infile.readline()
	for line in infile:
		genes.append(line.split("\t")[0])

with open("chr18.txt", "r") as infile:
	print(infile.readline())
	for line in infile:
		if line.split(" ")[0].replace('"','') in genes:
			print(line.strip())
