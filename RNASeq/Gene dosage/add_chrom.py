import sys

g2c = {}

with open("ensembl112mart.txt", "r") as infile:
	infile.readline()
	for line in infile:
		if line.strip() == "": continue
		words = line.strip().split("\t")
		g2c[words[0]] = words[-1]

with open(sys.argv[1], "r") as infile:
	print(infile.readline().strip()+" Chrom")
	for line in infile:
		line = line.strip()
		gene = line.split(" ")[0].replace('"','')
		if gene in g2c.keys():
			chrom = g2c[gene]
		else:
			chrom = "NA"
		print(line+" "+chrom)
