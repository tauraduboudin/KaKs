#!/software/bin/python2.7
from Bio.SeqIO import parse
from os import system

infile = "/scratch/cluster/monthly/croux/guillaume/alignements/res.txt"

input = open(infile, "r")

alignements = {}

for i in input:
	i = i.strip().split("\t")
	tmp = i[0].split("|")
	gene = tmp[0]
	species = tmp[1]
	if gene not in alignements:
		alignements[gene] = {}
	alignements[gene][species] = [ int(j) for j in i[2:5] ] # [ i[2], i[3], i[4] ]
input.close()

"/scratch/cluster/monthly/croux/guillaume/alignements/cleaned/"

cnt = 0
for i in alignements:
	cnt += 1
	if cnt%1000 == 0:
		print(cnt)
	infile = "/scratch/cluster/monthly/croux/guillaume/alignements/" + i + "_aln.fas"
	outfile = "/scratch/cluster/monthly/croux/guillaume/alignements/cleaned/" + i + "_aln.fas"
	if alignements[i]['mercurialis'][0] == 0 and alignements[i]['ricinus'][0] == 0:
		cmd = "cp " + infile + " " + outfile 
		system(cmd)
	else:
		test = False
		if alignements[i]['mercurialis'][1] == 0 and alignements[i]['ricinus'][1] == 0:
			phase = 1
			test = True
		if alignements[i]['mercurialis'][2] == 0 and alignements[i]['ricinus'][2] == 0 and test == False:
			phase = 2
			test = True
		if test == True:
			input = parse(infile, "fasta")
			res = ""
			for j in input:
				res = res + ">{0}\n{1}\n".format(j.id, j.seq[phase:])
			output = open(outfile, "w")
			output.write(res)
			output.close()

