#!/software/bin/python2.7
from Bio.SeqIO import parse
from Bio import pairwise2
from Bio.Align.Applications import MuscleCommandline

# all CDS of mercurialis 
# /scratch/cluster/monthly/gcossard/Hydrexpr_Kallisto/Hydrexpr_reads/CDS_listofbams2.txt.fas

# /scratch/cluster/monthly/gcossard/SNPcall_EXPR/AllCDSMannua_v_AllCDSRicinus.txt
# sortie de BLAST: /scratch/cluster/monthly/gcossard/SNPcall_EXPR/ALLCDS_v14.fasta sur /scratch/cluster/monthly/gcossard/Ricinus_data/TIGR_castorWGS_release_0.1.cds.fsa

blastFile = "/scratch/cluster/monthly/gcossard/SNPcall_EXPR/AllCDSMannua_v_AllCDSRicinus.txt"
ricinusFile = "/scratch/cluster/monthly/gcossard/Ricinus_data/TIGR_castorWGS_release_0.1.cds.fsa"
mercuFile = "/scratch/cluster/monthly/croux/guillaume/ALLCDS_v14.fasta"

ricinus = {}
infile = parse(ricinusFile, "fasta")
for i in infile:
	ricinus[i.id] = i.seq
infile.close()

mercu = {}
infile = parse(mercuFile, "fasta")
for i in infile:
	gene = i.id
	if gene not in mercu:
		mercu[gene] = {}
	mercu[gene] = i.seq
infile.close()

blast = {}
infile = open(blastFile, "r")
for i in infile:
	if i[0] != "#":
		i = i.strip().split("\t")
		scaffoldMercu = i[0]
		scaffoldRicinus = i[1]
		startMercu = int(i[3]) - 1
		stopMercu = int(i[4]) - 1
		startRicinus = int(i[5]) - 1
		stopRicinus = int(i[6]) - 1
		if scaffoldMercu not in blast:
			blast[scaffoldMercu] = {}
			blast[scaffoldMercu]["mercu"] = mercu[scaffoldMercu][startMercu:stopMercu]
			blast[scaffoldMercu]["ricinus"] = ""
		blast[scaffoldMercu]["ricinus"] += ricinus[scaffoldRicinus][startRicinus:stopRicinus]
infile.close()

# for pairwise
#res = ""
#for i in blast.keys()[0:3]:
#	alignement = pairwise2.align.localms(blast[i]["mercu"], blast[i]["ricinus"], 3, -1, -2, -2, one_alignment_only=1)
#	keep = []
#	for j in range(len(alignement[0][0])):
#		if alignement[0][0][j] != "-":
#			keep.append(j)
#	mercuTmp = ""
#	ricinusTmp = ""
#	for j in keep:
#		mercuTmp += alignement[0][0][j]
#		ricinusTmp += alignement[0][1][j]
#	ricinusTmp = ricinusTmp.replace("-", "N")
#	res = res + ">{0}\n{1}\n>{2}\n{3}\n".format(i + "|mercurialis|referenceMercu|Allele_1", mercuTmp, i + "|ricinus|referenceRicin|Allele_1", ricinusTmp) # g22102|mercurialis|accepted_H9|Allele_1
#	outfile = "/scratch/cluster/monthly/croux/guillaume/alignedCDS_mercu_ricinus.fas"
#	output = open(outfile, "w")
#	output.write(res)
#	output.close()

for i in blast:
	if len(blast[i]["mercu"]) != 0 and len(blast[i]["ricinus"]) != 0:
		tmp = ">{0}\n{1}\n>{2}\n{3}\n".format(i + "|mercurialis|referenceMercu|Allele_1", blast[i]["mercu"], i + "|ricinus|referenceRicin|Allele_1", blast[i]["ricinus"])
		output = open("tmp.fasta", "w")
		output.write(tmp)
		output.close()
		outfile = "/scratch/cluster/monthly/croux/guillaume/alignements/{0}_aln.fas".format(i)
		align = MuscleCommandline(input="tmp.fasta", out=outfile)
		align()
		input = parse(outfile, "fasta")
		for j in input:
			if "mercurialis" in j.id:
				nameMercu = j.id
				seqMercu = j.seq
			if "ricinus" in j.id:
				nameRicinus = j.id
				seqRicinus = j.seq
		input.close()
		keep = []
		for j in range(len(seqMercu)):
			if seqMercu[j] != "-":
				keep.append(j)
		seqMercuCorrected = ""
		seqRicinusCorrected = ""
		for j in keep:
			seqMercuCorrected += seqMercu[j]
			seqRicinusCorrected += seqRicinus[j]
		seqRicinusCorrected = seqRicinusCorrected.replace("-", "N")
		tmp = ">{0}\n{1}\n>{2}\n{3}\n".format(i + "|mercurialis|referenceMercu|Allele_1", seqMercuCorrected, i + "|ricinus|referenceRicin|Allele_1", seqRicinusCorrected)
		output = open(outfile, "w")
		output.write(tmp)
		output.close()


