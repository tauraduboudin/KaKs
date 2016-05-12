#!/usr/bin/env Rscript
# measure dN and dS at 15330 alignements between annua and ricinus
library(seqinr)
setwd("/scratch/cluster/monthly/croux/guillaume/alignements/cleaned")
listOfGenes = dir(pattern="fas")

res = NULL

for(i in listOfGenes){
	x = read.alignment(i, "fasta")
	L = nchar(x$seq[[1]][1])
	for(j in 0:2 ){
		if((L - j)%%3 == 0){
			L = L - j
			break
		}
	}
	x$seq[[1]][1] = substr(x$seq[[1]][1], 1, L)
	x$seq[[2]][1] = substr(x$seq[[2]][1], 1, L)
	y = kaks(x)
	res = rbind(res, c(strsplit(i, "_")[[1]][1], L, y$ka, y$ks))
}

colnames(res) = c("gene", "L", "dN", "dS")
write.table(res, "/scratch/cluster/monthly/croux/guillaume/alignements/cleaned/dNdS_annua_ricinus.txt", col.names=T, row.names=F, quote=F)

