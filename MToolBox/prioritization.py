#!/usr/bin/env python

"""
Written by Maria Angela Diroma and Claudia Calabrese 2014
"""

import sys

from collections import Iterable
from operator import itemgetter

raw_file=open(sys.argv[1],"r").readlines()
file=[line.strip('\n').split('\t') for line in raw_file]
output=open("prioritized_variants.txt", "w")

header=['Variant Allele','Sample','Locus','Nt Variability','Codon Position','Aa Change','Aa Variability','tRNA Annotation','Disease Score','RNA predictions','Mitomap Associated Disease(s)','Mitomap Homoplasmy','Mitomap Heteroplasmy','Somatic Mutations','SM Homoplasmy','SM Heteroplasmy','ClinVar','OMIM link', 'dbSNP ID', 'Mamit-tRNA link','PhastCons20Way','PhyloP20Way','AC/AN 1000 Genomes','1000 Genomes Homoplasmy','1000 Genomes Heteroplasmy\n']
header='\t'.join(header)
output.write(header)

#dictionary with key=variant, values=annotations 
clean={}
for x in file:
	if x[1] not in clean:
		clean[x[1]]=x[1:]
	
#dictionary with key=variant, values=sample id
dic_var_sam={}
for item in file:
	key = item[1]
	dic_var_sam.setdefault(key, []).append(item[0])

#remove possible duplicates of sample IDs (if best haplogroup is not unique)
dic_var_sam2={}
for k,v in dic_var_sam.items():
	dic_var_sam2[k]=list(set(v))

prio=[]
for k in dic_var_sam2:
	samples=dic_var_sam2[k]
	tmp=[k,samples]
	tmp.extend(clean[k][1:])
	tmp[3]=float(tmp[3])
	prio.append(tmp)

#sorting (ascending nt variability) of all the variants
sort_prio=sorted(prio, key=itemgetter(3))

#unique variants in output
for x in sort_prio:
	samples=",".join(x[1])
	x[3]=str(x[3])
	annot=str("\t".join(x[2:]))
	output.write(x[0]+"\t"+samples+"\t"+annot+"\n")
	

output.close()
