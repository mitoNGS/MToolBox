#!/usr/bin/env python

"""
Written by Maria Angela Diroma 2014
"""

import sys

output_file = open("summary_tmp.txt", "w")

header1=['Sample','mtDNA Coverage','Per base depth','Best predicted haplogroup(s)','N. of homoplasmic variants','N. of heteroplasmic variants HF>=thrsld','N. of heteroplasmic variants HF<thrsld','N. of variants','N. of prioritized variants\n']
header2=['Sample','Best predicted haplogroup(s)','N. of variants','N. of prioritized variants\n']
coverage = None
heteroplasmy = None

try:
	coverage=open(sys.argv[1], "r").readlines()
	heteroplasmy=open(sys.argv[2], "r").readlines()
except:
	pass


haplo=open("mt_classification_best_results.csv", "r").readlines()
haplo=[line.strip('\n').split(',') for line in haplo]
#variants=open("prioritized_variants.txt", "r").read()
var_number=open("variant_number.txt", "r").readlines()
var_number=[line.strip('\n').split('\t') for line in var_number]

dic_haplo={}
#key=sampleID, value=best haplogroup(s)
for i in haplo:
	dic_haplo[i[0]]=i[1]

del dic_haplo['SampleID']	

dic_var={}
#key=sampleID, value=number of variants
for i in var_number:
	dic_var[i[0]]=i[1]

dic_prio={}
#key=sampleID, value=number of prioritized variants
for i in var_number:
	dic_prio[i[0]]=i[2]

if coverage!=None:
	header1='\t'.join(header1)
	output_file.write(header1)		
	coverage=[line.strip('\n').split(' ') for line in coverage]
	heteroplasmy=[line.strip('\n').split(' ') for line in heteroplasmy]
	dic_cov={}
	#key=sampleID, value=mtDNA coverage
	for i in coverage:
		dic_cov[i[1]]=i[7]
	
	dic_dp={}
	#key=sampleID, value=per base depth
	for i in coverage:
		dic_dp[i[1]]=i[12]
	
	dic_homo={}
	#key=sampleID, value=n. of homoplasmic variants
	for i in heteroplasmy:
		dic_homo[i[0]]=i[1]
		
	dic_low_hetero={}
	#key=sampleID, value=n. of low heteroplasmic variants
	for i in heteroplasmy:
		dic_low_hetero[i[0]]=i[2]

	dic_high_hetero={}
	#key=sampleID, value=n. of high heteroplasmic variants
	for i in heteroplasmy:
		dic_high_hetero[i[0]]=i[3]		
		
	for k in dic_dp:
		dpt=dic_dp[k]
		dpt=str(dpt[:-1])
		output_file.write(str(k)+"\t"+str(dic_cov[k])+"\t"+str(dpt)+"\t"+str(dic_haplo[k])+"\t"+str(dic_homo[k])+"\t"+str(dic_low_hetero[k])+"\t"+str(dic_high_hetero[k])+"\t"+str(dic_var[k])+"\t"+str(dic_prio[k])+"\n")

else:
	header2='\t'.join(header2)
	output_file.write(header2)
	for i in dic_haplo:
		haplogroup=dic_haplo[i]
		output_file.write(str(i)+"\t"+str(haplogroup)+"\t"+str(dic_var[i])+"\t"+str(dic_prio[i])+"\n")

	
#output_file.write("\n"+"=============================="+"\n\n"+"Prioritized variants"+"\n\n"+variants)

output_file.close()

