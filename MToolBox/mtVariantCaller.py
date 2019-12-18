#!/usr/bin/env python

"""
Written by Claudia Calabrese - claudia.calabrese23@gmail.com
	and Domenico Simone - dome.simone@gmail.com
"""

import sys, os, glob, math
#print sys.version
import re
import ast
from collections import OrderedDict
import vcf

import pandas as pd
import scipy as sp
#######################################################

#defines global variables for Indels
def varnames(i):
	#global CIGAR, readNAME, seq, qs, refposleft, mate, refposright
	CIGAR=i[5]
	readNAME=i[0]
	seq=i[9]
	qs=i[10]
	refposleft=int(i[3])-1
	mate=int(i[1])
	#check strand
	if mate & 16 == 16:
		strand = '-'
	else:
		strand = '+'
	return CIGAR, readNAME, seq, qs, refposleft, strand #, refposright

#defines global variables for MT-table parsing
def varnames2(b,c,i):
	#global Position, Ref, Cons, Cov, A,C,G,T
	global Position, Ref, Cov, A,C,G,T,A_f,C_f,G_f,T_f,A_r,C_r,G_r,T_r
	Position=int((i[0]).strip())
	Ref=(i[1]).strip()
	Cov=int((i[3]).strip())
	A=b[0]
	C=b[1]
	G=b[2]
	T=b[3]
	A_f=c[0]
	C_f=c[1]
	G_f=c[2]
	T_f=c[3]
	A_r=c[4]
	C_r=c[5]
	G_r=c[6]
	T_r=c[7]
	return Position, Ref, Cov, A,C,G,T,A_f,C_f,G_f,T_f,A_r,C_r,G_r,T_r
#Heteroplasmic fraction quantification
def heteroplasmy(cov, Covbase):
	try:
		if Covbase >= cov: 
			Heteroplasmy=float(cov)/float(Covbase)
			het=round(Heteroplasmy, 3)
			return het
		else:
			return 1.0
	except ZeroDivisionError:
		het=1.0
		return het
#defines mathematical operations
#sum
def sum(left):
	s=0
	for i in left:
		s+=int(i)
	return s
#median
def median(l):
	try:
		if len(l)%2 != 0:
			median= sorted(l)[((len(l)+1)/2)-1]
		else:
			m1 = sorted(l)[(len(l)/2)+1-1]
			m2 = sorted(l)[(len(l)/2)-1]
			median= (float(m1)+float(m2))/2
		return median
	except ZeroDivisionError:
		return 0
#mean
def mean(list):
	try:
		s=sum(list)
		m=float(s)/float(len(list))
		return m
	except ZeroDivisionError:
		m=0
		return m
#defines function for value errors
def error(list):
	try:
		list.remove('')
	except ValueError:
		pass
#defines the function searching for and filtering indels within the read sequence	
def SearchINDELsintoSAM(readNAME,mate,CIGAR,seq,qs,refposleft,tail=5):
	m=re.compile(r'[a-z]', re.I)
	if 'I' in CIGAR:
		pcigar=CIGAR.partition('I')
		if 'S' in pcigar[0]:
			softclippingleft = int(re.findall(r'(\d+)S', pcigar[0])[0]) # calc softclipped bases left
		else:
			softclippingleft=0 # no softclipped bases left
		if 'S' in pcigar[-1]:
			softclippingright=int(re.findall(r'(\d+)S', pcigar[-1])[0]) # calc softclipped bases right
		else:
			softclippingright=0 # no softclipped bases right
		if 'H' in pcigar[0]:
			hardclippingleft=int(re.findall(r'(\d+)H', pcigar[0])[0]) # calc hardclipped bases left
		else:
			hardclippingleft=0 # no hardtclipped bases left
		if 'H' in pcigar[-1]:
			hardclippingright=int(re.findall(r'(\d+)H', pcigar[-1])[0]) # calc hardclipped bases right
		else:
			hardclippingright=0 # no hardclipped bases right
		left=m.split(pcigar[0])
		right=m.split(pcigar[-1])
		error(right)
		right = map(lambda x:int(x), right)
		left = map(lambda x:int(x),left)
		numIns=int(left[-1])
		left.pop(-1)
		s=sum(left) # number of 5' flanking positions to Ins
		lflank=s-(softclippingleft+hardclippingleft) # exclude softclipped bases from left pos calculation
		sr=sum(right)-(hardclippingright+softclippingright) # number of 3' flanking positions to Ins
		rLeft=refposleft+lflank # modified
		Ins=seq[s:s+numIns]
		qsInsASCI=qs[s:s+numIns]
		qsInsASCI.split()
		qsIns=[]
		type='Ins'
		if lflank >= tail and sr>= tail:
			for x in qsInsASCI:
				numeric=(ord(x)-33)
				qsIns.append(numeric)
		else:
			qsIns='delete'
		res=[]
		res.append(type)
		res.append(readNAME)
		res.append(mate)
		res.append(int(rLeft))
		res.append(Ins)
		res.append(qsIns)
		return res
	elif 'D' in CIGAR:
		pcigar=CIGAR.partition('D')
		if 'S' in pcigar[0]:
			softclippingleft=int(re.findall(r'(\d+)S', pcigar[0])[0]) # calc softclipped bases left
		else:
			softclippingleft=0 # no softclipped bases left
		if 'S' in pcigar[-1]:
			softclippingright=int(re.findall(r'(\d+)S', pcigar[-1])[0]) # calc softclipped bases right
		else:
			softclippingright=0 # no softclipped bases right
		if 'H' in pcigar[0]:
			hardclippingleft=int(re.findall(r'(\d+)H', pcigar[0])[0]) # calc hardclipped bases left
		else:
			hardclippingleft=0 # no hardtclipped bases left
		if 'H' in pcigar[-1]:
			hardclippingright=int(re.findall(r'(\d+)H', pcigar[-1])[0]) # calc hardclipped bases right
		else:
			hardclippingright=0 # no hardclipped bases right
		left=m.split(pcigar[0])
		right=m.split(pcigar[-1])
		error(right)
		right = map(lambda x:int(x),right)
		left = map(lambda x:int(x),left)
		nDel=int(left[-1])
		left.pop(-1)
		s=sum(left)
		sr=sum(right)-(hardclippingright+softclippingright) # number of 3' flanking positions to Del
		lflank=s-(softclippingleft+hardclippingleft) # exclude 5' flanking softclipped/hardclipped bases from left pos calculation
		rLeft=(refposleft+lflank)
		lowlimit=rLeft+1
		rRight=lowlimit+nDel
		Del=range(lowlimit,rRight)
		type='Del'
		qsDel=[]
		if lflank>=tail and sr>=tail:
			qsLeft=qs[(lflank-5):lflank]
			qsRight=qs[lflank:(lflank+5)]
			qsL=[]
			qsR=[]
			for x in qsLeft:
				qsL.append(ord(x)-33)
			for x in qsRight:
				qsR.append(ord(x)-33)
			medL=median(qsL)
			medR=median(qsR)
			qsDel.append(medL)
			qsDel.append(medR)
		else:
			qsDel='delete'
		res=[]
		res.append(type)
		res.append(readNAME)
		res.append(mate)
		res.append(int(rLeft))
		res.append(Del)
		res.append(qsDel)
		return res
	else:
		pass 

#defines function searching for point mutations. It produces both the consensus base and variant(s) as output 
def findmutations(A, C, G, T, Position, Ref, Cov, minrd, A_f, C_f, G_f, T_f, A_r, C_r, G_r, T_r):
	d_strand={'A':str(A_f)+';'+str(A_r),'C':str(C_f)+';'+str(C_r),'G':str(G_f)+';'+str(G_r),'T':str(T_f)+';'+str(T_r)}
	oo = []
	var=[]
	bases=[]
	strand=[]
	if A>=minrd:
		var.append(A)
		bases.append('A')
	if C>=minrd:
		var.append(C)
		bases.append('C')
	if G>=minrd:
		var.append(G)
		bases.append('G')
	if T>=minrd:
		var.append(T)
		bases.append('T')
	if len(var)>=2:
		if Ref in bases:
			indexRef=bases.index(Ref)
			bases.remove(Ref)
			var.remove(var[indexRef])
			for x in bases:
				strand.append(d_strand[x])
		if len(strand)==0:
			strand=['und;und']
		o=[Position, Ref, Cov, bases, var, strand]
		return o
	elif len(var)==1 and Ref not in bases:
		strand.append(d_strand[bases[0]])
		if len(strand)==0:
			strand=['und;und']
		o=[Position, Ref, Cov, bases, var, strand]
		return o
	else:	
		return oo

#Wilson confidence interval lower bound
def CIW_LOW(het, Covbase):
	'''The function calculates the heteroplasmic fraction and the related
	confidence interval with 95% of coverage probability,
	considering a Wilson score interval when n<=40 
	CIw= [1/(1+(1/n)*z^2)] * [p + (1/2n)*z^2 +- z(1/n *(p*q) + ((1/(4n^2))*z^2))^1/2]
	'''
	p=het
	n=Covbase
	z=1.96
	q=1-het
	num=p*q
	squarez=z*z
	squaren=n*n
	wilsonci_low=round((p+(z*z)/(2*n)-z*(math.sqrt(p*q/n+(z*z)/(4*(n*n)))))/(1+z*z/n),3)
	if wilsonci_low<0.0:
		return 0.0
	else:
		return wilsonci_low
		
#Wilson confidence interval upper bound
def CIW_UP(het, Covbase):
	'''The function calculates the heteroplasmic fraction and the related
	confidence interval with 95% of coverage probability,
	considering a Wilson score interval when n<=40 
	CIw= [1/(1+(1/n)*z^2)] * [p + (1/2n)*z^2 +- z(1/n *(p*q) + ((1/(4n^2))*z^2))^1/2]
	'''
	p=het
	n=Covbase
	z=1.96
	q=1-het
	num=p*q
	squarez=z*z
	squaren=n*n
	wilsonci_up=round((p+(z*z)/(2*n)+z*(math.sqrt(p*q/n+(z*z)/(4*(n*n)))))/(1+z*z/n),3)
	if wilsonci_up>1.0:
		return 1.0
	else:
		return wilsonci_up

#Agresti-Coull confidence interval lower bound		
def CIAC_LOW(cov,Covbase):
	'''The function calculates the heteroplasmic fraction and the related confidence interval 
	for heteroplasmic fraction with 95% of coverage probability,considering the 
	Agresti-Coull interval when n>40'''
	z=1.96
	if cov > Covbase:
		Covbase = cov
	n=Covbase
	X=cov+(z*z)/2
	N=n+(z*z)
	P=X/N
	Q=1-P
	agresticoull_low=round(P-(z*(math.sqrt(P*Q/N))),3)
	if agresticoull_low<0.0:
		return 0.0
	else:
		return agresticoull_low

#Agresti-Coull confidence interval upper bound	
def CIAC_UP(cov,Covbase):
	'''The function calculates the heteroplasmic fraction and the related confidence interval 
	for heteroplasmic fraction with 95% of coverage probability,considering the 
	Agresti-Coull interval when n>40'''
	z=1.96
	if cov > Covbase:
		Covbase = cov
	n=Covbase
	X=cov+(z*z)/float(2)
	N=n+(z*z)
	P=X/N
	Q=1-P
	agresticoull_up=round(P+(z*(math.sqrt(P*Q/N))),3)
	if agresticoull_up>1.0:
		return 1.0
	else:	
		return agresticoull_up

#IUPAC dictionary
dIUPAC={'R':['A','G'],'Y':['C','T'],'S':['G','C'],'W':['A','T'],'K':['G','T'],'M':['A','C'],'B':['C','G','T'],'D':['A','G','T'],'H':['A','C','T'],'V':['A','C','G'],'N':['A','C','G','T']}
#searches for IUPAC codes and returns the ambiguity
#returns '' if nucleotide in reference is N
def getIUPAC(ref_var, dIUPAC):
	iupac_code = ['']
	for i in dIUPAC.iteritems():
		i[1].sort()
		if ref_var == i[1]:
			iupac_code= [i[0]]
	return iupac_code
			

def mtvcf_main_analysis(mtable, sam, name2, tail, Q, minrd):
	mtable=[i.split('\t') for i in mtable]
	mtable.remove(mtable[0])
	sam=sam.readlines()
	sam=[i.split('\t') for i in sam]
	#assigns a null value to fields. NB: refpos is the leftmost position of the subject seq minus 1.
	CIGAR=''
	readNAME=''
	seq=''
	qs=''
	refposleft=''
	mate=''
	#assembly mtDNA ref sequence from MT-table
	mtDNA=[]
	for i in mtable:
		mtDNA.append((i[1]).strip())
	mtDNAseq="".join(mtDNA)
	#create list of the total depth per position across the mtDNA
	Coverage=[]
	for i in mtable:
		Coverage.append((i[3]).strip())
	#-----------------------------------------------------------------------------------------
	#apply functions to sam file and write outputs into a dictionary
	dic={}
	dic['Ins']=[]
	dic['Del']=[]
        print "minimum QS is = {0}\n".format(Q)
        print "minimum read depth per position is = {0}\n".format(minrd)
	print "\n\nsearching for indels in {0}.. please wait...\n\n".format(name2)
	for i in sam:
		[CIGAR, readNAME, seq, qs, refposleft, mate] = varnames(i)
		# varnames()
		if 'I' in CIGAR or 'D' in CIGAR:
			r=SearchINDELsintoSAM(readNAME,mate,CIGAR,seq, qs,refposleft,tail=tail)
			dic[r[0]].append(r[1:])
	#############
	rposIns={}
	rposDel={}
	for i in dic['Ins']:
		if i[2] not in rposIns:
			if i[-1] == 'delete':
				pass
			else:
				rposIns[i[2]]=[]
				rposIns[i[2]].append(i[3:]+[i[1]]) #append qs, alt allele and strand
		else:
			if i[-1] == 'delete':
				pass
			else:
				rposIns[i[2]].append(i[3:]+[i[1]])

	for i in dic['Del']:
		if i[2] not in rposDel:
			if i[-1]=='delete':
				pass
			else:
				rposDel[i[2]]=[]			
				rposDel[i[2]].append(i[3:]+[i[1]])
		else:
			if i[-1]=='delete':
				pass
			else:
				rposDel[i[2]].append(i[3:]+[i[1]])
	############
	dicqsDel={}
	dicqsIns={}
	#########
	for i in rposIns:
		dicqsIns[i]=[]
		for x in rposIns.get(i):
			for j in range(len(x[-2])):
                                if int(x[-2][j])>=Q:
					pass
				else:
					x[-2][j]='-'
			if '-' in x[-2]:
				pass
			else:
				dicqsIns[i].append(x)
	################
	for i in rposDel:
		dicqsDel[i]=[]
		for x in rposDel.get(i):
			for j in range(len(x[-2])):
                                if int(x[-2][j])>=Q:
					pass
				else:
					x[-2][j]='-'
			if '-' in x[-2]:
				pass
			else:
				dicqsDel[i].append(x)
	##########
	dicIns={}
	dicDel={}
	for i in dicqsIns:
		dicIns[i]=[]
		b=[]
		a=dicqsIns.get(i)
		for j in a:
			b.append(str(j[0]))
		s=set(b)
		for x in s:
			if b.count(x)>=minrd:
				for z in a:
					if x in z:
						dicIns[i].append(z)
	for i in dicqsDel:
		dicDel[i]=[]
		b=[]
		a=dicqsDel.get(i)
		for j in a:
			b.append(str(j[0]))
		s=set(b)
		for x in s:
			l=b.count(x)
			if l>=minrd:
				for z in a:
					if x==str(z[0]):
						dicDel[i].append(z)
	Final={}
	for i in dicIns:
		Final[i]=[]
		qs1=[]
		strand1=[]
		bases1=[]
		bases2=[]
		a=dicIns.get(i)
		l=len(dicIns.get(i))
		depth=[]
		if l>0:
			for x in a:
				bases2.append(x[0])
			b=set(bases2)
			for z in b:
				if z != '':
					n=bases2.count(z)
					if n>=minrd:
						qs2=[]
						strand=[]
						for x in a:
							if str(x[0])==z:
								qs2.extend(x[-2])
								strand.extend(x[-1])
						fwd = strand.count('+')
						rv = strand.count('-')
						strand2 = [str(fwd)+';'+str(rv)]
						bases1.append(z)
						qs1.append(median(qs2))
						depth.append(n)
						strand1.append(strand2)
			r=['ins', bases1, qs1, depth, strand1]
			Final[i].append(r)
	for i in dicDel:
		if i in Final:
			qs1=[]
			strand1=[]
			bases1=[]
			bases2=[]
			a=dicDel.get(i)
			l=len(a)
			depth=[]
			if l>0:
				for x in a:
					bases2.append(str(x[0]))
				b=set(bases2)
				for z in b:
					if z != '':
						n=bases2.count(z)
						if n>=minrd:
							qs2=[]
							strand=[]
							for x in a:
								if str(x[0])==z:
									qs2.extend(x[-2])
									strand.extend(x[-1])
							fwd = strand.count('+')
							rv = strand.count('-')
							strand2 = [str(fwd)+';'+str(rv)]
							strand1.append(strand2)
							bases1.append(z)
							qs1.append(median(qs2))
							depth.append(n)
				r=['del', bases1, qs1, depth, strand1]
				Final[i].append(r)
		else:
			Final[i]=[]
			qs1=[]
			strand1=[]
			bases1=[]
			bases2=[]
			a=dicDel.get(i)
			l=len(dicDel.get(i))
			depth=[]
			if l>0:
				for x in a:
					bases2.append(str(x[0]))
				b=set(bases2)
				for z in b:
					if z != '':
						n=bases2.count(z)
						if n>=minrd:
							qs2=[]
							strand=[]
							for x in a:
								if str(x[0])==z:
									qs2.extend(x[-2])
									strand.extend(x[-1])
							fwd = strand.count('+')
							rv = strand.count('-')
							strand2 = [str(fwd)+';'+str(rv)]
							strand1.append(strand2)
							bases1.append(z)
							qs1.append(median(qs2))
							depth.append(n)
				r=['del', bases1, qs1, depth, strand1]
				Final[i].append(r)
	ref=sorted(Final)
	Indels={}
	Indels[name2]=[]
	for i in ref:
		if len(Final.get(i))>0:
			for x in Final.get(i):
				if x[0]=='ins' and x[1] != []: #is not empty
					bases=x[1]
					qs=x[2]
					cov=x[3]
					strand=x[4]
					Refbase=mtDNAseq[int(i)-1]
					Variant=map(lambda x:Refbase+x,bases)
					InsCov=map(lambda x:int(x),cov)
					Covbase=int(Coverage[int(i)-1])
					if max(InsCov) > Covbase:
						Covbase = max(InsCov)
					QS=map(lambda x:round(float(x),2),qs)
					hetfreq=map(lambda x:heteroplasmy(x,Covbase),InsCov)
					if Covbase <=40:
						het_ci_low=map(lambda x: CIW_LOW(x, Covbase), hetfreq)
						het_ci_up=map(lambda x: CIW_UP(x, Covbase), hetfreq)
					else:
						het_ci_low=map(lambda x: CIAC_LOW(x,Covbase), InsCov)
						het_ci_up=map(lambda x: CIAC_UP(x,Covbase), InsCov)
					ins=[i, Refbase, Covbase, Variant, InsCov, strand, QS, hetfreq, het_ci_low, het_ci_up,'ins']
					Indels[name2].append(ins)
				else:
					if x[1] != [] :#is not empty
						Refbase=[]
						cov=x[3]
						strand=x[4]
						DelCov=map(lambda x:int(x),cov)
						qs=x[2]
						deletions=[]
						Covbase=[]
						for j in xrange(len(x[1])):
							dels=ast.literal_eval(x[1][j])
							delflank=dels[0]-2
							delfinal=dels[-1]
							covlist=Coverage[delflank:delfinal]
							convert=map(lambda x:int(x), covlist)
							Covbase.append(median(convert))
						maxcovbase=max(Covbase)
						Covbase=int(maxcovbase)
						if max(DelCov) > Covbase:
							Covbase = max(DelCov)
						hetfreq=map(lambda x:heteroplasmy(x,Covbase),DelCov)
						if Covbase <=40:
							het_ci_low=map(lambda x: CIW_LOW(x, Covbase), hetfreq)
							het_ci_up=map(lambda x: CIW_UP(x, Covbase), hetfreq)	
						else:
							het_ci_low=map(lambda x: CIAC_LOW(x,Covbase), DelCov)
							het_ci_up=map(lambda x: CIAC_UP(x,Covbase), DelCov)
						for j in xrange(len(x[1])):
							dels=ast.literal_eval(x[1][j])
							delflank=dels[0]-2
							delfinal=dels[-1]
							deletions.append(mtDNAseq[delflank])
							Refbase.append(mtDNAseq[delflank:delfinal])
						dele=[(dels[0]-1), Refbase, Covbase, deletions, DelCov, strand, qs, hetfreq, het_ci_low, het_ci_up, 'del']
						Indels[name2].append(dele)
	Subst={}
	Subst[name2] = []
	print "\n\nsearching for mismatches in {0}.. please wait...\n\n".format(name2)
	for i in mtable:
		b=ast.literal_eval((i[-2]).strip())
		c=ast.literal_eval((i[-1]).strip())
		varnames2(b,c,i)
		a=findmutations(A,C,G,T,Position, Ref, Cov,minrd,A_f,C_f,G_f,T_f,A_r,C_r,G_r,T_r)
		if len(a) > 0:
			hetfreq=map(lambda x:heteroplasmy(x,Cov),a[-2])
			if Cov<=40:
				het_ci_low=map(lambda x: CIW_LOW(x,Cov),hetfreq)
				het_ci_up=map(lambda x: CIW_UP(x,Cov),hetfreq)
			else:
				het_ci_low=map(lambda x: CIAC_LOW(x,Cov), a[-2])
				het_ci_up=map(lambda x: CIAC_UP(x,Cov), a[-2])
			a.append('PASS')
			a.append(hetfreq)
			a.append(het_ci_low)
			a.append(het_ci_up)
			a.append('mism')
			Subst[name2].append(a)
	Indels[name2].extend(Subst[name2])
	return Indels # it's a dictionary

### END OF MAIN ANALYSIS

#The dictionary with all the samples variations found
#applies the analysis only to OUT folders with OUT.sam, mt-table.txt and fasta sequence files.


def get_consensus_single(i, hf_max=0.8, hf_min=0.2):
	consensus_value=[]
	if len(i) != 0:
		for var in i:
			if var[-1] == 'mism' and max(var[7]) > hf_max:
				index=var[7].index(max(var[7]))
				basevar=var[3][index]
				res=[var[0], [basevar], 'mism']
				consensus_value.append(res)
			elif var[-1] == 'mism' and max(var[7]) >= hf_min and max(var[7]) <= hf_max:
				basevar=sp.array([var[1]]+var[3])
				#keep only basevar >= hf_min for IUPAC
				ref_hf = 1-sp.sum(var[7])
				hf_var = [ref_hf]
				hf_var.extend(var[7])
				hf_var = sp.array(hf_var)
				ii = sp.where(hf_var >= hf_min)[0]
				basevar = basevar[ii].tolist()
				basevar.sort()
				a=getIUPAC(basevar, dIUPAC)
				res=[var[0], a, 'mism']
				consensus_value.append(res)
			elif var[-1] == 'mism' and max(var[7]) < hf_min: #put the reference allele in consensus
				res=[var[0], [var[1]], 'mism']
			elif var[-1] == 'ins' and max(var[7]) > hf_max:
				index=var[7].index(max(var[7]))
				basevar=var[3][index]
				res=[var[0], [basevar], 'ins']
				consensus_value.append(res)
			elif var[-1] == 'del' and max(var[7]) > hf_max:
				index=var[7].index(max(var[7]))
				basevar=var[3][index]
				del_length=len(var[1][0]) - len(basevar)
				start_del=var[0]+1
				end_del=start_del+del_length
				res=[var[0], range(start_del,end_del), 'del']
				consensus_value.append(res)
			else:
				pass
		return consensus_value

def get_consensus(dict_of_dicts, hf_max, hf_min):
	"""Dictionary of consensus variants, for fasta sequences"""
	Consensus = {}
	for i in dict_of_dicts:
		Consensus[i] = get_consensus_single(dict_of_dicts[i],hf_max,hf_min)
	return Consensus

def VCFoutput(dict_of_dicts, reference='RSRS', name='sample'):
    print "Reference sequence used for VCF: %s" % reference
    VCF_RECORDS = []
    present_pos = set()
    # for each sample in dict_of_dicts
    for sample in dict_of_dicts.keys():
        #gets variants found per sample
        val = dict_of_dicts[sample]
        for variant in val:
            #variant[5] is SDP
            # if the v. position was never encountered before, is heteroplasmic and is a deletion
            if variant[0] not in present_pos and max(variant[7])<1 and variant[-1]=='del':
                allelecount=[1]*len(variant[1])
                aplotypes=map(lambda x: x+1, range(len(allelecount)))
                r = vcf.parser._Record(CHROM='chrMT', POS=variant[0], ID='.', REF=variant[1], ALT=variant[3], QUAL='.', FILTER='PASS', INFO=OrderedDict([('AC',allelecount),('AN', len(variant[1])+1)]), FORMAT='GT:DP:HF:CILOW:CIUP:SDP', sample_indexes={sample:''},samples=[])
                #print variant[7], variant[6]
                r._sample_indexes[sample]=[[0]+aplotypes, variant[2],variant[7], variant[8], variant[9], variant[5]]
                #print r._sample_indexes
                r.samples.append(sample)
                if len(variant[3])>1:
                    r.REF=r.REF*len(variant[3])
                VCF_RECORDS.append(r)
                present_pos.add(r.POS)
                r.TYPEVAR=[variant[-1]]*len(variant[3])
            # if the v. position was never encountered before, is heteroplasmic and is not a deletion
            elif variant[0] not in present_pos and max(variant[7])<1 and variant[-1]!='del':
                allelecount=[1]*len(variant[3])
                aplotypes=map(lambda x: x+1, range(len(allelecount)))
                r = vcf.parser._Record(CHROM='chrMT', POS=variant[0], ID='.', REF=[variant[1]], ALT=variant[3], QUAL='.', FILTER='PASS', INFO=OrderedDict([('AC',allelecount),('AN', len(variant[3])+1)]), FORMAT='GT:DP:HF:CILOW:CIUP:SDP', sample_indexes={sample:''},samples=[])
                #print variant[6], variant[7]
                r._sample_indexes[sample]=[[0]+aplotypes, variant[2],variant[7], variant[8], variant[9], variant[5]]
                r.samples.append(sample)
                if len(variant[3])>1:
                    r.REF=r.REF*len(variant[3])
                VCF_RECORDS.append(r)
                present_pos.add(r.POS)
                r.TYPEVAR=[variant[-1]]*len(variant[3])
                #print r.POS, sample
            # if the v. position was never encountered before,is homoplasmic and is a deletion
            elif variant[0] not in present_pos and max(variant[7])>=1 and variant[-1]=='del':
                allelecount=[1]*len(variant[1])
                r = vcf.parser._Record(CHROM='chrMT', POS=variant[0], ID='.', REF=variant[1], ALT=variant[3], QUAL='.', FILTER='PASS', INFO=OrderedDict([('AC',allelecount),('AN',len(variant[1]))]), FORMAT='GT:DP:HF:CILOW:CIUP:SDP', sample_indexes={sample:''}, samples=[])
                genotype=range(1,len(variant[1])+1) #returns a list of genotypes where one is homopolasmic [1,2,...]
                r._sample_indexes[sample]=[genotype,variant[2], variant[7], variant[8], variant[9], variant[5]]
                r.samples.append(sample)
                if len(variant[3])>1:
                    r.REF=r.REF*len(variant[3])
                VCF_RECORDS.append(r)
                present_pos.add(r.POS)
                r.TYPEVAR=[variant[-1]]*len(variant[3])
            # if the v. position was never encountered before,is homoplasmic and is not a deletion
            elif variant[0] not in present_pos and max(variant[7])>=1 and variant[-1]!='del':
                allelecount=[1]*len(variant[3])
                r = vcf.parser._Record(CHROM='chrMT', POS=variant[0], ID='.', REF=[variant[1]], ALT=variant[3], QUAL='.', FILTER='PASS', INFO=OrderedDict([('AC',allelecount),('AN',len(variant[3]))]), FORMAT='GT:DP:HF:CILOW:CIUP:SDP', sample_indexes={sample:''}, samples=[])
                genotype = range(1,len(variant[3])+1) #returns a list of genotypes where one is homopolasmic [1,2,...]
                r._sample_indexes[sample]=[genotype,variant[2], variant[7], variant[8],variant[9], variant[5]]
                r.samples.append(sample)
                if len(variant[3])>1:
                    r.REF=r.REF*len(variant[3])
                VCF_RECORDS.append(r)
                present_pos.add(r.POS)
                r.TYPEVAR=[variant[-1]]*len(variant[3])
                #print r.POS, sample
            #If the v.position was encountered before
            elif variant[0] in present_pos and max(variant[7])<1:
		if variant[-1] == 'ins':
			print variant[0], variant[-1]
                for i in VCF_RECORDS:
                    if variant[0] == i.POS:
                        #print i
                        #when there are multiple variants for a position of the same individual
                        if sample in i.samples and type(variant[1]) == type(list()):
                            for x in xrange(len(variant[3])):
                                if variant[3][x] in i.ALT and variant[1][x] in i.REF:
                                    index=i.ALT.index(variant[3][x])
                                    i.INFO['AC'][index]+=1
                                    i.INFO['AN'] += 1
                                    aplotype=index+1
                                    i._sample_indexes[sample][0].append(aplotype)
                                    i._sample_indexes[sample][2].append(variant[7][x])
                                    i._sample_indexes[sample][3].append(variant[8][x])
                                    i._sample_indexes[sample][4].append(variant[9][x])
                                    i._sample_indexes[sample][5].append(variant[5][x])
                                elif variant[3][x] in i.ALT and variant[1][x] not in i.REF:
                                    i.INFO['AC'].append(1)
                                    i.ALT.append(variant[3][x])
                                    i.REF.append(variant[1][x])
                                    i.INFO['AN'] += 1
                                    index=len(i.ALT)-1 #the alt allele added to i.ALT is the last index
                                    aplotype=len(i.INFO['AC'])
                                    i.TYPEVAR.append(variant[-1])
                                    i._sample_indexes[sample][0].append(aplotype)
                                    i._sample_indexes[sample][2].append(variant[7][x])
                                    i._sample_indexes[sample][3].append(variant[8][x])
                                    i._sample_indexes[sample][4].append(variant[9][x])
                                    i._sample_indexes[sample][5].append(variant[5][x])
                                else:
                                    i.REF.append(variant[1][x])
                                    #print i.REF, variant[1], i.ALT
                                    i.ALT.append(variant[3][x])
                                    i.INFO['AC'].append(1)
                                    i.INFO['AN'] += 1
                                    index=i.ALT.index(variant[3][x])
                                    i.TYPEVAR.append(variant[-1])
                                    aplotype=index+1
                                    i._sample_indexes[sample][0].append(aplotype)
                                    i._sample_indexes[sample][2].append(variant[7][x])
                                    i._sample_indexes[sample][3].append(variant[8][x])
                                    i._sample_indexes[sample][4].append(variant[9][x])
                                    i._sample_indexes[sample][5].append(variant[5][x])
                                    #print i
                        #for multiple variants of a position in different individuals
                        elif sample not in i.samples and type(variant[1]) == type(list()):
                            i.INFO['AN'] += 1
                            i.samples.append(sample)
                            for x in xrange(len(variant[3])):
                                if variant[3][x] in i.ALT and variant[1][x] in i.REF:
                                    index=i.REF.index(variant[1][x])
                                    i.INFO['AC'][index] += 1
                                    i.INFO['AN']+=1
                                    aplotype=index+1
                                    genotype=[aplotype]
                                elif variant[3][x] in i.ALT and variant[1][x] not in i.REF:
                                    i.INFO['AC'].append(1)
                                    i.ALT.append(variant[3][x])
                                    i.REF.append(variant[1][x])
                                    i.INFO['AN'] += 1
                                    index=i.REF.index(variant[1][x])
                                    aplotype=len(i.INFO['AC'])
                                    genotype=[aplotype]
                                    i.TYPEVAR.append(variant[-1])
                                else:
                                    i.INFO['AC'].append(1)
                                    i.ALT.append(variant[3][x])
                                    i.REF.append(variant[1][x])
                                    i.INFO['AN'] += 1
                                    index=i.ALT.index(variant[3][x])
                                    aplotype=index+1
                                    genotype=[aplotype]
                                    i.TYPEVAR.append(variant[-1])
                            i._sample_indexes.setdefault(sample,[[0]+genotype, variant[2], variant[7], variant[8], variant[9], variant[5]])
                        elif sample in i.samples and type(variant[1]) != type(list()):
                            for allele in variant[3]:
                                if allele not in i.ALT:
                                    i.REF.append(variant[1])
                                    i.ALT.append(allele)
                                    i.INFO['AC'].append(1)
                                    i.INFO['AN']+=1
                                    index=i.ALT.index(allele)
                                    aplotype=index+1
                                    hf_index=variant[3].index(allele)
                                    if type(i._sample_indexes[sample][0]) == type(list()):
                                    	i._sample_indexes[sample][0].append(aplotype)
                                    else:
                                        hap = i._sample_indexes[sample][0]
                                        i._sample_indexes[sample][0] = [hap].append(aplotype)
                                    #print i._sample_indexes[sample], variant, hf_index, i.POS
                                    i._sample_indexes[sample][2].append(variant[7][hf_index])
                                    i._sample_indexes[sample][3].append(variant[8][hf_index])
                                    i._sample_indexes[sample][4].append(variant[9][hf_index])
                                    if i.POS == 3107 and variant[10] == 'mism':
                                        i._sample_indexes[sample][5].append(variant[5][0])
                                    elif 'und' in variant[5][0] and len(variant[5])!= len(variant[3]): #this is the case when IUPAC ambiguity is in mt-table
                                       variant[5] = variant[5]*len(variant[3]) #not a very nice trick but this will be repleaced in the new v. caller 
                                       i._sample_indexes[sample][5].append(variant[5][hf_index])
                                    else:
                                        i._sample_indexes[sample][5].append(variant[5][hf_index])
                                    i.TYPEVAR.append(variant[-1])
                                else:
                                    index=i.ALT.index(allele)
                                    i.INFO['AC'][index]+=1
                                    i.INFO['AN']+=1
                                    aplotype=index+1
                                    i._sample_indexes[sample][0].append(aplotype)
                                    hf_index=variant[3].index(allele)
                                    i._sample_indexes[sample][2].append(variant[7][hf_index])
                                    i._sample_indexes[sample][3].append(variant[8][hf_index])
                                    i._sample_indexes[sample][4].append(variant[9][hf_index])
                                    i._sample_indexes[sample][5].append(variant[5][hf_index])
                        else:
                            i.INFO['AN'] += 1
                            i.samples.append(sample)
                            genotype=[]
                            #print i._sample_indexes
                            i._sample_indexes.setdefault(sample, [[0], variant[2], variant[7], variant[8],variant[9], variant[5]])
                            for allele in variant[3]:
                                if allele not in i.ALT:
                                    i.REF.append(variant[1])
                                    i.ALT.append(allele)
                                    i.INFO['AC'].append(1)
                                    i.INFO['AN']+=1
                                    index=i.ALT.index(allele)
                                    aplotype=index+1
                                    genotype.append(aplotype)
                                    i._sample_indexes[sample][0].append(aplotype)
                                    i.TYPEVAR.append(variant[-1])
                                else:
                                    index=i.ALT.index(allele)
                                    i.INFO['AC'][index] +=1
                                    i.INFO['AN'] += 1
                                    aplotype=index+1
                                    genotype.append(aplotype)
                                    i._sample_indexes[sample][0].append(aplotype)
            else:
                #for homoplasmic variants in a position encountered before
                for i in VCF_RECORDS:
                    if i.POS == variant[0]:
                        for x in range(len(variant[3])):
                            if variant[3][x] not in i.ALT:
                                i.INFO['AC'].append(1)
                                i.INFO['AN'] += 1
                                i.ALT.append(variant[3][x])
                                i.REF.append(variant[1][x])
                                i.samples.append(sample)
                                index=i.ALT.index(variant[3][x])
                                genotype=index+1
                                i.TYPEVAR.append(variant[-1])
                                if sample in i._sample_indexes:
                                    i._sample_indexes[sample][0].append(genotype)
                                    i._sample_indexes[sample][2].append(variant[7][x])
                                    i._sample_indexes[sample][3].append(variant[8][x])
                                    i._sample_indexes[sample][4].append(variant[9][x])
                                    i._sample_indexes[sample][5].append(variant[5][x])
                                else:
                                    i._sample_indexes.setdefault(sample,[genotype, variant[2], variant[7], variant[8], variant[9], variant[5]])
                                #if a deletion, add a further reference base
                                #print i
                                if type(variant[1])== type(list()):
                                    i.REF.extend(variant[1])
                                else:
                                    i.REF.append(variant[1])
                                #print i
                            else:
                                index = i.ALT.index(variant[3][x])
                                i.INFO['AC'][index] += 1
                                i.INFO['AN']+=1
                                i.samples.append(sample)
                                genotype=index+1
                                if genotype in i._sample_indexes[sample][0]: #if the alt change is already there then add 1 to the genotype
                                    genotype = genotype+1
                                    i.ALT.append(variant[3][x])
                                    i.REF.append(variant[1][x])
                                i.TYPEVAR.append(variant[-1])
                                if sample in i._sample_indexes:
                                    i._sample_indexes[sample][0].append(genotype)
                                    i._sample_indexes[sample][2].append(variant[7][x])
                                    i._sample_indexes[sample][3].append(variant[8][x])
                                    i._sample_indexes[sample][4].append(variant[9][x])
                                    i._sample_indexes[sample][5].append(variant[5][x])
                                else:                                
                                    i._sample_indexes.setdefault(sample,[genotype, variant[2], variant[7], variant[8], variant[9], variant[5]])
    for r in VCF_RECORDS:
        if len(r.REF)>1:
            setref=set(r.REF)
            if len(setref)>1:
                for x in xrange(len(r.TYPEVAR)):
                    ord=sorted(r.REF, key=lambda x:len(x))
                    if r.TYPEVAR[x] == 'ins':
                        r.ALT[x] = ord[-1]+r.ALT[x][1:]
                    elif r.TYPEVAR[x] == 'del':
                        ndel = len(r.REF[x][1:])
                        altdel=ord[-1][0:(len(ord[-1])-ndel)]
                        r.ALT[x] = altdel
                    else:
                        r.ALT[x] = r.ALT[x]+ord[-1][1:]
                count_repeated_alt_allele = sum(sp.unique(r.ALT,return_counts=True)[1]>1)
                if count_repeated_alt_allele > 0: #add exception when multiple reference alleles are needed to explain two or more variants of the same type (e.g. GCACA,GCA --> G,G)
                    bv = sp.unique(r.ALT,return_counts=True)[1]>1
                    n_repeated_alleles = sp.unique(r.ALT,return_counts=True)[1][bv].tolist()
                    n_repeated_alleles = n_repeated_alleles[0]
                    ref = sp.unique(ord).tolist()
                    r.REF =[ref[-1]]
                    c = -1
                    for x in range(n_repeated_alleles):
                        c = c - x #add references with backward iteration
                        r.REF.append(ref[c])
                    r.REF = ','.join(r.REF)
                    r.REF = [r.REF]
                else:
                    r.REF=[ord[-1]]
            else:
                r.REF=[r.REF[0]]
    #gets the index of each sample and assign the definitive genotype to each individual        
    for index, sample in enumerate(dict_of_dicts.keys()):
        for r in VCF_RECORDS:
            if sample in r.samples:
                r._sample_indexes[sample].append(index)
                if type(r._sample_indexes[sample][0])==type(list()):
                    genotype=map(lambda x:str(x), r._sample_indexes[sample][0])
                    aplotypes="/".join(genotype)
                    r._sample_indexes[sample][0]=aplotypes
                    

    #counts also alleles identical to the reference base
    INDEX=OrderedDict()
    for index, sample in enumerate(dict_of_dicts.keys()):
        INDEX.setdefault(sample, index)

    for samples in INDEX:
        for r in VCF_RECORDS:
            if samples not in r._sample_indexes.keys():
                r._sample_indexes.setdefault(samples,[0, INDEX[samples]])
                r.INFO['AN']+=1
    #VCF header
    header=OrderedDict()
    c=0
    for x in dict_of_dicts:
        header.setdefault(x,c)
        c+=1

    header="\t".join(header.keys())

    #writes variant call in the VCF file
    #out=open("VCF_file.vcf","w")
    out = open(name+'.vcf','w')
    out.write('##fileformat=VCFv4.0\n##reference=chr%s\n' % reference)
    out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    out.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Reads covering the REF position">\n')
    out.write('##FORMAT=<ID=HF,Number=.,Type=Float,Description="Heteroplasmy Frequency of variant allele">\n')
    out.write('##FORMAT=<ID=CILOW,Number=.,Type=Float,Description="Value defining the lower limit of the confidence interval of the heteroplasmy fraction">\n')
    out.write('##FORMAT=<ID=CIUP,Number=.,Type=Float,Description="Value defining the upper limit of the confidence interval of the heteroplasmy fraction">\n')
    out.write('##FORMAT=<ID=SDP,Number=.,Type=String,Description="Strand-specific read depth of the ALT allele">\n')
    out.write('##INFO=<ID=AC,Number=1,Type=Integer,Description="Allele count in genotypes">\n')
    out.write('##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">\n')      
    out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+header+'\n')

    for position in sorted(present_pos):
        for r in VCF_RECORDS:
            if position == r.POS:
                if len(r.INFO['AC'])>1:
                    alleles=map(lambda x: str(x), r.INFO['AC'])
                    alleles=",".join(alleles)
                    AC='AC='+alleles
                    AN='AN='+str(r.INFO['AN'])
                else:
                    AC='AC='+str(r.INFO['AC'][0])
                    AN='AN='+str(r.INFO['AN'])
                samples_per_position=[]
                r._sample_indexes=sorted(r._sample_indexes.items(),key=lambda x: x[1][-1])
                #print r._sample_indexes
                for items in r._sample_indexes:
                    #print items
                    if len(items[1])>2:
                        if len(items[1][2])>1:
                            het=map(lambda x:str(x), items[1][2])
                            heteroplasmy=",".join(het)
                            CILOW=map(lambda x:str(x), items[1][3])
                            CIUP=map(lambda x:str(x), items[1][4])
                            #print CILOW,CIUP
                            confidence_interval_low=",".join(CILOW)
                            confidence_interval_up=",".join(CIUP)
                            strandp=[]
                            for x in items[1][5]:
                                if isinstance(x,list):
                                    strandp.append(x[0])
                                else:
                                    strandp.append(str(x))
                            SDP=",".join(strandp)
                            individual=str(items[1][0])+':'+str(items[1][1])+':'+heteroplasmy+':'+confidence_interval_low+':'+confidence_interval_up+':'+SDP
                        else:
                            heteroplasmy=str(items[1][2][0])                        
                            confidence_interval_low=str(items[1][3][0])
                            confidence_interval_up=str(items[1][4][0])
                            if isinstance(items[1][5][0],list):
                                SDP=items[1][5][0][0]
                            else:
                                SDP=str(items[1][5][0])
                            individual=str(items[1][0])+':'+str(items[1][1])+':'+heteroplasmy+':'+confidence_interval_low+':'+confidence_interval_up+':'+SDP
                        samples_per_position.append(individual)
                    else:
                        individual=str(items[1][0])
                        samples_per_position.append(individual)
                samples="\t".join(samples_per_position)
                if len(r.ALT)>1:
                    var=",".join(r.ALT)
                    out.write(r.CHROM+'\t'+str(r.POS)+'\t'+r.ID+'\t'+r.REF[0]+'\t'+var+'\t'+r.QUAL+'\t'+r.FILTER+'\t'+AC+';'+AN+'\t'+r.FORMAT+'\t'+samples+'\n')
                else:
                    out.write(r.CHROM+'\t'+str(r.POS)+'\t'+r.ID+'\t'+r.REF[0]+'\t'+r.ALT[0]+'\t'+r.QUAL+'\t'+r.FILTER+'\t'+AC+';'+AN+'\t'+r.FORMAT+'\t'+samples+'\n')
                    
    out.close()


def FASTAoutput(Consensus, mtDNAseq, names):
    path = os.getcwd()
    fasta_dict2={}
    for name2 in names:
        fasta_dict2[name2]=[]
    for name2 in fasta_dict2:
        for i in xrange(len(mtDNAseq)):
            index=i
            val=(index, mtDNAseq[i])
            fasta_dict2[name2].append(val)
    for name2 in Consensus:
        for variants in Consensus[name2]:
            if variants[-1]=='ins':
                var_pos=(int(variants[0])-1)+(float(len(variants[1]))/10)
                tupla=(var_pos, variants[1])
                fasta_dict2[name2].append(tupla)
            elif variants[-1]=='del':
                for x in variants[1]:
                    for j in fasta_dict2[name2]:
                        if x == j[0]:
                            fasta_dict2[name2].remove(j)
            else:
                for j in fasta_dict2[name2]:
                    if variants[0] == j[0]:
                        index = fasta_dict2[name2].index(j)
                        fasta_dict2[name2][index] = (variants[0], variants[1])
    for name2 in fasta_dict2:
        for dirname, dirnames, filenames in os.walk('.'):
            for subdirname in dirnames:
                if subdirname.startswith('OUT') and subdirname == names[name2]:
                    fasta_dir=glob.glob(os.path.join(path+'/'+subdirname))[0]				
                    fasta_out=open(fasta_dir+'/'+name2+'.fasta', "w")
                    fasta_out.write('>'+name2+'_complete_mitochondrial_sequence\n')
                    seq=[]
                    fasta_dict2[name2]=sorted(fasta_dict2[name2])
                    for tuples in fasta_dict2[name2]:
                        seq.append(tuples[1][0])
                    seq_def=''.join(seq)
                    fasta_out.write(seq_def)
                    fasta_out.close()

if __name__ == '__main__':
	print "This script is used only when called by assembleMTgenome.py."
	pass

