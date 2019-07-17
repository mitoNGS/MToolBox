#!/usr/bin/python

# -*- coding: UTF-8 -*-

#	Written by Ernesto Picardi - e.picardi@biologia.uniba.it
#	Edited by Claudia Calabrese - claudia.calabrese23@gmail.com
#	Edited by Domenico Simone - dome.simone@gmail.com


import getopt, sys, os

def usage():
	print """Map Nanopore FASTQ onto mtDNA
		Options:
		-a		Input Fastq (Nanopore only)
		-r		mtDNA reference name
		-q		Mapping Quality [30]
		-g		GMAP executable [/usr/local/bin/gmap]
		-D		GMAP database location [/usr/local/share]
		-M		GMAP database for mtDNA [chrRSRS]
		-H		GMAP database for complete human genome [hg19RSRS]
		-t		GMAP threads [8]
		-o		Out folder
		"""

try:
	opts, args = getopt.getopt(sys.argv[1:], "ha:r:q:g:D:M:H:t:o:")
except getopt.GetoptError, err:
	print str(err)
	usage()
	sys.exit()

if len(opts)==0:
	usage()
	sys.exit()
fastq1=None

gmapexe='/usr/local/bin/gmap'
gmapdb='/usr/local/share/gmapdb'
mtdb='chrRSRS'
humandb='hg19RSRS'
mapq=30
thread=8
mtref_name='chrM'
folder=os.path.join(os.getcwd(),'OUTfolder2')

for o,a in opts:
	if o == "-h":
		usage()
		sys.exit()
	elif o == "-a": fastq1 = a
	elif o == "-r": mtref_name = a
	elif o == "-g": gmapexe = a
	elif o == "-D": gmapdb = a
	elif o == "-M": mtdb = a
	elif o == "-H": humandb = a
	elif o == "-q": mapq = int(a)
	elif o == "-t": thread = int(a)
	elif o == "-o": folder = a
	else:
		assert False, "unhandled option"
def rev(seq):
	d={'A':'T','T':'A','C':'G','G':'C','N':'N'}
	s=''.join([d[x] for x in seq])
	return s[::-1]


if not os.path.exists(folder): os.mkdir(folder)

RG_tag = '--read-group-id=sample --read-group-name=sample --read-group-library=sample --read-group-platform=sample'

map1cmd='%s -D %s -d %s --nosplicing -f samse %s --nofails -n 1 -O -t %i %s > %s 2> %s' %(gmapexe,gmapdb,mtdb,RG_tag,thread,fastq1,os.path.join(folder,'outmt.sam'),os.path.join(folder,'logmt.txt'))

print 'Mapping onto mtDNA...'
print map1cmd
os.system(map1cmd)

print 'Extracting FASTQ from SAM...'
mtoutsam=os.path.join(folder,'outmt.sam')
dics={}
f=open(mtoutsam)
for i in f:
	# original version
	# if i.strip()=='': continue
	if i.strip()=='' or i.startswith('@'): continue
	l=(i.strip()).split('\t')
	if l[2]=='*': continue
	if dics.has_key(l[0]): dics[l[0]].append(l)
	else: dics[l[0]]=[l]
f.close()
single,pair1,pair2=[],[],[]

for i in dics:
	ll=dics[i]
	if len(ll)==1:
		strand,seq,qual=int(ll[0][1]) & 16,ll[0][9],ll[0][10]
		if strand==16: seq,qual=rev(seq),qual[::-1]
		entry='\n'.join(['@'+ll[0][0],seq,'+',qual])+'\n'
		single.append(entry)
	else:
		strand,seq,qual=int(ll[0][1]) & 16,ll[0][9],ll[0][10]
		if strand==16: seq,qual=rev(seq),qual[::-1]
		entry='\n'.join(['@'+ll[0][0],seq,'+',qual])+'\n'
		pair1.append(entry)
		strand,seq,qual=int(ll[1][1]) & 16,ll[1][9],ll[1][10]
		if strand==16: seq,qual=rev(seq),qual[::-1]
		entry='\n'.join(['@'+ll[1][0],seq,'+',qual])+'\n'
		pair2.append(entry)

sig=0
if len(single)!=0:
	mtoutfastq=os.path.join(folder,'outmt.fastq')
	out=open(mtoutfastq,'w')
	out.writelines(single)
	out.close()
	sig=1
else:
	print 'no mitochondrial reads found\n'
	sys.exit(1)

if len(pair1)!=0:
	print 'something wrong in here\n...paired end reads not expected for Nanopore.\nAre you sure you are using the right MToolBox version?\n'

if sig:
	print 'Mapping on complete human genome...single reads'
	map2cmd='%s -D %s -d %s --nosplicing -f samse --nofails --no-sam-headers -x 1 -O -t %i %s --split-output %s 2> %s' %(gmapexe,gmapdb,humandb,thread,mtoutfastq,os.path.join(folder,'outhumanS'),os.path.join(folder,'loghumanS.txt'))
	os.system(map2cmd)
	print map2cmd

	print 'Filtering reads...'

	hgoutsam=os.path.join(folder,'outhumanS.uniq')
	dicsingle={}
	f=open(hgoutsam)
	for i in f:
		if i.strip()=='': continue
		l=(i.strip()).split('\t')
		#if l[2]=='*': continue
		if l[2] == mtref_name and int(l[4]) >= mapq: #filter reads mapping on mitochondrial DNA only and above the mapping quality cutoff
			if dicsingle.has_key(l[0]):
				print 'skip read %s because duplicated in outhumanS.uniq' %(l[0])
				#dicsingle[l[0]].append(l)
			else:
				dicsingle[l[0]]=[l]
	f.close()

	print 'Filtering reads...done'

	finalsam=os.path.join(folder,'OUT2.sam')
	out=open(finalsam,'w')

	for k,v in dicsingle.iteritems():
		print v[0]
		o = '\t'.join(v[0])+'\n'
		out.write(o)
	out.close()

print 'Outfile saved on %s.' %(finalsam)
print 'Done.'

