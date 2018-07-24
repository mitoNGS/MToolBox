#!/usr/bin/python

"""
	Written by Ernesto Picardi - e.picardi@biologia.uniba.it
	"""

import getopt, sys, os

def filter_alignments(outmt, outhumanS, outhumanP, OUT):
    sig=1
    pai=1
    #single,pair1,pair2=[],[],[]

    # for i in dics:
    # 	ll=dics[i]
    # 	if len(ll)==1:
    # 		strand,seq,qual=int(ll[0][1]) & 16,ll[0][9],ll[0][10]
    # 		if strand==16: seq,qual=rev(seq),qual[::-1]
    # 		entry='\n'.join(['@'+ll[0][0],seq,'+',qual])+'\n'
    # 		single.append(entry)
    # 	else:
    # 		strand,seq,qual=int(ll[0][1]) & 16,ll[0][9],ll[0][10]
    # 		if strand==16: seq,qual=rev(seq),qual[::-1]
    # 		entry='\n'.join(['@'+ll[0][0],seq,'+',qual])+'\n'
    # 		pair1.append(entry)
    # 		strand,seq,qual=int(ll[1][1]) & 16,ll[1][9],ll[1][10]
    # 		if strand==16: seq,qual=rev(seq),qual[::-1]
    # 		entry='\n'.join(['@'+ll[1][0],seq,'+',qual])+'\n'
    # 		pair2.append(entry)
    print('Reading Results...')
    if sig:
        hgoutsam = outhumanS
        #hgoutsam=os.path.join(folder, outhumanS)
        dicsingle={}
        f=open(hgoutsam)
        for i in f:
            if i.strip()=='': continue
            l=(i.strip()).split('\t')
            if l[2]=='*': continue # the read is not mapped
            # keeping multiple mappings
            if dicsingle.has_key(l[0]):
                dicsingle[l[0]].append(l)
            else:
                dicsingle[l[0]]=[l]
        f.close()
    if pai:
        hgoutsam2 = outhumanP
        #hgoutsam2=os.path.join(folder,outhumanP.sam)
        dicpair={}
        f=open(hgoutsam2)
        for i in f:
            if i.strip()=='': continue
            l=(i.strip()).split('\t')
            if l[2]=='*': continue
            if dicpair.has_key(l[0]):
                dicpair[l[0]].append(l)
            else:
                dicpair[l[0]]=[l]
        f.close()

    print('Extracting FASTQ from SAM...')
    mtoutsam = outmt
    #mtoutsam=os.path.join(folder,outmt)
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
    
    finalsam = OUT
    #finalsam = os.path.join(folder,'OUT.sam')
    out=open(finalsam,'w')
    out.write("@SQ	SN:%s	LN:16569\n" % mtdb)
    out.write("@RG	ID:sample	PL:sample	PU:sample	LB:sample	SM:sample\n")

    
    print('Filtering reads...')
    #good=[]
    for i in dics:
        ll=dics[i]
        if len(ll)==1: # if the read has one mapping I assume it's SE
            if dicsingle.has_key(i):
                r=dicsingle[i] # i is a list of lists (splitted sam lines)
                #print ll
                #print r
                if len(r)==1: 
                    # check if read aligned on MT when aligned against nuclear+MT
                    # fields checked: RNAME, POS
                    if r[0][2]==ll[0][2] and ll[0][3]==r[0][3]:
                        #good.append('\t'.join(ll[0])+'\n')
                        finalsam.write('\t'.join(ll[0])+'\n')
            else:
                #good.append('\t'.join(ll[0])+'\n')
                finalsam.write('\t'.join(ll[0])+'\n')
        else:
            if dicpair.has_key(i):
                r=dicpair[i]
                if len(r) == 2:
                    if r[0][2]==ll[0][2] and ll[0][3]==r[0][3] and r[1][2]==ll[1][2] and ll[1][3]==r[1][3]:
                        finalsam.write('\t'.join(ll[0])+'\n')
                        finalsam.write('\t'.join(ll[1])+'\n')
                        #good.append('\t'.join(ll[0])+'\n')
                        #good.append('\t'.join(ll[1])+'\n')
            else:
                finalsam.write('\t'.join(ll[0])+'\n')
                finalsam.write('\t'.join(ll[1])+'\n')
                # good.append('\t'.join(ll[0])+'\n')
                # good.append('\t'.join(ll[1])+'\n')

    #finalsam=os.path.join(folder,'OUT.sam')
    #out=open(finalsam,'w')
    # add SAM Header and Read Group for GATK Indels realignment
    # out.write("@SQ	SN:%s	LN:16569\n" % mtdb)
    # out.write("@RG	ID:sample	PL:sample	PU:sample	LB:sample	SM:sample\n")
    #out.writelines(good)
    out.close()
    
    print('Outfile saved on %s.' %(finalsam))
    print('Done.')
