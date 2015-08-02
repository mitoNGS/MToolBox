# eventually to be included in consts.py (just for the lulz)

"""
script workflow:
- get best haplogroup
- find the corresponding macrohaplo in the MHCS list
    
    as to Apr 3rd 2013, the list is
    ['A', 'B', 'D', 'E', 'F', 'G', 'H', 'H1', 'H2', 'HV', 'I', 'J', 'J1', 'J2', 'K', 'L', 'L0', 'L1', 'L2', 'L3', 'M', 'N', 'R', 'R*', 'R0', 'T', 'T2', 'U', 'V', 'W', 'W7476', 'X']
    
    haplogroups for which 2nd character control must be performed:
    - H [H, H1, H2]
    - J [J, J1, J2]
    - L [L, L0, L1, L2, L3]
    - R [R, R*, R0]
    - T [T, T2]
    
	haplogroups not having their own MHCSs that should be matched with a upper level MHCS:
	recode_haplo = {'C':'M', 'Q':'M', 'Z':'M', 'O':'N', 'S':'N', 'Y':'N', 'P':'R', 'HV':'R0', 'V':'R0'}    
    
    =======================================================================================================================
    
	as to Feb 8th 2013, the list is:
	['A', 'AF', 'AM', 'AS', 'B', 'D', 'EU', 'F', 'G', 'H', 'H1', 'H2', 'I', 'J1', 'J2', 'K', 'L', 'L3', 'M', 'N', 'OC', 'R', 'R0', 'RSRS', 'T', 'T2', 'U', 'W', 'WHOLE', 'X', 'rCRS']
	upon discarding of "continental" MHCSs
	['A', 'B', 'D', 'F', 'G', 'H', 'H1', 'H2', 'I', 'J1', 'J2', 'K', 'L', 'L3', 'M', 'N', 'R', 'R0', 'T', 'T2', 'U', 'W', 'X']
	haplogroups for which 2nd character control must be performed:
	- H [H, H1, H2]
	- J [J1, J2]
	- L [L, L3]
	- R [R, R0]
	- T [T, T2]
	
	haplogroups not having their own MHCSs that should be matched with a upper level MHCS:
	recode_haplo = {'C':'M', 'E':'M', 'Q':'M', 'Z':'M', 'O':'N', 'S':'N', 'Y':'N', 'P':'R', 'HV':'R0', 'V':'R0'}
"""
import NGclassify
import datatypes
import os

# haplogroups not having their own MHCSs that should be matched with a upper level MHCS:
recode_haplo = {'C':'M', 'E':'M', 'Q':'M', 'Z':'M', 'O':'N', 'S':'N', 'Y':'N', 'P':'R', 'HV':'R0', 'V':'R0'}

def align_sequence(sequence, rif=None):
    if rif is None:
        rif = datatypes.Sequence('RSRS', consts.RCRS)
    seq_diff = NGclassify.SequenceDiff()
    print "Aligning sequence %s" % sequence.name
    seq_diff.gen_diff(rif, datatypes.Sequence(sequence.name, str(sequence.seq)))
    print "-"*30
    return seq_diff

def subparse2mhcs(best_haplo, mhcs_dict):
	# print "best_haplo is", best_haplo
	# print "best_haplo starts with", best_haplo[0]
	# print "best_haplo starts with", best_haplo[:2]
	if best_haplo[0] == 'H':
		if best_haplo[:2] == 'H1':
			name = 'H1'
			seq = mhcs_dict['H1']
		elif best_haplo[:2] == 'H2':
			name = 'H2'
			seq = mhcs_dict['H2']
		elif best_haplo[:2] == 'HV':
			name = 'HV'
			seq = mhcs_dict['HV']
		else:
			name = 'H'
			seq = mhcs_dict['H']
	elif best_haplo[0] == 'J':
		if best_haplo[:2] == 'J1':
			name = 'J1'
			seq = mhcs_dict['J1']
		elif best_haplo[:2] == 'J2':
			name = 'J2'
			seq = mhcs_dict['J2']
		else:
			name = 'J'
			seq = mhcs_dict['J']
	elif best_haplo[0] == 'L':
		if best_haplo[:2] == 'L3':
			name = 'L3'
			seq = mhcs_dict['L3']
		else:
			name = 'L'
			seq = mhcs_dict['L']
	elif best_haplo[0] == 'R':
		if best_haplo[:2] == 'R0':
			name = 'R0'
			seq = mhcs_dict['R0']
		else:
			name = 'R'
			seq = mhcs_dict['R']
	elif best_haplo[0] == 'T':
		if best_haplo[:2] == 'T2':
			name = 'T2'
			seq = mhcs_dict['T2']
		else:
			name = 'T'
			seq = mhcs_dict['T']
	else:
		name = None
		seq = None
	return name, seq

def parse2mhcs_dict(infile, mhcs_dict = {}):
	"""
    Parse the MHCS tab, whose format is
        
    <HAPLOGROUP>        <SEQ>
    
	Takes a list like
	[['seq1.id', 'seq1.seq'], ['seq2.id', 'seq2.seq']]
	and returns a dict like
	{'seq1.id' : 'seq1.seq', 'seq2.id' : 'seq2.seq'}
	"""
	inhandle = open(infile, 'r').readlines()
	inhandle = [i.strip().split() for i in inhandle]
	for row in inhandle:
		mhcs_dict[row[0].replace('_MHCS', '')] = row[1]
	return mhcs_dict

def which_mhcs_lite(best_haplo, mhcs_dict):
    """
    Takes as input the best haplogroup assignment and the mhcs sequence dictionary
    and returns the ID of the correspondent mhcs
    """
    f = best_haplo[0]
    if f in ['A', 'B', 'D', 'F', 'G', 'I', 'K', 'M', 'N', 'U', 'W', 'X']:
        # get mhcs code as it is
        name = f
        # seq = mhcs_dict[f]
    elif f in recode_haplo.keys():
        # haplogroup not having their own mhcs;
        # must be matched with a upper level mhcs
        name = recode_haplo[f]
        # seq = mhcs_dict[name]
    elif f in ['H', 'J', 'L', 'R', 'T']:
        print "I'm looking for", best_haplo
        # specific criteria for getting MHCS
        name, seq = subparse2mhcs(best_haplo,mhcs_dict)
    else:
        print "\nWHICH CAZZO\n"
    return name

def which_mhcs(best_haplo, mhcs_dict):
    """
    Takes as input the best haplogroup assignment and the mhcs sequence dictionary
    and returns a datatypes.Sequence object with data of the mhcs to be aligned
    """
    f = best_haplo[0]
    if f in ['A', 'B', 'D', 'F', 'G', 'I', 'K', 'M', 'N', 'U', 'W', 'X']:
        # get mhcs code as it is
        name = f
        seq = mhcs_dict[f]
    elif f in recode_haplo.keys():
        # haplogroup not having their own mhcs;
        # must be matched with a upper level mhcs
        name = recode_haplo[f]
        seq = mhcs_dict[name]
    elif f in ['H', 'J', 'L', 'R', 'T']:
        # specific criteria for getting mhcs
        name, seq = subparse2mhcs(best_haplo,mhcs_dict)
    return datatypes.Sequence(name, seq)

"""
# mhcs_dict is generated here --> add to consts module
mhcs_dict = parse2mhcs_dict('mhcs.tab')

# example of implementation
# query is the sequence to be analyzed
query = datatypes.Sequence('query', 'AACCCAAGTCAATAGAAGCCGGCGTAAAGACTGTTTTAGATGACCCCCTCCCCAATAAAGCTA')
p = which_mhcs('HV2a1', mhcs_dict)
print 'MHCS seq:', p.__dict__['name']
print p.__dict__['seq'][:50]

# SeqDiff object
print "SeqDiff object"
a = align_sequence(query, p)

# find SeqDiff object
print "find SeqDiff object"
# drop useless changes (gaps)
a.find_segment()
# print dir(a)
print "="*30
print a.pprint(outfile=open("ciao.txt", 'w'))
"""