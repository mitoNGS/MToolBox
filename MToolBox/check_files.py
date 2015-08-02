#!/usr/bin/env python

# USAGE
# check_files.py -t $mt_classifier_OPTS -a $assembleMTgenome_OPTS -m $mapExome_OPTS

import getopt, argparse, os, sys

####
# printed messages
sys_message = "The default file/directory %s will be used."
yes_message = "Specified file/directory %s is ok."
no_message = "Specified file/directory %s does not exist and the default file %s was not found. Abort."
####

def check_file(user_file, sys_file):
	f = True
	#	mtdb_fasta = user_file
	#	try:
	#		if os.path.exists(user_file):
	#	except
	
	if os.path.exists(user_file):
		print yes_message % user_file
	#print >>sys.stdout, yes_message % user_file
	else:
		if os.path.exists(sys_file):
			print sys_message % sys_file
		#print >>sys.stdout, sys_message % sys_file
		else:
			f = False
			print no_message % (user_file, sys_file)
	#print >>sys.stdout, no_message % (user_file, sys_file)
	return f

# Parse check_files command line
parser = argparse.ArgumentParser(description="Check files.")
parser.add_argument("--mapExome_OPTS", help="mapExome options", \
					default="-g /usr/local/bin/gsnap -D /usr/local/share/gmapdb/ -M chrRSRS -H hg19RSRS")

parser.add_argument("--assembleMTgenome_OPTS", help="assembleMTgenome options", \
					default="-r /usr/local/share/genomes/ -f chrRSRS.fa -a hg19RSRS.fa -s /usr/local/bin/samtools")

parser.add_argument("--mt_classifier_OPTS", help="mt-classifier options", \
					default="-m /usr/local/bin/muscle")

f = parser.parse_args()
#print f
#print "\nmapExome argument user-defined"
user_mapExome_opts, user_mapExome_args = getopt.getopt(f.mapExome_OPTS.split(), "ha:b:c:g:D:M:H:t:o:")
d_user_mapExome = dict(user_mapExome_opts)
#print d_user_mapExome
#print "\nmapExome argument defaults"
default_mapExome_OPTS = parser.get_default("mapExome_OPTS")
#print default_mapExome_OPTS
default_mapExome_opts, default_mapExome_args = getopt.getopt(default_mapExome_OPTS.split(), "ha:b:c:g:D:M:H:t:o:")
d_def_mapExome = dict(default_mapExome_opts)
#print d_def_mapExome

#print "\nassembleMTgenome argument user-defined"
user_assembleMTgenome_opts, user_assembleMTgenome_args = getopt.getopt(f.assembleMTgenome_OPTS.split(), "hf:i:q:c:d:o:g:a:r:s:FCUPNA:D:z:t:")
d_user_assembleMTgenome = dict(user_assembleMTgenome_opts)
#print d_user_assembleMTgenome
#print "\nassembleMTgenome argument defaults"
default_assembleMTgenome_OPTS = parser.get_default("assembleMTgenome_OPTS")
default_assembleMTgenome_opts, default_assembleMTgenome_args = getopt.getopt(default_assembleMTgenome_OPTS.split(), "hf:i:q:c:d:o:g:a:r:s:FCUPNA:D:z:t:")
d_def_assembleMTgenome = dict(default_assembleMTgenome_opts)
#print d_def_assembleMTgenome

#print "\nmt-classifier argument user-defined"
user_mt_classifier_opts, user_mt_classifier_args = getopt.getopt(f.mt_classifier_OPTS.split(), "hm:b:i:s:")
d_user_mt_classifier = dict(user_mt_classifier_opts)
#print d_user_mt_classifier
#print "\nmt_classifier argument defaults"
default_mt_classifier_OPTS = parser.get_default("mt_classifier_OPTS")
default_mt_classifier_opts, default_mt_classifier_args = getopt.getopt(default_mt_classifier_OPTS.split(), "hm:b:i:s:")
d_def_mt_classifier = dict(default_mt_classifier_opts)
#print d_def_mt_classifier


print "\nChecking mapExome parameters..."
for i in d_user_mapExome:
	if i in (d_def_mapExome.keys()):
		if i in ("-M", "-H"): # particular cases: gsnap databases
			if check_file(d_user_mapExome['-D']+d_user_mapExome[i], d_user_mapExome['-D']+d_user_mapExome[i]) == True: pass
			else: sys.exit(1)
		else: # all the other cases
			if check_file(d_user_mapExome[i], d_def_mapExome[i]) == True: pass
			else: sys.exit(1)
print "OK."

print "\nChecking assembleMTgenome parameters..."
for i in d_user_assembleMTgenome:
	if i in (d_def_assembleMTgenome.keys()):
		# newly added to accommodate the possibility to change reference sequence
		if i in ("-f", "-a"):
			if check_file(d_user_assembleMTgenome['-r']+d_user_assembleMTgenome[i], d_user_assembleMTgenome['-r']+d_user_assembleMTgenome[i]) == True:
				pass
				print d_user_assembleMTgenome['-r']+d_user_assembleMTgenome[i]
			else: sys.exit(1)			
		if check_file(d_user_assembleMTgenome[i], d_def_assembleMTgenome[i]) == True: pass
		else: sys.exit(1)
print "OK."

print "\nChecking mt-classifier parameters..."
for i in d_user_mt_classifier:
	if i in (d_def_mt_classifier.keys()):
		if check_file(d_user_mt_classifier[i], d_def_mt_classifier[i]) == True: pass
		else: sys.exit(1)
print "OK."

