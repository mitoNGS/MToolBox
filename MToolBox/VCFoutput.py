#!/usr/bin/env python

import getopt, os, ast, sys
from mtVariantCaller import VCFoutput

def usage():
	print """Produces VCF file output from VCF_dict_tmp file.
		Version 1.1 - Written by Domenico Simone and Claudia Calabrese - 2013-2014
		
		Options:
		-r		Reference sequence [RCRS]
		-s		VCF name
		
		If launched stand-alone needs the VCF_dict_tmp file to be positioned in the working directory
		"""

reference_sequence="RCRS"
sample_vcf_name='sample'
try:
	opts, args = getopt.getopt(sys.argv[1:], "h:r:s:")
except getopt.GetoptError, err:
	print str(err)
	usage()
	sys.exit()

for o,a in opts:
	if o == "-r":
		reference_sequence = a
	elif o == "-s":
		sample_vcf_name = str(a)


VCF_dict = ast.literal_eval(open('VCF_dict_tmp', 'r').read())
#VCFoutput(VCF_dict, reference=reference_sequence)
VCFoutput(VCF_dict, reference=reference_sequence, name=sample_vcf_name)
