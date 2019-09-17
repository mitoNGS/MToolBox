#!/usr/bin/env python

import fileinput
import sys, os

def usage():
	print ''' 

This script is compatible only with MToolBox v.1.2.

The script filters the MToolBox vcf file based on Heteroplasmy and Read Depth thresholds

Usage:
filter_HF.py <sample_name>  <vcf_file> <HF_threshold[float]> <DP_threshold[float]> <out_type[vcf|txt]> <outfilename> <convert_to_homoplamy[Yes|No]> \n<vcf_file> can also be .gz file\n\n<convert_to_homoplasmy> is boolean and takes Yes or No values and converts HF >= 0.9 to GT=1/1. Useful for haplogroup prediction with other methods (e.g. haplogrep)\n\n'''


def check_inconsistency(genotype, HF, CIUP, CILOW, SDP, sample_name):
	if len(genotype) != len(HF):
		diff = len(genotype) - len(HF)
		HF.extend(['und']*diff)
		sys.stderr.write('Unable to determine HF for pos {0} in sample {1}\n'.format(line[1], sample_name))
	if len(genotype) != len(CIUP):
		diff = len(genotype) - len(CIUP)
		CIUP.extend(['und']*diff)
		sys.stderr.write('Unable to determine CILOW for pos {0} in sample {1}\n'.format(line[1], sample_name))
	if len(genotype) != len(CILOW):
		diff = len(genotype) - len(CILOW)
		CILOW.extend(['und']*diff)
		sys.stderr.write('Unable to determine CIUP for pos {0} in sample {1}\n'.format(line[1], sample_name))
	if len(genotype) != len(SDP):
		diff = len(genotype) - len(SDP)
		SDP.extend(['und;und']*diff)
		sys.stderr.write('Unable to determine SDP for pos {0} in sample {1}\n'.format(line[1], sample_name))
	return HF, CIUP, CILOW, SDP
	

if __name__ == "__main__":
	if len(sys.argv[1:]) < 7:
		sys.stderr.write('ERROR: argument missing\n')
		usage()
		sys.exit(1)

	samplename,vcf,HFt,DPt,out_type,outfile,homo_convert= sys.argv[1:]
	HFt = float(HFt)
	DPt = float(DPt)
	out = open(outfile,'w')
	homo_convert = str(homo_convert)
	if homo_convert not in ['Yes','No']:		
		sys.stderr.write('Values accepted for <convert_to_homoplasmy> are [Yes|No].\nExit!\n')
		sys.exit(1)
	if 'gz' in vcf or 'gzip' or 'bz2' in vcf:
		ifile = fileinput.input(vcf,openhook=fileinput.hook_compressed)
	else:	
		ifile = fileinput.input(vcf)
	for line in ifile:
		if line.startswith('##'):
			if out_type == 'vcf':
				command_string = "##contig=<ID=chrMT,length=16569>\n##filter_VCF_command=filter_vcf.py {0} {1} {2} {3} {4} {5}\n".format(vcf,HFt,DPt,out_type,outfile,homo_convert)
				out.write(line)
			else:
				pass
		else:
			if line.startswith('#CHROM') and out_type == 'vcf':
				out.write(command_string)
				line = line.split('\t')
				line[-1] = samplename+'\n'
				line = '\t'.join(line)
				out.write(line)
			elif line.startswith('#CHROM') and out_type == 'txt':
				header='CHROM\tPOS\tID\tREF\tALT\tDP\tHF\tCIL\tCIU\tSDP_R\tSDP_R\t'+samplename
				out.write(header+'\n')
			else:
				line = line.split('\t')
				geno,DPv,HFv_l,CIL,CIU,SDP = line[-1].split(':')
				geno = geno.split('/')
				if '0' in geno:
					geno.remove('0')
				HFv_l = HFv_l.split(',')
				CIL = CIL.split(',')
				CIU = CIU.split(',')
				SDP = SDP.split(',')
				ALT = line[4].split(',')
				HFv_l,CIL,CIU,SDP = check_inconsistency(geno, HFv_l, CIL, CIU, SDP,samplename)
				c =0
				while c < (len(geno)):
					HFv = float(HFv_l[c])
					CILv = float(CIL[c])
					CIUv = float(CIU[c])
					DPv = float(DPv)
					ALTv = str(ALT[c])
					try:
						SDP_F = int(SDP[c].split(';')[0])
						SDP_R = int(SDP[c].split(';')[1])
					except:
						#case when SDP is undetermined
						SDP_F = SDP[c].split(';')[0].strip()
						SDP_R = SDP[c].split(';')[1].strip()
					if DPv >= float(DPt) and HFv >= float(HFt):
						if out_type == 'txt':
							res='\t'.join(map(lambda x:str(x),[line[0],line[1],line[2],line[3],ALTv,DPv,HFv,CILv,CIUv,SDP_F,SDP_R,samplename]))
							out.write(res+'\n')
						else:
							if HFv == 1:
								res='\t'.join(map(lambda x:str(x),[line[0],line[1],line[2],line[3],ALTv,'.','PASS','AC=2,AN=2','GT','1/1']))
							elif HFv >= 0.9 and homo_convert == 'Yes':
								res='\t'.join(map(lambda x:str(x),[line[0],line[1],line[2],line[3],ALTv,'.','PASS','AC=2,AN=2','GT','1/1']))
							else:
								res='\t'.join(map(lambda x:str(x),[line[0],line[1],line[2],line[3],ALTv,'.','PASS','AC=1,AN=2','GT','0/1']))
							out.write(res+'\n')
					else:
						pass
					c += 1

	out.close()
