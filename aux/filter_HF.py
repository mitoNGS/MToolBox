#!/usr/bin/env python

import fileinput
import sys, os

def usage():
	print ''' 
This script filter the MToolBox vcf file based on Heteroplasmy threshold
Usage:
filter_HF.py <sample_name>  <vcf_file> <HF_threshold[float]> <DP_threshold[float]> <out_type[vcf|txt]> <outfile.txt>\n<vcf_file> can also be .gz file\n\n '''


if __name__ == "__main__":
	if len(sys.argv[1:]) < 6:
		sys.stderr.write('ERROR: argument missing\n')
		usage()
		sys.exit(1)

	samplename,vcf,HFt,DPt,out_type,outfile = sys.argv[1:]
	HFt = float(HFt)
	DPt = float(DPt)
	out = open(outfile,'w')
	if 'gz' in vcf or 'gzip' or 'bz2' in vcf:
		ifile = fileinput.input(vcf,openhook=fileinput.hook_compressed)
	else:	
		ifile = fileinput.input(vcf)
	for line in ifile:
		if line.startswith('##'):
			if out_type == 'vcf':
				command_string = "##contig=<ID=chrMT,length=16569>\n##filter_VCF_command=filter_vcf.py {0} {1} {2} {3} {4}\n".format(vcf,HFt,DPt,out_type,outfile)
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
				header='CHROM\tPOS\tID\tREF\tALT\tDP\tHF\tCIL\tCIU\t'+samplename
				out.write(header+'\n')
			else:
				line = line.split('\t')
				geno,DPv,HFv_l,CIL,CIU = line[-1].split(':')
				geno = geno.split('/')
				if '0' in geno:
					geno.remove('0')
				HFv_l = HFv_l.split(',')
				CIL = CIL.split(',')
				CIU = CIU.split(',')
				ALT = line[4].split(',')
				c =0
				while c < (len(geno)):
					HFv = float(HFv_l[c])
					CILv = float(CIL[c])
					CIUv = float(CIU[c])
					DPv = float(DPv)
					ALTv = str(ALT[c])
					if DPv >= float(DPt) and HFv >= float(HFt):
						if out_type == 'txt':
							res='\t'.join(map(lambda x:str(x),[line[0],line[1],line[2],line[3],ALTv,DPv,HFv,CILv,CIUv,samplename]))
							out.write(res+'\n')
						else:
							if HFv == 1:
								res='\t'.join(map(lambda x:str(x),[line[0],line[1],line[2],line[3],ALTv,'.','PASS','AC=2,AN=2','GT','1/1']))
							else:
								res='\t'.join(map(lambda x:str(x),[line[0],line[1],line[2],line[3],ALTv,'.','PASS','AC=1,AN=2','GT','0/1']))
							out.write(res+'\n')
					else:
						pass
					c += 1

	out.close()
