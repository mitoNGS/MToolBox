MToolBox - README
=================

The MToolBox pipeline is a bioinformatics pipeline aimed at the analysis of mitochondrial DNA (mtDNA) in high throughput sequencing studies. It includes:

- an extended version of a previously published computational strategy for mtDNA genome assembly (PMID: 22669646). The pipeline has been integrated with the detection of insertions and deletions (indels) and the assessment of the heteroplasmic fraction (HF) and related confidence interval (CI) for each mt variant. HF and CI are integrated as genotype specific meta-information in a Variant Call Format (version 4.0) file;
- the mt-classifier tool, for haplogroup prediction, mt variant functional annotation and prioritization


SYSTEM REQUIREMENTS
===================

- UNIX-based OS
- Python2.7 (www.python.org)
- GSNAP (https://github.com/julian-gehring/GMAP-GSNAP) installed in /usr/local/bin/gsnap (otherwise you must specify the actual path using mapExome.py options)
- samtools (https://sourceforge.net/projects/samtools/files/samtools/) installed in /usr/local/bin/samtools (otherwise you must specify the actual path using assembleMTgenome.py options)
- MUSCLE (http://www.drive5.com/muscle/downloads.htm) installed in /usr/local/bin/muscle (otherwise you should specify the actual path using mt-classifier.py options)


MTOOLBOX PROGRAMS
=================

MToolBox.sh is the shell script invoking all the following programs:
- Sam2Fastq.jar to convert either SAM or BAM files to FASTQ format. This module is included in the ext_tools folder;
- mapExome.py to map reads on RSRS and hg19. It invokes GSNAP;
- SortSam.jar, MarkDuplicates.jar, SamFormatConverter.jar (PicardTools suite; http://picard.sourceforge.net) for SAM/BAM manipulation at several stages of the pipeline and elimination of PCR duplicates. These modules are included in the ext_tools folder;
- GenomeAnalysisTK.jar (GATK suite; PMID: 20644199) for indels realignment. This module is included in the ext_tools folder;
- assembleMTgenome.py to assemble the mitochondrial genome and perform variant calling and heteroplasmy quantification, by invoking:
- mpileup (samtools)
- mtVariantCaller.py module;
- VCFoutput.py to report variants in the VCF file (version 4.0). It invokes the vcf module (elease 0.60 PyVCF; https://githubcom/jamescasbon/PyVCF/);
- mt-classifier.py for haplogroup prediction;
- variants_functional_annotation.py for functional annotation and prioritization.


OTHER FILES REQUIRED 
====================

By default, MToolBox adopts the RSRS (Reconstructed Sapiens Reference Sequence, PMID: 22482806) as mitochondrial reference genome and hg19 as nuclear reference genome. Alternatively, the user can choose to use the rCRS (revised Cambridge Reference Sequence).
For user's convenience, fasta files and related gsnap indexes for the two reference sequences can be downloaded at https://sourceforge.net/projects/mtoolbox/.

List of fasta files and index .fai files in the genome_fasta folder (https://sourceforge.net/projects/mtoolbox/files/genome_fasta/) :
- chrRSRS.fa
- chrRSRS.fa.fai
- hg19RSRS.fa
- hg19RSRS.fa.fai

- chrRCRS.fa
- chrRCRS.fa.fai
- hg19RCRS.fa
- hg19RCRS.fa.fai

The default directory for these files is /usr/local/share/genomes/, otherwise you must specify the actual path using the assembleMTgenome.py options -r (for further details please refer to "RUNNING MTOOLBOX" section). 

List of compressed files of gsnap indexes in the genome_index folder (https://sourceforge.net/projects/mtoolbox/files/genome_index/):
- chrRSRS.tar.gz
- hg19RSRS.tar.gz 

The default directory for decompressed folders is /usr/local/share/gmapdb/, otherwise you must specify the actual path using the mapExome.py options -M and -H, respectively (for further details please refer to "RUNNING MTOOLBOX" section).

The MToolBox folder includes the MITOMAP_HMTDB_known_indels.vcf file, containing 127 known indels annotated in MITOMAP and HMTDB and the related intervals_file.list used by GATK GenomeAnalysisTK.jar module.


INSTALLATION
============

Download MToolBox.tar.gz. Decompress the file and copy the package folder in a folder of your choice. Add this folder to your system PATH with the following command:

export PATH=/dir/MToolBox/:$PATH


RUNNING MTOOLBOX
================

Basic execution of MToolBox can be run as follows:

	MToolBox.sh -i <input_format>

Please note that you MUST specify the -i option (sam|bam|fastq|fasta).
The command must be executed inside the folder containing your input files. Please note that only one of the supported formats can be used within a single MToolBox run.
For a complete list of MToolBox.sh options please run as follows:
	
	MToolBox.sh -h

IMPORTANT: if you wish to run again MToolBox in the same folder, it is preferable that you delete all the files produced during the previous execution.

MTOOLBOX OUTPUTS
================

MToolBox default outputs are:
- VCF_file.vcf contains all the mitochondrial variant positions against RSRS and other meta-information; 
- mt_classification_best_results.csv reports for each sequence the best haplogroup prediction. If the sorting results in more than one best haplogroup prediction with equal probability, the output will enclose all of them. 
- OUT_<sample_name> folders containing:
- outmt.sam: reads mapped onto the human mitochondrial DNA;
- logmt.txt: GSNAP log file for mitochondrial DNA mapping;
- outmt.fastq: fastq file of single reads extracted from outmt.sam file;
- outmt1.fastq: fastq file of paired reads extracted from outmt.sam file;
- outmt2.fastq: fastq file of paired reads extracted from outmt.sam file;
- outhumanS.sam: single reads mapped onto the entire human genome;
- loghumanS.txt: GSNAP log file for human genome mapping (single reads);
- outhumanP.sam: paired reads mapped onto the entire human genome;
- loghumanP.txt: GSNAP log file for human genome mapping (paired reads); - OUT.sam: alignments of reads uniquely mapped on mitochondrial genome;  
- OUT2.sam: alignments of reads uniquely mapped on mitochondrial genome, after processing with IndelRealigner and/or MarkDuplicates. This file will be generated anyway, even if these two processes have been disabled;
- mtDNAassembly-table.txt: the main table describing the assembly position by position;
- mtDNAassembly-Contigs.fasta: a fasta file including all reconstructed contigs;
- mtDNAassembly-coverage.txt: a text file including the coverage per contig and per mitochondrial known annotation;
-logassemble.txt, which is the log file of the assembleMTgenome.py script 
- sorted.csv contains a table with each haplogroup whose prediction is > 90%. It contains the following fields:
	
	N, the number of SNPs in the fragment sequence vs RSRS;
	Nph, the number of SNPs (among N) mapped in Phylotree;
	Nph_tot, the number of SNPs defining the haplogroup in the whole genome;
	Nph_exp, the number of SNPs defining the haplogroup in the fragment region;
	P_Hg, the prediction percentage value for the haplogroup (Nph/Nph_exp*100);
	Missing sites, the mutation events that are not present in the query genome but were expected from its respective path to the RSRS. These mutations may also point to a sequencing 	error;
- merged_diff.csv file reports the SNPs between the query genome and each of the three sequences RSRS, rCRS and hg_MHCS (Macro-Haplogroup Consensus Sequence);
- <sample_name>.csv contains a table where, for all the haplogroups present in the Phylotree Build 15, are reported the same data as in the file sequence_name>.sorted.csv, except for the Missing Sites field;
- annotation.csv is a further elaboration of the file merged_diff.csv, providing, for each mt variant allele between the query genome and each of the three sequences RSRS, rCRS and hg_MHCS, several annotations:
	 Nt Position:	Nucleotide position in mitochondrial genome
	Locus	Mitochondrial gene locus;
	Nt Variability: SiteVar variability value; 
	Codon Position: Nucleotide position within the codon;
	Aa Change: Amino Acid Change;
	Aa variability: MitVarProt amino acid variability value;
	% Disease: Percentage of predictor methods providing a mutation as 'Disease';
	% Neutral:Percentage of predictor methods providing a mutation as 'Neutral';
	% Unclassified: Percentage of predictor methods providing a mutation as 'Unclassified';
	MutPred Score:	MutPred Pathogenicity Score (0.000-1.000);
	PolyPhen-2 HumDiv Pred: Polyphen-2 HumDiv predictions (Benign, Possibly damaging, Probably damaging, Unknown);
	PolyPhen-2 HumDiv Prob: Polyphen-2 HumDiv probabilities (0.000-1.000);
	PolyPhen-2 HumVar Pred: Polyphen-2 HumVar predictions (Benign, Possibly damaging, Probably damaging, Unknown);
	PolyPhen-2 HumVar Prob: Polyphen-2 HumVar probabilities (0.000-1.000);
	PANTHER Pred:	PANTHER predictions (Neutral, Disease, Unclassified) by SNPs&GO software;
	PANTHER Prob:	PANTHER probabilities (0.000-1.000) by SNPs&GO software;
	PhD-SNP Pred:	PhD-SNP predictions (Neutral, Disease, Unclassified) by SNPs&GO software;
	PhD-SNP Prob:	PhD-SNP probabilities (0.000-1.000) by SNPs&GO software;
	SNPs&GO Pred:	SNPs&GO predictions (Neutral, Disease, Unclassified) by SNPs&GO software;
	SNPs&GO Prob:	SNPs&GO probabilities (0.000-1.000) by SNPs&GO software;
	Mitomap Associated Disease(s):	MITOMAP annotation of disease-associated mutations;
	Mitomap Homoplasmy:	MITOMAP annotation of homoplasmy condition;
	Mitomap Heteroplasmy: MITOMAP annotation of heteroplasmy condition;
	Somatic Mutations: MITOMAP annotation of cell or tissue type for somatic mutations;
	SM Homoplasmy:	MITOMAP annotation of homoplasmy condition in somatic mutations;
	SM Heteroplasmy:MITOMAP annotation of heteroplasmy condition in somatic mutations;
	OMIM Link to OMIM entry;
	Mamit-tRNA: Link to Mamit-tRNA site annotation;
	AC/AN: 1000 Genomes	Ratio between allele count and allele number of possibly pathogenic variants found in 1000 Genomes;
	1000 Genomes Homoplasmy: annotation of homoplasmy status in 1000 Genomes variants;
	1000 Genomes Heteroplasmy	: annotation of the heteroplasmy status in 1000 Genomes variants.

	WARNING! Please note that the heteroplasmy fractions and the related confidence interval will be reported only for variants found against the reference sequence chosen for read mapping.

NOTE ON FILE NAMES
==================

The basename for output folder and files will be parsed from the input filename, for each sample. 

BAM|SAM files: BAM or SAM files MUST be renamed as <sample_name>.ext, eg:

mv HG00096.chrom20.ILLUMINA.bwa.GBR.low_coverage.20101123.bam HG00096.bam

and "HG00096" will be the output basename.

FASTQ files: FASTQ files MUST be renamed as <sample_name>.R1.fastq, <sample_name>.R2.fastq for PAIRED-END data and <sample_name>.fastq for SINGLE END data.

IMPORTANT: Please note that MToolBox cannot recognize more than one PAIRED-END couple of fastq files (R1+R2) and one SINGLE-END fastq file per sample.


MTOOLBOX PROGRAMS OPTIONS
=========================

Besides the mandatory -i option, the execution of MToolBox can be fine-tuned by tweaking parameters of the mapExome.py and assembleMTgenome.py scripts by running MToolBox as follows:

	MToolBox.sh -i <input_format> -r <reference_sequence> -m "<mapExome_options>" -a "<assembleMTgenome_options>" -c "<mt-classifier_options>"


Most relevant options:

-i (MToolBox.sh) to choose the input file(s) format (mandatory);
-r (MToolBox.sh) to choose the mitochondrial reference sequence for read mapping. Please note that the selected reference sequence will be used as reference for the VCF output file. Allowed values are RCRS and RSRS. Default is RSRS.
-M (MToolBox.sh) to enable duplicate read removal by MarkDuplicates;
-I (MToolBox.sh) to enable mapped reads realignment around indels annotated in MITOMAP and HMTDB by GenomeAnalysisTK.jar;

-t (mapExome.py) to set the number of threads used by gsnap. Default is 8;

-t (assembleMTgenome.py) to set the minimum distance from the read end required to retain an indel for variant calling. Default is 5. Please note that only values â‰¥ 5 are allowed;
-z (assembleMTgenome.py) to set the minimum heteroplasmy threshold for variants to be reported in the FASTA consensus sequence. Default is 0.80;

For the full list of MToolBox options, please run as follows:

	MToolBox.sh -h

For the full list of mapExome and assembleMTgenome options, please run as follows:

mapExome.py -h
assembleMTgenome.py -h
mt-classifier.py -h

Usage example:

MToolBox.sh -M -I -m "-t 20 -M /path/to/gsnapdb/chrRSRS/ -H /path/to/gsnapdb/hg19RSRS/" -a "-t 10"

this command will enable duplicate read removal (-M), enable mapped reads realignment (-I); also, the number of threads used by gsnap will be set to 20 (-m "-t 20..") and the custom path for gsnap mitochondrial and nuclear indexes will be set (-M and -H respectively). Finally, the minimum nucleotide distance from read end to retain indels will be set to 10 (-a "-t 10").


CONTACTS AND CITATION
=====================

Contacts:
dome.simone@gmail.com
claudia.calabrese23@gmail.com

If you use MToolBox, please cite:
Calabrese C, Simone D, Diroma MA, Santorsola M, Guttˆ C, Gasparre G, Picardi E, Pesole G, Attimonelli M. MToolBox: a highly automated pipeline for heteroplasmy annotation and prioritization analysis of human mitochondrial variants in high-throughput sequencing. Bioinformatics. 2014 Jul 14. pii: btu483. Epub ahead of print] PubMed PMID: 25028726

To get support please 
1) Join the MToolBox google group:
https://groups.google.com/forum/?hl=IT#!forum/mtoolbox-users

or
2) Create a topic in the MToolBox General Discussion here:
https://sourceforge.net/p/mtoolbox/discussion/general/

CHANGELOG - 5 June 2014
==========================

- changed the method for estimation of the heteroplasmy confidence interval (CI). For sites with coverage depth <= 40, the heteroplasmy CI is estimated with the Wilson score interval; for larger coverage depth values, the Agresti-Coull interval is used.
- added the possibility to use fasta inputs to perform haplogroup prediction and functional annotation.
- added the possibility to use the revised Cambridge Reference Sequence (rCRS) as reference sequence for read mapping. By using rCRS as reference sequence, the VCF output will be rCRS-based. 

CHANGELOG - 19 July 2014
==========================

- an error encountered during the bam to fastq extraction has been fixed. Empty unpaired fastq files are now removed.  

CHANGELOG - 29 September 2014
==========================

- an error fixed in the -t parameter of assembleMTgenome.

CHANGELOG - 6 November 2014
==========================
Update to MToolBox version 0.2.1:

 - an error encountered during the sam to fastq extraction has been fixed in the MToolBox.sh file.
 - RCRS, hg19RCRS, RSRS and hg19RSRS gmap indexed databases have been regenerated using the -c option for circularized chromosomes.
 - update to Phylotree Build 16 for haplogroup prediction. 

CHANGELOG - 23 January 2015
==========================
Update to MToolBox version 0.2.2:

-an error encountered during the analysis of hard clipping mapped reads has been fixed in the mtVariantCaller.py
-the mtVariantCaller.py has been improved to better manage sites with multiple alleles.
-hidden files included in the package and generating problems with the mt-classifier.py have been eliminated.

CHANGELOG - 25 February 2015
==========================

-An error occurred during the generation of circularized mitochondrial chromosome, used for the gsnap db generation. Hg19RCRS/hg19RSRS and chrRCRS/chrRSRS gsnap indexed databases have been replaced with those using the linearized mitochondrial chromosome. We apologize with the MToolBox users for this inconvenient. 

CHANGELOG - 28 February 2015
==========================
Update to MToolBox version 0.3 with the following new options and changes:

- fastq.gz is a further possible input format file. Installation of zlib libraries is therefore required. 
- users can specify the path of the working directories using Ðp (path to input folder) and Ðo (path to output folder) options.
- users can specify a list of files to be used as input through -l option. It accepts a text file containing one sample name for each line. This list should be named as "list.txt" and placed in the input folder. Alternatively, users can provide comma-separated names with the same option. It is mandatory to report in such list the filename extension (e.g. mysample.sam or mysample.R1.fastq).
- users can use -X option to allow the extraction from a BAM file of mitochondrial reads mapped onto a mitochondrial reference sequence. This option can be useful when using Whole Genome or Exome sequencing BAM files containing a huge amount ofnuclear reads. The option works only with the BAM format.

new fields added to the annotation.csv output file:
- Disease Score: "% Disease", "% Neutral" and "% Unclassified" fields have been replaced with an overall Disease Score, generated as a weighted average of pathogenicity prediction scores for non-synonymous variants, derived from a training dataset of 53 non synonymous variants selected among mitochondrial diseases or cancer associated mutations. Weights have been calculated by taking into account the right prediction and the best probability to predict a truly pathogenic mtDNA variants generated by the pathogenicity prediction algorithms currently implemented in MToolBox.
- MutPred Pred: MutPred prediction (Low pathogenicity, High pathogenicity).
- dbSNP ID: Variant ID in dbSNP.

for users convenience, new scripts to help the generation of reports about the annotated and prioritized variants have been added to the suite of tools provided by MToolBox:
 - prioritization.py, which generates the prioritized_variants.txt file, reporting annotation only for prioritized variants for each sample analyzed, defined as variants recognized by the three reference sequences (rCRS, RSRS and MHCS), sorted per increasing variability.
 - summary.py, which generates the summary.txt file, reporting statistics about the coverage of reconstructed mitochondrial genomes, number of homoplasmic and heteroplasmic variants (for NGS data), haplogroup prediction and number of prioritized variants.

CHANGELOG - 8 June 2015
==========================
Update to MToolBox version 0.3.1 with the following change:

- A bug in the MToolBox.sh script has been fixed. Empty paired-end fastq files generated during SAM/BAM to fastq conversion are now removed.
