MToolBox - README
=================

MToolBox is a highly automated bioinformatics pipeline to reconstruct and analyze human mitochondrial DNA from high throughput sequencing data. MToolBox includes an updated computational strategy to assemble mitochondrial genomes from Whole Exome and/or Genome Sequencing (PMID: 22669646) and an improved fragment-classify tool (PMID:22139932) for haplogroup assignment, functional and prioritization analysis of mitochondrial variants. MToolBox provides pathogenicity scores, profiles of genome variability and disease-associations for mitochondrial variants. MToolBox provides also a Variant Call Format file (version 4.0) featuring, for the first time, allele-specific heteroplasmy.
The MToolBox pipeline includes:

- an extended version of a previously published computational strategy for mtDNA genome assembly (PMID: 22669646). The pipeline has been integrated with the detection of insertions and deletions (indels) and the assessment of the heteroplasmic fraction (HF) and related confidence interval (CI) for each mt variant. HF and CI are integrated as genotype specific meta-information in a Variant Call Format (version 4.0) file;
- the mt-classifier tool, for haplogroup prediction, mt variant functional annotation and prioritization


SYSTEM REQUIREMENTS
===================

- UNIX-based OS
- Python2.7 (www.python.org)
- samtools (https://sourceforge.net/projects/samtools/files/samtools/) installed in /usr/local/bin/samtools (otherwise you must specify the actual path using assembleMTgenome.py options). Please be aware that the current version of MToolBox does not support samtools versions >= 1.0. 
- MUSCLE (http://www.drive5.com/muscle/downloads.htm) installed in /usr/local/bin/muscle (otherwise you should specify the actual path using mt-classifier.py options)
- GSNAP (https://github.com/julian-gehring/GMAP-GSNAP; newest version of GMAP-GSNAP at http://research-pub.gene.com/gmap/) installed in /usr/local/bin/gsnap (otherwise you must specify the actual path using mapExome.py options). 
  WARNING: The GSNAP indexes provided in the sourceforge page of MToolBox (http://sourceforge.net/projects/mtoolbox/files/genome_index/) have been generated with GSNAP version 2015-12-31.v7. We cannot guarantee that such indexes will work with later or previous versions of GSNAP and we strongly recommend to regenerate the indexes with the aligner version the user may want to choose.

WARNING: due to the file size of the nuclear GSNAP databases generated with the GSNAP version 2015-12-31.v7 (~20GB) we were unable to upload the hg19RCRS/hg19RSRS databases to the MToolBox sourceforge webpage at https://sourceforge.net/projects/mtoolbox/files/genome_index/. MToolBox users may follow the guidelines to build GSNAP databases available at https://github.com/mitoNGS/MToolBox/blob/master/how_to_build_db.md.
Further information about how to to generate GMAP/GSNAP database can be found at the http://research-pub.gene.com/gmap/src/README

MTOOLBOX SCRIPTS
================

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
- variants_functional_annotation.py for functional annotation;
- prioritization.py for variant prioritization;
- summary.py to report statistics about coverage of reconstructed genomes, predicted haplogroups, number of homoplasmic, heteroplasmic and prioritized variants


OTHER FILES REQUIRED 
====================

By default, MToolBox adopts the RSRS (Reconstructed Sapiens Reference Sequence, PMID: 22482806) as mitochondrial reference genome and hg19 as nuclear reference genome. Alternatively, the user can choose to use the rCRS (revised Cambridge Reference Sequence).
For user's convenience, fasta files and related gsnap indexes for the two reference sequences can be downloaded at https://sourceforge.net/projects/mtoolbox/.

List of fasta files and index .fai files in the genome_fasta folder (https://sourceforge.net/projects/mtoolbox/files/genome_fasta/) :
- chrRSRS.fa
- chrRSRS.fa.fai
- hg19RSRS.fa
- hg19RSRS.fa.fai
- chrM.fa
- chrM.fa.fai
- hg19RCRS.fa
- hg19RCRS.fa.fai

The default directory for these files is /usr/local/share/genomes/, otherwise you must specify the actual path using the assembleMTgenome.py options -r (for further details please refer to "RUNNING MTOOLBOX" section). 

List of compressed files of gsnap indexes in the genome_index folder (https://sourceforge.net/projects/mtoolbox/files/genome_index/):
- chrRSRS.tar.gz
- hg19RSRS.tar.gz [currently unavailable]
- chrM.tar.gz
- hg19RCRS.tar.gz [currently anavailable]

The default directory for decompressed folders is /usr/local/share/gmapdb/, otherwise you must specify the actual path using the mapExome.py options -M and -H, respectively (for further details please refer to "RUNNING MTOOLBOX" section).

The MToolBox folder includes the MITOMAP_HMTDB_known_indels.vcf file, containing 127 known indels annotated in MITOMAP and HMTDB and the related intervals_file.list used by GATK GenomeAnalysisTK.jar module.
The MToolBox folder also includes 2 tab-separated files, variant_info.txt (previously called patho_table.txt) and site_info.txt (previously called sitevar_modified.txt) containing variant-specific and site-specific information, respectively, used in the annotation step.


INSTALLATION
============

Simply add the MToolBox path to your system PATH with the following command:

export PATH="/path/to/MToolBox/:$PATH"


RUNNING MTOOLBOX
================
Basic execution of MToolBox can be run as follows:

MToolBox.sh -i <input_format>


Please note that you MUST specify the -i option (sam|bam|fastq|fasta).
The basic command must be executed inside the folder containing your input files, otherwise see MTOOLBOX PROGRAM OPTIONS section to specify input/output folders paths. Please note that only one of the supported formats can be used within a single MToolBox run. MToolBox also accepts compressed fastq files (FASTQ.GZ). 
For a complete list of MToolBox.sh options please run as follows:

	
MToolBox.sh -h


IMPORTANT: if you wish to run again MToolBox in the same folder, it is preferable that you delete all the files produced during the previous execution.


MTOOLBOX PROGRAMS OPTIONS
=========================

Besides the mandatory -i option, the execution of MToolBox can be fine-tuned by tweaking parameters of the mapExome.py and assembleMTgenome.py scripts by running MToolBox as follows:

MToolBox.sh -i <input_format> -r <reference_sequence> -m "<mapExome_options>" -a "<assembleMTgenome_options>" -c "<mt-classifier_options>"


Most relevant options:

-i (MToolBox.sh) to choose the input file(s) format (mandatory).

-p (MToolBox.sh) to set the specific path to the input folder.

-o (MToolBox.sh) to set the specific path to the output folder.

-l (MToolBox.sh) . Indicate a list of file names to be analyzed. Argument of the option can be a list of comma separated file names or a txt file with one file name per line. 

-X (MToolBox.sh) to enable extraction of mitochondrial reads from the input file, avoiding realignment of all the reads. Useful only for BAM input files

-r (MToolBox.sh) to choose the mitochondrial reference sequence for read mapping. Please note that the selected reference sequence will be used as reference for the VCF output file. Allowed values are RCRS and RSRS. Default is RSRS.

-M (MToolBox.sh) to enable duplicate read removal by MarkDuplicates.

-I (MToolBox.sh) to enable mapped reads realignment around indels annotated in MITOMAP and HMTDB by GenomeAnalysisTK.jar.

-t (mapExome.py) to set the number of threads used by gsnap. Default is 8.

-t (assembleMTgenome.py) to set the minimum distance from the read end required to retain an indel for variant calling. Default is 5. Please note that only values >= 5 are allowed.

-z (assembleMTgenome.py) to set the minimum heteroplasmy threshold for variants to be reported in the FASTA consensus sequence. Default is 0.80.

For the full list of MToolBox options, please run as follows:

MToolBox.sh -h

For the full list of mapExome and assembleMTgenome options, please run as follows:

mapExome.py -h
assembleMTgenome.py -h
mt-classifier.py -h

Usage examples:

MToolBox.sh -i fastq -M -I -m "-t 20 -M /path/to/gsnapdb/chrRSRS/ -H /path/to/gsnapdb/hg19RSRS/" -a "-t 10"

this command, launched in the input folder, will take all the fastq files as input (-i fastq), enable duplicate read removal (-M), enable mapped reads realignment (-I); also, the number of threads used by gsnap will be set to 20 (-m "-t 20..") and the custom path for gsnap mitochondrial and nuclear indexes will be set (-M and -H respectively). Finally, the minimum nucleotide distance from read end to retain indels will be set to 10 (-a "-t 10").

MToolBox.sh -i bam -l mysample1.bam,mysample2.bam -p /path/to/input/folder/ -X -a "-z 0.9" -o /path/to/output/folder/

this command will analyze 2 bam files (-i bam -l mysample1.bam,mysample2.bam) in the input folder (-p /path/to/input/folder/). Only reads previously mapped to the mitochondrial genome will be considered in the analysis (-X), avoiding the re-mapping of all the reads cointained in the input files. Only variants with heteroplasmic fraction HF>0.9 will be reported in the FASTA consensus sequence (-a "-z 0.9"). Finally, all the results will be written in the specified output folder (-o /path/to/output/folder/).


MTOOLBOX OUTPUTS
================

MToolBox default outputs are:
- VCF_file.vcf contains all the mitochondrial variant positions against RSRS and other meta-information. 
- mt_classification_best_results.csv reports for each sequence the best haplogroup prediction. If the sorting results in more than one best haplogroup prediction with equal probability, the output will enclose all of them.
- prioritized_variants.txt contains annotation only for prioritized variants for each sample analyzed, defined as variants recognized by the three reference sequences (rCRS, RSRS and MHCS), sorted per increasing nucleotide variability.
- summary_<date_time>.txt reporting a brief summary of selected options, predicted haplogroups, number of total and prioritized variants for each sample and, for NGS data only, coverage of reconstructed genomes, number of homoplasmic and heteroplasmic variants.
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
- logassemble.txt, which is the log file of the assembleMTgenome.py script 
- sorted.csv contains a table with each haplogroup whose prediction is > 90%. It contains the following fields:
	1. N = the number of SNPs in the fragment sequence vs RSRS;
	2. Nph = the number of SNPs (among N) mapped in Phylotree;
	3. Nph_tot = the number of SNPs defining the haplogroup in the whole genome;
	4. Nph_exp = the number of SNPs defining the haplogroup in the fragment region;
	5. P_Hg = the prediction percentage value for the haplogroup (Nph/Nph_exp*100);
	6. Missing sites = the mutation events that are not present in the query genome but were expected from its respective path to the RSRS. These mutations may also point to a sequencing error;
- merged\_diff.csv file reports the SNPs between the query genome and each of the three sequences RSRS, rCRS and hg\_MHCS (Macro-Haplogroup Consensus Sequence);
- \<sample\_name>.csv contains a table where, for all the haplogroups present in the Phylotree Build 15, are reported the same data as in the file \<sequence\_name>.sorted.csv, except for the Missing Sites field;
- annotation.csv is a further elaboration of the file \merged\_diff.csv, providing, for each mt variant allele between the query genome and each of the three sequences RSRS, rCRS and hg\_MHCS, several annotations:
	1. Sample = sample name;
	2. Variant Allele = nucleotide position in mitochondrial genome followed by the variant allele;
	3. HF = Heteroplasmic Fraction, as reported in the VCF file;
	4. CI_lower;CI_upper = lower and upper limits of the confidence interval of the heteroplasmic fraction;
	5. RSRS = if "yes", the variant is recognized by RSRS;
	6. MHCS = if "yes", the variant is recognized by the Macro-Haplogroup Consensus Sequence;
	7. rCRS = if "yes", the variant is recognized by rCRS;
	8. Haplogroup = the best predicted haplogroup;
	9. Other Haplogroups = if "+", the variant defines other haplogroups beside the sample specific haplogroup;
	10. Locus = mitochondrial gene locus;
	11. Nt Variability = SiteVar variability value;
	12. Codon Position = nucleotide position within the codon;
	13. Aa Change = amino Acid Change;
	14. Aa variability = MitVarProt amino acid variability value;
	15. tRNA annotation = specific information regarding mitochondrial tRNA genes (position in tRNA; tRNA type; cloverleaf secondary region; mature nucleotide; involvement of the specific position in tRNA folding);
	16. Disease score = an overall Disease Score, generated as a weighted average of pathogenicity prediction scores for non-synonymous variants. For details, see the related publication (PMID: 26621530);
	17. RNA predictions = score added for 49 variants in rRNA genes (PMID: 24092330) and 207 variants in tRNA genes (PMID: 21882289; PMID:23696415). Scores were correlated on a scale from 0 to 1. Threshold for rRNAs=0.51. Threshold for tRNAs= 0.31. Low pathogenicity under the fixed thresholds;
	18. MutPred pred = MutPred predictions (High pathogenicity, Low pathogenicity);
	19. MutPred Score = MutPred Pathogenicity Score (0.000-1.000);
	20. PolyPhen-2 HumDiv Pred = Polyphen-2 HumDiv predictions (Benign, Possibly damaging, Probably damaging, Unknown);
	21. PolyPhen-2 HumDiv Prob = Polyphen-2 HumDiv probabilities (0.000-1.000);
	22. PolyPhen-2 HumVar Pred = Polyphen-2 HumVar predictions (Benign, Possibly damaging, Probably damaging, Unknown);
	23. PolyPhen-2 HumVar Prob = Polyphen-2 HumVar probabilities (0.000-1.000);
	24. PANTHER Pred = PANTHER predictions (Neutral, Disease, Unclassified) by SNPs&GO software;
	25. PANTHER Prob = PANTHER probabilities (0.000-1.000) by SNPs&GO software;
	26. PhD-SNP Pred = PhD-SNP predictions (Neutral, Disease, Unclassified) by SNPs&GO software;
	27. PhD-SNP Prob = PhD-SNP probabilities (0.000-1.000) by SNPs&GO software;
	28. SNPs&GO Pred = SNPs&GO predictions (Neutral, Disease, Unclassified) by SNPs&GO software;
	29. SNPs&GO Prob = SNPs&GO probabilities (0.000-1.000) by SNPs&GO software;
	30. Mitomap Associated Disease(s) = MITOMAP annotation of disease-associated mutations;
	31. Mitomap Homoplasmy = MITOMAP annotation of homoplasmy condition;
	32. Mitomap Heteroplasmy = MITOMAP annotation of heteroplasmy condition;
	33. Somatic Mutations = MITOMAP annotation of cell or tissue type for somatic mutations;
	34. SM Homoplasmy = MITOMAP annotation of homoplasmy condition in somatic mutations;
	35. SM Heteroplasmy = MITOMAP annotation of heteroplasmy condition in somatic mutations;
	36. ClinVar = ClinVar annotation of associated disease(s);
	37. OMIM Link = link to OMIM entry;
	38. dbSNP ID = rs ID reported in dbSNP database;
	39. Mamit-tRNA link = link to Mamit-tRNA site annotation;
	40. PhastCons20Way = PhastCons conservation score calculated on 20 vertebrates using hg38+rCRS as reference sequence;
	41. PhyloP20Way = PhyloP conservation score calculated on 20 vertebrates using hg38+rCRS as reference sequence;
	42. AC/AN 1000 Genomes = Ratio between allele count and allele number of possibly pathogenic variants found in 1000 Genomes using MToolBox;
	43. 1000 Genomes Homoplasmy = annotation of homoplasmy status in 1000 Genomes variants;
	44. 1000 Genomes Heteroplasmy = annotation of the heteroplasmy status in 1000 Genomes variants.

	WARNING! Please note that the heteroplasmy fractions and the related confidence interval will be reported only for variants found against the reference sequence chosen for read mapping.


NOTE ON FILE NAMES
==================

The basename for output folder and files will be parsed from the input filename, for each sample. 

BAM|SAM files: BAM or SAM files MUST be renamed as \<sample\_name\>.ext, eg:

mv HG00096.chrom20.ILLUMINA.bwa.GBR.low_coverage.20101123.bam HG00096.bam

and "HG00096" will be the output basename.

FASTQ files: FASTQ files MUST be renamed as \<sample\_name\>.R1.fastq, \<sample\_name\>.R2.fastq for PAIRED-END data and \<sample\_name\>.fastq for SINGLE END data. FASTQ compressed input files could be accepted with *.fastq.gz extension.

IMPORTANT: Please note that MToolBox cannot recognize more than one PAIRED-END couple of fastq files (R1+R2) and one SINGLE-END fastq file per sample.


MTOOLBOX GUI
============

For small datasets, you could use MToolBox at MSeqDR website https://mseqdr.org/mtoolbox.php (PMID: 25542617).


CONTACTS AND CITATION
=====================

Contacts:
dome.simone [at] gmail.com ; 
claudia.calabrese23 [at] gmail.com


If you use MToolBox, please cite:
Calabrese C, Simone D, Diroma MA, Santorsola M, GuttÃ  C, Gasparre G, Picardi E, Pesole G, Attimonelli M. MToolBox: a highly automated pipeline for heteroplasmy annotation and prioritization analysis of human mitochondrial variants in high-throughput sequencing. Bioinformatics. 2014 Nov 1;30(21):3115-7. doi: 10.1093/bioinformatics/btu483. Epub 2014 Jul 14. PubMed PMID: 25028726; PubMed Central PMCID: PMC4201154.

If you use MToolBox GUI, please also cite:

Falk MJ, Shen L, Gonzalez M, Leipzig J, Lott MT, Stassen AP, Diroma MA, Navarro-Gomez D, Yeske P, Bai R, Boles RG, Brilhante V, Ralph D, DaRe JT, Shelton R, Terry SF, Zhang Z, Copeland WC, van Oven M, Prokisch H, Wallace DC, Attimonelli M, Krotoski D, Zuchner S, Gai X; MSeqDR Consortium Participants; MSeqDR Consortium participants. Mitochondrial Disease Sequence Data Resource (MSeqDR): a global grass-roots consortium to facilitate deposition, curation, annotation, and integrated analysis of genomic data for the mitochondrial disease clinical and research communities. Mol Genet Metab. 2015 Mar;114(3):388-96. doi: 10.1016/j.ymgme.2014.11.016. Epub 2014 Dec 4. Review. PubMed PMID: 25542617; PubMed Central PMCID: PMC4512182.

If you would like to refer to our variant prioritization method, please also cite:

Santorsola M, Calabrese C, Girolimetti G, Diroma MA, Gasparre G, Attimonelli M. A multi-parametric workflow for the prioritization of mitochondrial DNA variants of clinical interest. Hum Genet. 2016 Jan;135(1):121-36. doi:10.1007/s00439-015-1615-9. Epub 2015 Nov 30. PubMed PMID: 26621530; PubMed Central PMCID: PMC4698288.

Get support 
===========

Please join the MToolBox google group:
https://groups.google.com/forum/?hl=IT#!forum/mtoolbox-users


CHANGELOGS
==========

March 1, 2016
=================
Update to MToolBox version 0.3.3 with the following change:

- GSNAP databases available at
  https://sourceforge.net/projects/mtoolbox/files/genome_index/
 generated with the GSNAP version 2013-09-11 have been removed, due to an
inconsistency betw
een the mitochodrial reference rCRS sequence fasta and index database name.
New GSNAP nuclea
r-mitochondrial and mitochondrial databases have been uploaded, generated with
the GSNAP ver
sion 2015-12-31.v7. Please, be aware that these databases might not be
compatible with previ
ous GSNAP versions.

-  Mitochondrial rCRS reference fasta and GSNAP database default names used in
   the MToolBox.
sh script were changed as following:

chrRCRS.fa --> chrM.fa
chrRCRS (GSNAP db) --> chrM

to reflect changes in GSNAP databases uploaded at
https://sourceforge.net/projects/mtoolbox/
files/genome_index/.

September 21, 2015
==================

Update to MToolBox version 0.3.1a with the following change:

 - A bug in patho-table.txt has been fixed. 313 new stop-gain mutations and 6 new missense variants are now included.

New fields added to the annotation.csv output file:
 - tRNA annotation: specific information regarding mitochondrial tRNA genes (position in tRNA; tRNA type; cloverleaf secondary region; mature nucleotide; involvement of the specific position in tRNA folding)
 - RNA predictions: score added for 49 variants in rRNA genes (Smith PM et al, 2014, PMID:24092330) and 207 variants in tRNA genes (Yarham JW et al, 2011, PMID:21882289; Blakely EL et al, 2013, PMID:23696415). Scores were retrieved from literature and correlated on a scale from 0 to 1. Threshold for rRNAs=0.51. Threshold for tRNAs= 0.31. Low pathogenicity under the fixed thresholds.
 - ClinVar: ClinVar annotation of associated disease(s) (January 21, 2015 update)
 - PhastCons20Way: PhastCons conservation score calculated on 20 vertebrates using hg38+rCRS as reference sequence
 - PhyloP20Way: PhyloP conservation score calculated on 20 vertebrates using hg38+rCRS as reference sequence

Fields updated in the annotation.csv output file:
 - Nt variability: SiteVar variability value calculated on 22,691 complete healthy genomes in HmtDB database (May 2015 update)
 - Aa variability: MitVarProt variability value calculated on 22,691 complete healthy genomes in HmtDB database (May 2015 update)
 - Mitomap associated disease(s), Mitomap Homoplasmy, Mitomap Heteroplasmy: July 20, 2015 update
 - Mitomap somatic mutations, SM Homoplasmy, SM Heteroplasmy and Mitomap associated disease(s) only RNA mutations: July 29, 2015 update
 - dbSNP ID: release 144, May 26, 2015
 - OMIM link: August 4, 2015 update


June 8, 2015
============

Update to MToolBox version 0.3.1 with the following change:

- A bug in the MToolBox.sh script has been fixed. Empty paired-end fastq files generated during SAM/BAM to fastq conversion are now removed.


February 28, 2015
=================

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


February 25, 2015
=================

 - An error occurred during the generation of circularized mitochondrial chromosome, used for the gsnap db generation. Hg19RCRS/hg19RSRS and chrRCRS/chrRSRS gsnap indexed databases have been replaced with those using the linearized mitochondrial chromosome. We apologize with the MToolBox users for this inconvenient. 


January 23, 2015
================

Update to MToolBox version 0.2.2:

 - an error encountered during the analysis of hard clipping mapped reads has been fixed in the mtVariantCaller.py
 - the mtVariantCaller.py has been improved to better manage sites with multiple alleles.
 - hidden files included in the package and generating problems with the mt-classifier.py have been eliminated.


November 6, 2014
================

Update to MToolBox version 0.2.1:

 - an error encountered during the sam to fastq extraction has been fixed in the MToolBox.sh file.
 - RCRS, hg19RCRS, RSRS and hg19RSRS gmap indexed databases have been regenerated using the -c option for circularized chromosomes.
 - update to Phylotree Build 16 for haplogroup prediction. 


September 29, 2014
==================

- an error fixed in the -t parameter of assembleMTgenome.


July 19, 2014
=============

- an error encountered during the bam to fastq extraction has been fixed. Empty unpaired fastq files are now removed.  


June 5, 2014
============

- changed the method for estimation of the heteroplasmy confidence interval (CI). For sites with coverage depth <= 40, the heteroplasmy CI is estimated with the Wilson score interval; for larger coverage depth values, the Agresti-Coull interval is used.
- added the possibility to use fasta inputs to perform haplogroup prediction and functional annotation.
- added the possibility to use the revised Cambridge Reference Sequence (rCRS) as reference sequence for read mapping. By using rCRS as reference sequence, the VCF output will be rCRS-based. 
