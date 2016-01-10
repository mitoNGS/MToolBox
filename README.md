###MToolBox 
MToolBox is a highly automated bioinformatics pipeline to reconstruct and analyze human mitochondrial DNA from high throughput sequencing data. MToolBox includes an updated computational strategy to assemble mitochondrial genomes from Whole Exome and/or Genome Sequencing (PMID: 22669646) and an improved fragment-classify tool (PMID:22139932) for haplogroup assignment, functional and prioritization analysis of mitochondrial variants. MToolBox provides pathogenicity scores, profiles of genome variability and disease-associations for mitochondrial variants. MToolBox provides also a Variant Call Format file (version 4.0) featuring, for the first time, allele-specific heteroplasmy.
The MToolBox pipeline includes:

- an extended version of a previously published computational strategy for mtDNA genome assembly (PMID: 22669646). The pipeline has been integrated with the detection of insertions and deletions (indels) and the assessment of the heteroplasmic fraction (HF) and related confidence interval (CI) for each mt variant. HF and CI are integrated as genotype specific meta-information in a Variant Call Format (version 4.0) file;
- the mt-classifier tool, for haplogroup prediction, mt variant functional annotation and prioritization


###SYSTEM REQUIREMENTS

- UNIX-based OS
- Python2.7 (www.python.org)
- GSNAP (https://github.com/julian-gehring/GMAP-GSNAP; newest version of GMAP-GSNAP at http://research-pub.gene.com/gmap/) installed in */usr/local/bin/gsnap* (otherwise you must specify the actual path using mapExome.py options). 
  **WARNING: The GSNAP indexes provided in the sourceforge page of MToolBox (http://sourceforge.net/projects/mtoolbox/files/genome_index/) have been generated with GSNAP version 2013-09-11. We cannot guarantee that such indexes will work with later versions of GSNAP and we strongly reccomend to regenerate the indexes with the aligner version the user may want to choose**.

  To generate GMAP/GSNAP database you can visit this link http://research-pub.gene.com/gmap/src/README or use the *build_gsnap_index.sh* script provided within the mitoNGS repo, to help the MToolBox users building GMAP/GSNAP databases. To check the usage of the script please type:
  ```
  build_gsnap_index.sh -h
 ```
  Usage example to build the GMAP/GSNAP database in the current working directory:
  ```
  build_gsnap_index.sh -D . -n hg19RCRS -f hg19RCRS.fa -k 12 
  ```

- samtools (https://sourceforge.net/projects/samtools/files/samtools/) installed in /usr/local/bin/samtools (otherwise you must specify the actual path using assembleMTgenome.py options)
- MUSCLE (http://www.drive5.com/muscle/downloads.htm) installed in */usr/local/bin/muscle* (otherwise you should specify the actual path using mt-classifier.py options)


###MTOOLBOX SCRIPTS

*MToolBox.sh* is the shell script invoking all the following programs:
- *Sam2Fastq.jar* to convert either SAM or BAM files to FASTQ format. This module is included in the ext_tools folder;
- *mapExome.py* to map reads on RSRS and hg19. It invokes GSNAP;
- *SortSam.jar*, MarkDuplicates.jar, SamFormatConverter.jar (PicardTools suite; http://picard.sourceforge.net) for SAM/BAM manipulation at several stages of the pipeline and elimination of PCR duplicates. These modules are included in the ext_tools folder;
- *GenomeAnalysisTK.jar* (GATK suite; PMID: 20644199) for indels realignment. This module is included in the ext_tools folder;
- *assembleMTgenome.py* to assemble the mitochondrial genome and perform variant calling and heteroplasmy quantification, by invoking:
- *mpileup* (samtools)
- *mtVariantCaller.py* module;
- *VCFoutput.py* to report variants in the VCF file (version 4.0). It invokes the vcf module (elease 0.60 PyVCF; https://githubcom/jamescasbon/PyVCF/);
- *mt-classifier.py* for haplogroup prediction;
- *variants_functional_annotation.py* for functional annotation and prioritization.


###OTHER FILES REQUIRED 

By default, MToolBox adopts the RSRS (Reconstructed Sapiens Reference Sequence, PMID: 22482806) as mitochondrial reference genome and hg19 as nuclear reference genome. Alternatively, the user can choose to use the rCRS (revised Cambridge Reference Sequence).
For user's convenience, fasta files and related gsnap indexes for the two reference sequences can be downloaded at https://sourceforge.net/projects/mtoolbox/.

List of fasta files and index .fai files in the genome_fasta folder (https://sourceforge.net/projects/mtoolbox/files/genome_fasta/) :
- *chrRSRS.fa*
- *chrRSRS.fa.fai*
- *hg19RSRS.fa*
- *hg19RSRS.fa.fai*
- *chrRCRS.fa*
- *chrRCRS.fa.fai*
- *hg19RCRS.fa*
- *hg19RCRS.fa.fai*

The default directory for these files is */usr/local/share/genomes/*, otherwise you must specify the actual path using the assembleMTgenome.py options ```-r``` (for further details please refer to "RUNNING MTOOLBOX" section). 

List of compressed files of gsnap indexes in the genome_index folder (https://sourceforge.net/projects/mtoolbox/files/genome_index/):
- *chrRSRS.tar.gz*
- *hg19RSRS.tar.gz*

The default directory for decompressed folders is */usr/local/share/gmapdb/*, otherwise you must specify the actual path using the *mapExome.py* options ```-M``` and ```-H```, respectively (for further details please refer to "RUNNING MTOOLBOX" section).

The MToolBox folder includes the *MITOMAP_HMTDB_known_indels.vcf* file, containing 127 known indels annotated in MITOMAP and HMTDB and the related *intervals_file.list* used by GATK *GenomeAnalysisTK.jar* module.


###INSTALLATION

Simply add the MToolBox path to your system PATH with the following command:
```
export PATH="/path/to/MToolBox/:$PATH"
```

###RUNNING MTOOLBOX

Basic execution of MToolBox can be run as follows:
```
MToolBox.sh -i <input_format>
```

Please note that you MUST specify the ```-i``` option (sam|bam|fastq|fasta).
The command must be executed inside the folder containing your input files. Please note that only one of the supported formats can be used within a single MToolBox run.
For a complete list of *MToolBox.sh* options please run as follows:

```	
MToolBox.sh -h
```

**IMPORTANT**: if you wish to run again MToolBox in the same folder, it is preferable that you delete all the files produced during the previous execution.

###MTOOLBOX OUTPUTS

MToolBox default outputs are:
- *VCF_file.vcf* contains all the mitochondrial variant positions against RSRS and other meta-information; 
- *mt_classification_best_results.csv* reports for each sequence the best haplogroup prediction. If the sorting results in more than one best haplogroup prediction with equal probability, the output will enclose all of them. 
- *OUT_<sample_name>* folders containing:
- *outmt.sam*: reads mapped onto the human mitochondrial DNA;
- *logmt.txt*: GSNAP log file for mitochondrial DNA mapping;
- *outmt.fastq*: fastq file of single reads extracted from outmt.sam file;
- *outmt1.fastq*: fastq file of paired reads extracted from outmt.sam file;
- *outmt2.fastq*: fastq file of paired reads extracted from outmt.sam file;
- *outhumanS.sam*: single reads mapped onto the entire human genome;
- *loghumanS.txt*: GSNAP log file for human genome mapping (single reads);
- *outhumanP.sam*: paired reads mapped onto the entire human genome;
- *loghumanP.txt*: GSNAP log file for human genome mapping (paired reads); - OUT.sam: alignments of reads uniquely mapped on mitochondrial genome;  
- *OUT2.sam*: alignments of reads uniquely mapped on mitochondrial genome, after processing with IndelRealigner and/or MarkDuplicates. This file will be generated anyway, even if these two processes have been disabled;
- *mtDNAassembly-table.txt*: the main table describing the assembly position by position;
- *mtDNAassembly-Contigs.fasta*: a fasta file including all reconstructed contigs;
- *mtDNAassembly-coverage.txt*: a text file including the coverage per contig and per mitochondrial known annotation;
- *logassemble.txt*, which is the log file of the assembleMTgenome.py script 
- *sorted.csv* contains a table with each haplogroup whose prediction is > 90%. It contains the following fields:
	1. *N* = the number of SNPs in the fragment sequence vs RSRS;
	2. *Nph* = the number of SNPs (among N) mapped in Phylotree;
	3. *Nph_tot* = the number of SNPs defining the haplogroup in the whole genome;
	4. *Nph_exp* = the number of SNPs defining the haplogroup in the fragment region;
	5. *P_Hg* = the prediction percentage value for the haplogroup (Nph/Nph_exp*100);
	6. *Missing sites* = the mutation events that are not present in the query genome but were expected from its respective path to the RSRS. These mutations may also point to a sequencing error;
- *merged\_diff.csv* file reports the SNPs between the query genome and each of the three sequences RSRS, rCRS and hg\_MHCS (*Macro-Haplogroup Consensus Sequence*);
- *\<sample\_name>.csv* contains a table where, for all the haplogroups present in the Phylotree Build 15, are reported the same data as in the file *\<sequence\_name>.sorted.csv*, except for the *Missing Sites* field;
- *annotation.csv* is a further elaboration of the file *\merged\_diff.csv*, providing, for each mt variant allele between the query genome and each of the three sequences RSRS, rCRS and hg\_MHCS, several annotations:
	1. *Nt Position* = the nucleotide position in mitochondrial genome
	2. *Locus* = Mitochondrial gene locus;
	3. *Nt Variability* = SiteVar variability value; 
	4. *Codon Position* = Nucleotide position within the codon;
	5. *Aa Change* = Amino Acid Change;
	6. *Aa variability* = MitVarProt amino acid variability value;
	7. *% Disease* = Percentage of predictor methods providing a mutation as 'Disease';
	8. *% Neutra* =  Percentage of predictor methods providing a mutation as 'Neutral';
	9. *% Unclassified* = Percentage of predictor methods providing a mutation as 'Unclassified';
	10. *MutPred Score* =	MutPred Pathogenicity Score (0.000-1.000);
	11. *PolyPhen-2 HumDiv Pred* = Polyphen-2 HumDiv predictions (Benign, Possibly damaging, Probably damaging, Unknown);
	12. *PolyPhen-2 HumDiv Prob* = Polyphen-2 HumDiv probabilities (0.000-1.000);
	13. *PolyPhen-2 HumVar Pred* = Polyphen-2 HumVar predictions (Benign, Possibly damaging, Probably damaging, Unknown);
	14. *PolyPhen-2 HumVar Prob* = Polyphen-2 HumVar probabilities (0.000-1.000);
	15. *PANTHER Pred* =	PANTHER predictions (Neutral, Disease, Unclassified) by SNPs&GO software;
	16. *PANTHER Prob* = PANTHER probabilities (0.000-1.000) by SNPs&GO software;
	17. *PhD-SNP Pred* = PhD-SNP predictions (Neutral, Disease, Unclassified) by SNPs&GO software;
	18. *PhD-SNP Prob* = PhD-SNP probabilities (0.000-1.000) by SNPs&GO software;
	19. *SNPs&GO Pred* = SNPs&GO predictions (Neutral, Disease, Unclassified) by SNPs&GO software;
	20. *SNPs&GO Prob* = SNPs&GO probabilities (0.000-1.000) by SNPs&GO software;
	21. *Mitomap Associated Disease(s)* = MITOMAP annotation of disease-associated mutations;
	22. *Mitomap Homoplasmy* = MITOMAP annotation of homoplasmy condition;
	23. *Mitomap Heteroplasmy* = MITOMAP annotation of heteroplasmy condition;
	24. *Somatic Mutations* = MITOMAP annotation of cell or tissue type for somatic mutations;
	25. *SM Homoplasmy* = MITOMAP annotation of homoplasmy condition in somatic mutations;
	26. *SM Heteroplasmy* = MITOMAP annotation of heteroplasmy condition in somatic mutations;
	26. *OMIM Link to OMIM entry*;
	26. *Mamit-tRNA* = Link to Mamit-tRNA site annotation;
	27. *AC/AN* = 1000 Genomes	Ratio between allele count and allele number of possibly pathogenic variants found in 1000 Genomes;
	28. *1000 Genomes Homoplasmy* = annotation of homoplasmy status in 1000 Genomes variants;
	29. *1000 Genomes Heteroplasmy* = annotation of the heteroplasmy status in 1000 Genomes variants.

	**WARNING**! Please note that the heteroplasmy fractions and the related confidence interval will be reported only for variants found against the reference sequence chosen for read mapping.

###NOTE ON FILE NAMES

The basename for output folder and files will be parsed from the input filename, for each sample. 

BAM|SAM files: BAM or SAM files MUST be renamed as *\<sample\_name\>.ext*, eg:
```
mv HG00096.chrom20.ILLUMINA.bwa.GBR.low_coverage.20101123.bam HG00096.bam
```
and "HG00096" will be the output basename.

FASTQ files: FASTQ files MUST be renamed as \<sample\_name\>.R1.fastq, \<sample\_name\>.R2.fastq for PAIRED-END data and \<sample\_name\>.fastq for SINGLE END data.

**IMPORTANT**: Please note that MToolBox cannot recognize more than one PAIRED-END couple of fastq files (R1+R2) and one SINGLE-END fastq file per sample.


###MTOOLBOX PROGRAMS OPTIONS


Besides the mandatory ```-i``` option, the execution of MToolBox can be fine-tuned by tweaking parameters of the *mapExome.py* and *assembleMTgenome.py* scripts by running MToolBox as follows:
```
MToolBox.sh -i <input_format> -r <reference_sequence> -m "<mapExome_options>" -a "<assembleMTgenome_options>" -c "<mt-classifier_options>"

```
**Most relevant options**:

```-i``` (*MToolBox.sh*) to choose the input file(s) format (mandatory).

```-r``` (*MToolBox.sh*) to choose the mitochondrial reference sequence for read mapping. Please note that the selected reference sequence will be used as reference for the VCF output file. Allowed values are RCRS and RSRS. Default is RSRS.

```-M``` (*MToolBox.sh*h) to enable duplicate read removal by MarkDuplicates.

```-I``` (*MToolBox.sh*) to enable mapped reads realignment around indels annotated in MITOMAP and HMTDB by GenomeAnalysisTK.jar.

```-t``` (*mapExome.py*) to set the number of threads used by gsnap. Default is 8.

```-t``` (*assembleMTgenome.py*) to set the minimum distance from the read end required to retain an indel for variant calling. Default is 5. Please note that only values >= 5 are allowed.

```-z``` (*assembleMTgenome.py*) to set the minimum heteroplasmy threshold for variants to be reported in the FASTA consensus sequence. Default is 0.80.

**For the full list of MToolBox options, please run as follows**:
```
MToolBox.sh -h
```
**For the full list of mapExome and assembleMTgenome options, please run as follows**:
```
mapExome.py -h
assembleMTgenome.py -h
mt-classifier.py -h
```
**Usage example**:
```
MToolBox.sh -M -I -m "-t 20 -M /path/to/gsnapdb/chrRSRS/ -H /path/to/gsnapdb/hg19RSRS/" -a "-t 10"
```
this command will enable duplicate read removal (-M), enable mapped reads realignment (-I); also, the number of threads used by gsnap will be set to 20 (-m "-t 20..") and the custom path for gsnap mitochondrial and nuclear indexes will be set (-M and -H respectively). Finally, the minimum nucleotide distance from read end to retain indels will be set to 10 (-a "-t 10").


###CONTACTS AND CITATION

Contacts:
dome.simone [at] gmail.com ; 
claudia.calabrese23 [at] gmail.com

##If you use MToolBox, please cite:
Calabrese C, Simone D, Diroma MA, Santorsola M, Guttà C, Gasparre G, Picardi E, Pesole G, Attimonelli M. MToolBox: a highly automated pipeline for heteroplasmy annotation and prioritization analysis of human mitochondrial variants in high-throughput sequencing. *Bioinformatics*. 2014 Jul 14. pii: btu483. Epub ahead of print] PubMed PMID: 25028726

##Get support 
Please join the MToolBox google group:
https://groups.google.com/forum/?hl=IT#!forum/mtoolbox-users
