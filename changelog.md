###CHANGELOG - 12 March 2016

- a bug in the inclusions of *insertions*/*deletions* into the consensus fasta sequence has been fixed in the `assembleMTgenome.py -z` option.

WARNING: We discourage MToolBox users to set the HF threshold (using the `assembleMTgenome.py -z` parameter) to values <0.5 since this might evidence some ambiguities between low heteroplasmy insertions / deletions and mismatches that cannot be handled by IUPAC ambiguity encoding. In this case, we decided to *arbitrarly* report in the fasta consensus sequence:
-the *insertion* event, with an ambiguity of three mutations types (insertion/deletion/mismatch )
-the *deletion* event, with an ambiguity of two mutations types (deletion and mismatch)

We would also like to remind the MToolBox users that the consensus fasta file (and related `annotation.csv` file) might not report the full set of mitochondrial insertions and deletions found in the sample under study if they follow below the chosen HF, due to the limit imposed by the IUPAC ambiguity encoding. However, the full set of variations found is still available in the `VCF_file.vcf`, disregarding of their HF.
 

###CHANGELOG - 6 March 2016

- an error in the samtools path specification has been fixed in the MToolBox.sh script. From now on the users MUST specify the full path to samtools, if they are not installed in /usr/local/bin/ default directory, using the MToolBox.sh `-s` option. We apologize with the MToolBox users for this inconvenience.


###CHANGELOG - 1 March 2016

Update to MToolBox version 0.3.3 with the following change:

- GSNAP databases available at https://sourceforge.net/projects/mtoolbox/files/genome_index/ generated with the GSNAP version 2013-09-11 have been removed, due to an inconsistency between the mitochodrial reference rCRS sequence fasta and index database name. New GSNAP nuclear-mitochondrial and mitochondrial databases have been uploaded, generated with the GSNAP version 2015-12-31.v7. Please, be aware that these databases might not be compatible with previous GSNAP versions. 

-  Mitochondrial rCRS reference fasta and GSNAP database default names used in the MToolBox.sh script were changed as following:

`chrRCRS.fa` --> `chrM.fa`

`chrRCRS` (GSNAP db) --> `chrM`

to reflect changes in GSNAP databases uploaded at https://sourceforge.net/projects/mtoolbox/files/genome_index/. 


###CHANGELOG - 21 September 2015

Update to MToolBox version 0.3.2 with the following change:

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

###CHANGELOG - 8 June 2015

Update to MToolBox version 0.3.1 with the following change:

- A bug in the MToolBox.sh script has been fixed. Empty paired-end fastq files generated during SAM/BAM to fastq conversion are now removed.

###CHANGELOG - 28 February 2015

Update to MToolBox version 0.3 with the following new options and changes:

 - fastq.gz is a further possible input format file. Installation of zlib libraries is therefore required. 
 - users can specify the path of the working directories using ```–p``` (path to input folder) and ```–o``` (path to output folder) options.
 - users can specify a list of files to be used as input through -l option. It accepts a text file containing one sample name for each line. This list should be named as "list.txt" and placed in the input folder. Alternatively, users can provide comma-separated names with the same option. It is mandatory to report in such list the filename extension (e.g. mysample.sam or mysample.R1.fastq).
 - users can use -X option to allow the extraction from a BAM file of mitochondrial reads mapped onto a mitochondrial reference sequence. This option can be useful when using Whole Genome or Exome sequencing BAM files containing a huge amount ofnuclear reads. The option works only with the BAM format.

new fields added to the annotation.csv output file:
 - Disease Score: "% Disease", "% Neutral" and "% Unclassified" fields have been replaced with an overall Disease Score, generated as a weighted average of pathogenicity prediction scores for non-synonymous variants, derived from a training dataset of 53 non synonymous variants selected among mitochondrial diseases or cancer associated mutations. Weights have been calculated by taking into account the right prediction and the best probability to predict a truly pathogenic mtDNA variants generated by the pathogenicity prediction algorithms currently implemented in MToolBox.
 - MutPred Pred: MutPred prediction (Low pathogenicity, High pathogenicity).
 - dbSNP ID: Variant ID in dbSNP.

for users convenience, new scripts to help the generation of reports about the annotated and prioritized variants have been added to the suite of tools provided by MToolBox:
 - prioritization.py, which generates the prioritized_variants.txt file, reporting annotation only for prioritized variants for each sample analyzed, defined as variants recognized by the three reference sequences (rCRS, RSRS and MHCS), sorted per increasing variability.
 - summary.py, which generates the summary.txt file, reporting statistics about the coverage of reconstructed mitochondrial genomes, number of homoplasmic and heteroplasmic variants (for NGS data), haplogroup prediction and number of prioritized variants.

###CHANGELOG - 25 February 2015

 - An error occurred during the generation of circularized mitochondrial chromosome, used for the gsnap db generation. Hg19RCRS/hg19RSRS and chrRCRS/chrRSRS gsnap indexed databases have been replaced with those using the linearized mitochondrial chromosome. We apologize with the MToolBox users for this inconvenient. 

###CHANGELOG - 23 January 2015

Update to MToolBox version 0.2.2:

 - an error encountered during the analysis of hard clipping mapped reads has been fixed in the mtVariantCaller.py
 - the mtVariantCaller.py has been improved to better manage sites with multiple alleles.
 - hidden files included in the package and generating problems with the mt-classifier.py have been eliminated.

###CHANGELOG - 6 November 2014

Update to MToolBox version 0.2.1:

 - an error encountered during the sam to fastq extraction has been fixed in the MToolBox.sh file.
 - RCRS, hg19RCRS, RSRS and hg19RSRS gmap indexed databases have been regenerated using the -c option for circularized chromosomes.
 - update to Phylotree Build 16 for haplogroup prediction. 

###CHANGELOG - 29 September 2014

- an error fixed in the -t parameter of assembleMTgenome.

###CHANGELOG - 19 July 2014

- an error encountered during the bam to fastq extraction has been fixed. Empty unpaired fastq files are now removed.  

###CHANGELOG - 5 June 2014

- changed the method for estimation of the heteroplasmy confidence interval (CI). For sites with coverage depth <= 40, the heteroplasmy CI is estimated with the Wilson score interval; for larger coverage depth values, the Agresti-Coull interval is used.
- added the possibility to use fasta inputs to perform haplogroup prediction and functional annotation.
- added the possibility to use the revised Cambridge Reference Sequence (rCRS) as reference sequence for read mapping. By using rCRS as reference sequence, the VCF output will be rCRS-based. 
