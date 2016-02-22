#!/bin/bash

######################MTOOLBOX CONFIG FILE###################################################################
##If the default installation of MToolBox was used (install.sh), the user should specify 
##only the MANDATORY parameters. Otherwsie, please remove the # before the OPTIONAL parameters 
## and assing them the correct value to set the environment to run MToolBox. 
##
##
##for more information please refert to...
##to ask for support, please write to ...
##or to claudia.calabrese23 at gmail.com/dome.simone at gmail.com
##
######################SET PATH TO MTOOLBOX EXECUTABLES########################################################
##
##OPTIONAL. If MToolBox default installation (install.sh) was used, samtoolsexe path is $MTOOLBOX_BIN/samtools-{samtools_version}/samtools. Otherwise please specify the FULL PATH to samtools executables.
##
#samtoolsexe=
##
##OPTIONAL. If MToolBox default installation (./install.sh) was used, samtools_version is 1.3
##
#samtools_version=
##
##OPTIONAL. If MToolBox default installation (./install.sh) was used. musclexe path is $MTOOLBOX_BIN/muscle3.8.31_i86linux64. Otherwise, please specify the FULL PATH to muscle executables.
##
#muscleexe=
##
#OPTIONAL. If MToolBox default installation (install.sh) was used, gsnapexe path is $MTOOLBOX_BIN/gmap/bin/gsnap. Otherwise, please specify the FULL PATH to gsnap executables.
##
#gsnapexe=
##
#####################SET FILE NAMES OF REFERENCE SEQUENCE AND DATABASE TO BE USED FOR THE MAPPING STEP#######
##
#OPTIONAL. If MToolBox default installation  (install.sh) was used, fasta_path $MTOOLBOX_DIR/genome_fasta/. Otherwise, please specify the FULL PATH to fasta and fasta.fai reference sequences to be used in the mapping step
##
#fasta_path=
##
#OPTIONAL. If MToolBox default installation  (install.sh) was used, gsnapdb is $MTOOLBOX_DIR/gmapdb/. Otherwise, please specify the FULL PATH to gsnap/gmap database executables
##
#gsnapdb=
##
##OPTIONAL. If MToolBox default installation (install.sh) was used, mtdb_fasta is MT.fa (corresponding to the rCRS reference sequence). To use the RSRS reference sequence, please set mtdb_fasta=chrRSRS.fa. Otherwise, please specify the file name of the mitochondrial reference fasta sequence you want to use.
##
#mtdb_fasta=
##
##OPTIONAL. If MToolBox default installation (install.sh) was used, hg19_fasta is hg19RCRS.fa (corresponding to the GCh37/hg19 reference multifasta including the rCRS sequence). To use the nuclear multifasta including the RSRS reference sequence, please set hg19_fasta=hg19RSRS.fa. Otherwise, please specify the file name of the nuclear reference sequence you want to use.
##
#hg19_fasta=
##
#OPTIONAL. If MToolBox dafault installation (istall.sh) was used: [MT | RSRS; DEFAULT: MT], where MT refers to the rCRS database. Otherwise, please specify the name of the mitochondrial gsnap database you want to use.
##
#mtdb=
##
#OPTIONAL. If MToolBox default installation (install.sh) was used: [hg19RCRS| hg19RSRS; DEFAULT: hg19RCRS], where hg19RCRS refers to hg19 nuclear genome + rCRS mitochondrial database. Otherwise, please specify the name of the nuclear (+mitochondrial) gsnap database you want to use.
##
#humandb=
##
######################SET PATH TO INPUT/OUTPUT and PRE/POST PROCESSING PARAMETERS############################
##
##OPTIONAL. Specify the FULL PATH of the input directory. Default is the current working directory
##
#input_path=
##
##OPTIONAL. Specify the FULL PATH of the output directory. Default is the current working directory
##
output_name=test_output
##
##OPTIONAL. Specify the FULL PATH to the list of files to be analyzed. Default is use all the files with the specified file format extension
##in the current working directory
##
#list=
##
##MANDATORY. Specify the input file format extension. [fasta | bam | sam | fastq | fastq.gz]
##
input_type=fasta
##
##MANDATORY. Specify the mitochondrial reference to be used for the mapping step with mapExome. [RCRS | RSRS; DEFAULT is RSRS} 
ref=RSRS
##
##OPTIONAL. Specify if duplicate removal by MarkDuplicates should be set on. [false | true; DEFAULT is false]
UseMarkDuplicates=false
##
##OPTIONAL. Specify if realignment around ins/dels should be set on. [false | true; DEFAULT is false]
UseIndelRealigner=false
##
##OPTIONAL: specify if to exctract only mitochondrial reads from bam file provided as input. [false | true; DEFAULT is false]
MitoExtraction=false
##
################################################################################################################
