#! /bin/bash

export MTOOLBOX_DIR=/path/to/MToolBox
export MTOOLBOX_BIN=/path/to/MToolBox/bin
export PATH=$MTOOLBOX_BIN/anaconda/bin:$PATH
export PYTHONUSERBASE=$MTOOLBOX_BIN/anaconda
export LD_LIBRARY_PATH=$MTOOLBOX_BIN/zlib/lib:$LD_LIBRARY_PATH:/usr/local/lib
export CFLAGS="-I$MTOOLBOX_BIN/zlib/include $CFLAGS"
export samtoolsexe=$MTOOLBOX_BIN/samtools-1.3/samtools
export muscleexe=$MTOOLBOX_BIN/muscle3.8.31_i86linux64
export gsnapexe=$MTOOLBOX_BIN/gmap/bin/gsnap
fasta_path=$MTOOLBOX_DIR/genome_fasta/
gsnapdb=$MTOOLBOX_DIR/gmapdb/
samtools_version=1.3
mtdb_fasta=chrRSRS.fa
hg19_fasta=hg19RSRS.fa
mtdb=chrRSRS
humandb=hg19RSRS
