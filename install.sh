#!/bin/bash

usage()
{
	USAGE="""
MToolBox: a tool for heteroplasmy annotation and accurate functional analysis of mitochondrial variants from high throughput sequencing data.
Written by Domenico Simone, Claudia Calabrese and Maria Angela Diroma 2013-2014.
	
Usage:
./install.sh -g <gsnap_version> -z <zlib_version> -m <muscle_file> -s <samtools_version> -k <kmer_to_build_gsnap_db>

	-g	default is 2015-12-31.v7
	-z	default is 1.2.8
	-m 	default is muscle3.8.31_i86linux64
	-s 	default is 1.3
	-k	default is 15

Help options:

	-h	show this help message
	"""
	echo "$USAGE"
}

gsnap_gmap_version=2015-12-31.v7
muscle_file=muscle3.8.31_i86linux64
samtools_version=1.3
zlib_version=1.2.8
kmer=15

while getopts ":hg:s:m:z:k:" opt; do
	case $opt in
		h)
			usage
			exit 1
			;;
		g)
			gsnap_gmap_version=$OPTARG
			;;
		s)
			samtools_version=$OPTARG
			;;
		m)
			muscle_file=$OPTARG
			;;
		z)
			zlib_version=$OPTARG
			;;
		k)
			kmer=$OPTARG
			;;
		\?)
			echo "Invalid option: -$OPTARG" >&2
			exit 1
			;;
		:)
			echo "Option -$OPTARG requires an argument." >&2
			exit 1
			;;
	esac
done

download()
{
    URL=$1
    OFILE=$2
    wget  --no-check-certificate -c -nc -L "$URL" -O $OFILE
}

muscle_install()
{
	echo "Installing muscle v $muscle_file ..."
	muscle_file=$muscle_file
	muscle_url=http://www.drive5.com/muscle/downloads3.8.31/${muscle_file}.tar.gz
	pushd .
	cd bin
	download $muscle_url $muscle_file.tar.gz
	tar -zxvf $muscle_file
	popd
	echo "Installing muscle $muscle_file ...done."
}

samtools_install()
{
	echo "Installing samtools $samtools_version"
	samtools_version=$samtools_version
	samtools_file=samtools-${samtools_version}.tar.bz2
	samtools_url=https://github.com/samtools/samtools/releases/download/1.3/$samtools_file
	pushd .
	cd bin
	download $samtools_url $samtools_file
	tar xjvf $samtools_file
	cd samtools-${samtools_version}
	./configure 
	make
	make prefix=../samtools/ install
	popd
	echo "Installing samtools $samtools_file ...done."
}

zlib_install() {
	echo "Installing zlib $zlib_version ..."
	zlib_version=$zlib_version
	zlib_file=zlib-${zlib_version}.tar.gz
	zlib_url=http://zlib.net/$zlib_file
	pushd .
	cd bin
	download $zlib_url $zlib_file
	tar xvzf $zlib_file
	cd zlib-${zlib_version} 
	./configure --prefix ../zlib
	make 
	make install
	popd
	echo  "Installing zlib $zlib_version ...done."
}

gsnap_gmap_install() {
	echo "Installing gmap $gsnap_gmap_version ..."
	gsnap_gmap_file=gmap-gsnap-${gsnap_gmap_version}.tar.gz
	gsnap_gmap_url=http://research-pub.gene.com/gmap/src/$gsnap_gmap_file
	pushd .
	cd bin
	download $gsnap_gmap_url $gsnap_gmap_file
	tar xvzf $gsnap_gmap_file &&
	rm $gsnap_gmap_file
	cd gmap*
	./configure --prefix=$MTOOLBOX_BIN/gmap
	make
	make install
	echo "Installing gmap $gsnap_gmap_version...done"
	popd
}

gmap_db() {
	nfasta_rcrs_url=http://sourceforge.net/projects/mtoolbox/files/genome_fasta/hg19RCRS.fa.gz
	nfasta_rcrs=hg19RCRS.fa.gz
	nfasta_rsrs_url=http://sourceforge.net/projects/mtoolbox/files/genome_fasta/hg19RSRS.fa.gz
	nfasta_rsrs=hg19RSRS.fa.gz
	rcrs_mfasta_url=http://sourceforge.net/projects/mtoolbox/files/genome_fasta/chrM.fa.gz
	rcrs_mfasta=chrM.fa.gz
	rsrs_mfasta_url=http://sourceforge.net/projects/mtoolbox/files/genome_fasta/chrRSRS.fa.gz
	rsrs_mfasta=chrRSRS.fa.gz
	download $nfasta_rcrs_url $nfasta_rcrs
	download $nfasta_rsrs_url $nfasta_rsrs
	download $rcrs_mfasta_url $rcrs_mfasta
	download $rsrs_mfasta_url $rsrs_mfasta
	if [ "$decompress_fasta"==True ]; 
	then
		echo "Building gmap db using $nfasta_rcrs..."
		$MTOOLBOX_BIN/gmap/bin/gmap_build -D $gmap_db -d $database_name_nfasta_rcrs -g $nfasta_rcrs -s numeric-alpha -k $kmer 
		echo "Building gmap db using $nfasta_rsrs..."
		$MTOOLBOX_BIN/gmap/bin/gmap_build -D $gmap_db -d $database_name_nfasta_rsrs -g $nfasta_rsrs -s numeric-alpha -k $kmer
		echo "Building gmap db using $rcrs_mfasta"
		$MTOOLBOX_BIN/gmap/bin/gmap_build -D $gmap_db -d $database_name_rcrs -g $rcrs_mfasta -s numeric-alpha -k $kmer
		echo "Building gmap db using $rsrs_mfasta"
		$MTOOLBOX_BIN/gmap/bin/gmap_build -D $gmap_db -d $database_name_rsrs -g $rsrs_mfasta -s numeric-alpha -k $kmer
	else
                echo "Building gmap db using $nfasta_rcrs..."
                $MTOOLBOX_BIN/gmap/bin/gmap_build -D $gmap_db -d $database_name_nfasta_rcrs $nfasta_rcrs -s numeric-alpha -k $kmer
                echo "Building gmap db using $nfasta_rsrs..."
                $MTOOLBOX_BIN/gmap/bin/gmap_build -D $gmap_db -d $database_name_nfasta_rsrs $nfasta_rsrs -s numeric-alpha -k $kmer
                echo "Building gmap db using $rcrs_mfasta"
                $MTOOLBOX_BIN/gmap/bin/gmap_build -D $gmap_db -d $database_name_rcrs $rcrs_mfasta -s numeric-alpha -k $kmer
                echo "Building gmap db using $rsrs_mfasta"
                $MTOOLBOX_BIN/gmap/bin/gmap_build -D $gmap_db -d $database_name_rsrs $rsrs_mfasta -s numeric-alpha -k $kmer
	fi
	echo "Building gmap db ...done"
}
#########################################
decompress_fasta=True
database_name_rcrs=chrM
database_name_rsrs=chrRSRS
database_name_nfasta_rcrs=hg19RCRS
database_name_nfasta_rsrs=hg19RSRS
DIRECTORY=bin
if [ -d "$DIRECTORY" ]; then
	echo "$DIRECTORY already exists..."
else
	echo "creating $DIRECTORY..."
	mkdir $DIRECTORY
fi

MTOOLBOX_BIN=$(readlink -f $DIRECTORY)

GENOME_FASTA=genome_fasta
if [ -d "$GENOME_FASTA" ]; then
        echo "$GENOME_FASTA already exists..."
else
        echo "creating $GENOME_FASTA..."
        mkdir $GENOME_FASTA
fi

gmap_db=gmapdb
if [ -d "$gmap_db" ]; then
        echo "$gmap_db already exists..."
else
        echo "creating $gmap_db..."
        mkdir $gmap_db
fi


echo "Installing zlib, samtools, muscle, gsnap/gmap...wait"
zlib_install
samtools_install
muscle_install
gsnap_gmap_install
echo "Installing zlib, samtools, muscle, gsnap/gmap...done"
echo "exporting zlib to generate gmap db..."
export LD_LIBRARY_PATH=$MTOOLBOX_BIN/zlib/lib:$LD_LIBRARY_PATH:/usr/local/lib
export CFLAGS="-I$MTOOLBOX_BIN/zlib/include $CFLAGS"
echo "Building gmap database..."
gmap_db
echo "Building gmap database ...done"
echo "moving fasta files into $GENOME_FASTA"

mv $rcrs_mfasta $GENOME_FASTA
mv $rsrs_mfasta $GENOME_FASTA
mv $nfasta_rcrs $GENOME_FASTA
mv $nfasta_rsrs $GENOME_FASTA

if [ "$decompress_fasta"==True ]
then
	cd $GENOME_FASTA
	echo "decompress nfasta..."
	for i in $(ls *.gz); do gzip -d $i; done &&
	cd ..
	echo "Done"
else
	echo "Done"
fi

echo "Generating fasta.fai indexes with samtools..."
cd $GENOME_FASTA
ls * > fasta.lst
for i in $(cat fasta.lst); do $MTOOLBOX_BIN/samtools/bin/samtools faidx $i; done
cd ..

MTOOLBOX_DIR=$(pwd)
SETUP_FILE=${MTOOLBOX_DIR}/MToolBox/setup.sh

cat <<EOF > $SETUP_FILE
#! /bin/bash

export MTOOLBOX_DIR=$MTOOLBOX_DIR
export MTOOLBOX_BIN=$MTOOLBOX_BIN
export LD_LIBRARY_PATH=\$MTOOLBOX_BIN/zlib/lib:\$LD_LIBRARY_PATH:/usr/local/lib
export CFLAGS="-I\$MTOOLBOX_BIN/zlib/include \$CFLAGS"
export samtoolsexe=\$MTOOLBOX_BIN/samtools-${samtools_version}/samtools
export muscleexe=\$MTOOLBOX_BIN/$muscle_file
export gsnapexe=\$MTOOLBOX_BIN/gmap/bin/gsnap
fasta_path=\$MTOOLBOX_DIR/$GENOME_FASTA
gsnapdb=\$MTOOLBOX_DIR/$gmap_db
mtdb_fasta=chrRSRS.fa
hg19_fasta=hg19RSRS.fa
mtdb=chrRSRS
humandb=hg19RSRS
EOF

exit 0
