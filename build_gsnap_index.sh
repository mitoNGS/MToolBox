#!/bin/bash

function usage() {

    cat<<EOF    
Usage: 
build_gsnap_index.sh -D <path_to_database> -n <database_name> -f <file.fa> -k <kmer> 

	Mandatory options:
	-D Destination directory for installation  [ex. /usr/local/share/gmapdb]
	-n database name [ex. hg19RCRS]
	-f single o multifasta file [ex. hg19RCRS.fa]
	-k kmer [ex. 12]

	Optional:
	-c assumes fasta file is gzipped [default=false]

EOF

}

decompress_gzip=false

while getopts ":hD:n:f:ck:" opt; do
	case $opt in
		h)
			usage
			exit 1
			;;
		D)
			database_path=$OPTARG
			;;
		n)
			database_name=$OPTARG
			;;
		f)
			fasta=$OPTARG
			;;
		c)
			decompress_gzip=true
			;;
		k)	
			kmer=$OPTARG
			;;			
		\?)
			echo "Invalid option: -$OPTARG" >&2
			usage
			exit 1
			;;
		:)
			echo "Option -$OPTARG requires an argument." >&2
			usage
			exit 1
			;;
	esac
done

echo 'Check if gmap is installed...'

GMAP=`which gmap`
if [ "$GMAP" == "" ]; then
	echo "gmap is not installed. Exit."
	exit 1
else
	echo 'OK!'
fi


echo 'Building gmap indexes...'

if $decompress_gzip
then
	echo "file $fasta is compressed..."
	gmap_build -D $database_path -d $database_name -g $fasta -s numeric-alpha -k $kmer
else
	gmap_build -D $database_path -d $database_name $fasta -s numeric-alpha -k $kmer
fi

echo 'Done'
