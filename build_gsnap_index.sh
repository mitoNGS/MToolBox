#!/bin/bash

function usage() {

    cat<<EOF    
Usage: gmap_build_index.sh -D <path_to_database> -n <database_name> -f <file.fa> -k <kmer> 

	-D Destination directory for installation  [ex. /usr/local/share/gmapdb]
	-n database name [ex. hg19RCRS]
	-f single o multifasta file [ex. hg19RCRS.fa]
	-k kmer [ex. 12]
EOF

}

while getopts ":hD:n:f:k:" opt; do
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

export PATH=$GMAP:$PATH

echo 'Building gmap indexes...'

gmap_build -D $database_path -d $database_name $fasta -s numeric-alpha -k $kmer


echo 'Done'
