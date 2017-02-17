#!/bin/bash


usage()
{
	USAGE="""
Auxiliary script to get summary.txt and prioritized_variant.txt files from MToolBox annotation.csv files.
Please export the MToolBox folder with executables before running this script!!!

Usage:
get_prioritzed_var_and_summary.sh  -i config.sh -o output

		-i	config file [MANDATORY]
		-t	Heteroplasmy threshold used to generate the fasta file [default 0.8]
		-d	Minimum distance of indels from read end used for variant calling [default 5]
		-o	out dir where OUT_filename directories are located [MANDATORY]
	Help options:
		-h	show this help message
	"""
	echo "$USAGE"
}



UseMarkDuplicates=false
UseIndelRealigner=false
MitoExtraction=false
HFthreshold=0.8
REdistance=5
while getopts ":hi:t:d:o:" opt; do
	case $opt in
		h)
			usage
			exit 1
			;;
		i)
			config=$OPTARG
			;;
		t)
			HFthreshold=$OPTARG
			;;
		d)
			REdistance=$OPTARG
			;;
		o)	
			outdir=$OPTARG
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

if [ -f "$config" ]
then
	echo -e '\nsetting up MToolBox environment variables...'
	source $config
	echo -e 'done\n'
else
	usage
	exit 1
fi

cd $outdir


echo "Output files in $outdir "

for i in $(ls -d OUT_*/*annotation.csv); do tail -n+2 $i | awk 'BEGIN {FS="\t"}; {if ($5 == "yes" && $6 == "yes" && $7 == "yes") {print $1"\t"$2"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$30"\t"$31"\t"$32"\t"$33"\t"$34"\t"$35"\t"$36"\t"$37"\t"$38"\t"$39"\t"$40"\t"$41"\t"$42"\t"$43"\t"$44}}' >> priority_tmp.txt; done

for i in $(ls -d OUT_*/*annotation.csv); do tail -n+2 $i | awk 'BEGIN {FS="\t"}; {if ($5 == "yes" && $6 == "yes" && $7 == "yes") count++} END {print $1"\t"NR"\t"count}' >> variant_number.txt; done
prioritization.py priority_tmp.txt

rm priority_tmp.txt

echo ""
echo "Prioritization analysis done."
echo ""

for i in $(ls -d OUT_*); do name=$(echo $i | sed 's/OUT_//g'); cd $i; coverage=$(cat *coverage.txt | grep "Assemble"); cd ..; echo "Sample:" "$name" "$coverage"; done >> coverage_tmp.txt

for i in $(ls -d OUT_*); do name=$(echo $i | sed 's/OUT_//g'); cd $i; heteroplasmy=$(echo "$HFthreshold"); homo_variants=$(awk 'BEGIN {FS="\t"}; {if ($3 == "1.0") count++} END {print count}' *annotation.csv); above_threshold=$(awk -v thrsld=$heteroplasmy 'BEGIN {FS="\t"};{if ( $3 >= thrsld && $3 < "1.0" ) count++} END {print count}' *annotation.csv); under_threshold=$(awk -v thrsld=$heteroplasmy 'BEGIN {FS="\t"};{if ( $3 < thrsld && $3 > "0" ) count++} END {print count}' *annotation.csv); cd ..; echo "$name" "$homo_variants" "$above_threshold" "$under_threshold"; done >> heteroplasmy_count.txt

summary.py coverage_tmp.txt heteroplasmy_count.txt
rm coverage_tmp.txt
rm heteroplasmy_count.txt
rm variant_number.txt

echo -e "Selected input format\t$(echo "$input_type")\nReference sequence chosen for mtDNA read mapping\t$(echo "$ref")\nReference sequence used for haplogroup prediction\tRSRS\nDuplicate read removal?\t$(echo "$UseMarkDuplicates")\nLocal realignment around known indels?\t$(echo "$UseIndelRealigner")\nMinimum distance of indels from read end\t$(echo "$REdistance")\nHeteroplasmy threshold for FASTA consensus sequence\t$(echo "$HFthreshold")\n\nWARNING: If minimum distance of indels from read end set < 5, it has been replaced with 5\n\n==============================\n\n$(cat summary_tmp.txt | sed "s/thrsld/$HFthreshold/g")\n\n==============================\n\nTotal number of prioritized variants\t$(awk 'END{print NR-1}' prioritized_variants.txt)" > summary_`date +%Y%m%d_%H%M%S`.txt

rm summary_tmp.txt
echo ""
echo "Analysis completed!"
cd ..

