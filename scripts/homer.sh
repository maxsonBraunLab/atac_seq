#!/usr/bin/bash

while getopts "i:g:" op
do
	case "$op" in
		i)  i="$OPTARG";;
		g)  g="$OPTARG";;
		\?) exit 1;;
	esac
done

# reset current homer results
if [ -d "data/homer" ]; then
	rm -r "data/homer"
fi

# check if input files exists.
if [ -f $i ];
then
	echo "contrast combinations exists"
else
	echo -e "$i does not exist. Exiting program."
	exit
fi

if [ -f $g ];
then
	echo "genome file exists"
else
	echo -e "$g does not exist. Exiting program."
	exit
fi

log_file_exists() {

	if [ -f $1 ]; then
		rm $1
		touch $1
	fi

	if [ ! -f $1 ]; then
		touch $1
	fi

}

# import contrast_combinations.txt
while IFS=$'\t' read condition_1 condition_2
do

	# define contrast
	contrast=$(echo "$condition_1-vs-$condition_2")

	# define up and downregulated input intervals and log files
	up_peaks=$(echo "data/deseq2/$contrast/$contrast-up_sig.bed")
	dn_peaks=$(echo "data/deseq2/$contrast/$contrast-down_sig.bed")
	up_log=$(echo "data/logs/homer-$contrast.log")
	dn_log=$(echo "data/logs/homer-$contrast.log")

	# check and make empty log files
	log_file_exists "$up_log"
	log_file_exists "$dn_log"

	# check input bed files before running HOMER on up and down peaks
	if [ -f "$up_peaks" ]; then
		echo "Upregulated peaks file detected."
	fi

	if [ -f "$dn_peaks" ]; then
		echo "Downregulated peaks file detected."
	fi

	up_peaks_count=$(wc -l "$up_peaks" | cut -d" " -f1)
	dn_peaks_count=$(wc -l "$dn_peaks" | cut -d" " -f1)

	if [ "$up_peaks_count" -lt 10 ]; then
		echo "ERROR: $up_peaks had less than 10 DE intervals."
	else
		echo "Running HOMER for $up_peaks_count $up_peaks"
		findMotifsGenome.pl "$up_peaks" "$g" data/homer/"$contrast"-up -size 200 > "$up_log" &
	fi

	if [ "$dn_peaks_count" -lt 10 ]; then
		echo "ERROR: $dn_peaks had less than 10 DE intervals."
	else
		echo "Running HOMER for $dn_peaks_count $dn_peaks"
		findMotifsGenome.pl "$dn_peaks" "$g" data/homer/"$contrast"-down -size 200 > "$dn_log" 2>&1 &
	fi

done < $i

wait