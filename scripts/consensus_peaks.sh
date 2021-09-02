#!/bin/bash

# import + check args
n_intersects=$1
all_inputs=$(echo "${@:2}")

if ! [[ "$n_intersects" =~ ^[0-9]+$ ]]; then
	echo "First argument must be an integer. Your input: $n_intersects"
	exit 1
fi

all_conditions=$(echo $all_inputs | tr ' ' '\n' | cut -d/ -f3 | cut -d_ -f1 | sort | uniq)

for condition in $all_conditions; do

	# file I/O
	tmp_output="data/macs2/$condition.tmp.bed"

	# list all replicates in one condition
	all_replicates=$(find data/macs2/ -name "*$condition*.broadPeak" | sort | tr '\n' ' ')

	# find widest peak that appear in at least n replicates + export to tmp file.
	cat $all_replicates | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge | \
	bedtools intersect -c -a stdin -b $all_replicates | \
	awk -v n=$n_intersects '$4 >= n' | awk -v OFS='\t' '{print $1,$2,$3}' > $tmp_output

done

# merge intervals for all tmp files and export as consensus peak.
all_temp_files=$(find data/macs2 -name "*.tmp.bed" | sort | tr '\n' ' ')
cat $all_temp_files | sort -k1,1 -k2,2n | bedtools merge
rm $all_temp_files