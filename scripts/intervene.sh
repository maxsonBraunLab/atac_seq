# intersect peaks for each sample and condition.

# we are looking for >= 1 intersections because that region would be found in 2 (out of 3) replicates.

bedtools intersect -c -a data/macs2/AFC_1/AFC_1_peaks.broadPeak -b data/macs2/AFC_2/AFC_2_peaks.broadPeak data/macs2/AFC_3/AFC_3_peaks.broadPeak | awk '{if ($10 >= 1) {print}}' | cut -f1-3 > data/intervene/AMLETO_FIRST.bed

bedtools intersect -c -a data/macs2/ACF_1/ACF_1_peaks.broadPeak -b data/macs2/ACF_2/ACF_2_peaks.broadPeak data/macs2/ACF_3/ACF_3_peaks.broadPeak | awk '{if ($10 >= 1) {print}}' | cut -f1-3 > data/intervene/CSF3R_FIRST.bed

bedtools intersect -c -a data/macs2/EC_1/EC_1_peaks.broadPeak -b data/macs2/EC_2/EC_2_peaks.broadPeak data/macs2/EC_3/EC_3_peaks.broadPeak | awk '{if ($10 >= 1) {print}}' | cut -f1-3 > data/intervene/CSF3R_ONLY.bed

bedtools intersect -c -a data/macs2/EE_1/EE_1_peaks.broadPeak -b data/macs2/EE_2/EE_2_peaks.broadPeak data/macs2/EE_3/EE_3_peaks.broadPeak | awk '{if ($10 >= 1) {print}}' | cut -f1-3 > data/intervene/EMPTY.bed

# plot intersections with intervene.

intervene venn -i data/intervene/AMLETO_FIRST.bed data/intervene/CSF3R_FIRST.bed data/intervene/CSF3R_ONLY.bed data/intervene/EMPTY.bed --save-overlaps --project amleto_intersect --figtype png --colors='#1abab6','#ffd47d','#f1446c','#a2aaaf' -o data/intervene/intervene_results