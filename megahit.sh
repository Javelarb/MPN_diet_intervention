#!/bin/bash
#$ -N MPN_assembly2
#$ -q bigmemory
#$ -pe openmp 80
#$ -ckpt restart
#$ -m beas

module purge
module load megahit/1.1.1

#Don't add '/' at end.
MEGAHIT=/dfs5/bio/javelarb/MPN_fiber_SEH/all_MPN_stuff/MPN_data/megahit2
QUALITY=/dfs5/bio/javelarb/MPN_fiber_SEH/all_MPN_stuff/MPN_data/fastqs

readlink -f $QUALITY/*.1.fastq.gz > READ_1.txt #Places reads file names in file
readlink -f $QUALITY/*.2.fastq.gz > READ_2.txt
sed -i 's/$/,/' READ_1.txt #adds commas to end of each file name
sed -i 's/$/,/' READ_2.txt
sed -i '$ s/.$//' READ_1.txt #deletes last comma of last line
sed -i '$ s/.$//' READ_2.txt
READ1=`cat READ_1.txt | tr -d '\n'` #sets variable to echo file contents without newline
READ2=`cat READ_2.txt | tr -d '\n'`

megahit \
-1 ${READ1} \
-2 ${READ2} \
-t 80 \
--min-count 3 \
--k-list 31,41,51,61,71,81,91,101,111 \
--kmin-1pass \
--min-contig-len 2500 \
--memory 0.90 \
--out-dir ${MEGAHIT} \
--continue

#megahit --continue --out-dir ${MEGAHIT} -t 80 --memory 0.85
