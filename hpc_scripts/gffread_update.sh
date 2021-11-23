#!/bin/bash --login
#$ -cwd
#$ -t 1-218
#$ -N build_genome2
#$ -pe smp.pe 12
module load apps/binapps/anaconda3/2020.07
module load tools/env/proxy2 

conda activate imperial

FOLDER=`awk "NR==$SGE_TASK_ID" folder_paths.txt`

# Enter the isolate folder
cd $FOLDER
# Extract prefix
PREFIX=$(echo $FOLDER | rev | cut -d"/" -f1 | rev)


conda activate gffread_env
gffread -J -y true_$PREFIX.pep.fa  -w true_$PREFIX.masked.transcripts.fa -g updated_masked_$PREFIX..contigs.fa updated_masked_genemark.gtf
original_size=$(grep -c ">" $PREFIX.masked.transcripts.fa)
new_size=$(grep -c ">" true_$PREFIX.masked.transcripts.fa)
if [ -e ~/scratch/imperial/gene_number.log ]
then
	if grep -q "${PREFIX}," ~/scratch/imperial/gene_number.log; then
	    rm ~/scratch/imperial/gene_number.log
	fi
fi
echo "${PREFIX},${original_size},${new_size}" >> ~/scratch/imperial/gene_number.log

