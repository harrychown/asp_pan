#!/bin/bash
source activate blastenv
cd /data/pan_transposon

GENOME=ARAF001.contigs.fa

: <<'END'
# 1) PERFORM ALIGNMENTS
makeblastdb -in $GENOME -dbtype nucl -out genome_db
# Search the ORF gainst genome
blastn -db genome_db -query orf.fasta -out orf_vs_genome -outfmt 6
# Search TIR against genome
blastn -db genome_db -query tir.fasta -out tir_vs_genome -outfmt 6

END

# 2) IDENTIFY ITR-ORF-ITR's
# Generate BED files from BLAST output
rm *.bed
grep -v '^#' orf_vs_genome| perl -ane 'if($F[8]<=$F[9]){print join("\t",$F[1],$F[8]-1,$F[9],$F[0],"0","+"),"\n";}else{print join("\t",$F[1],$F[9]-1,$F[8],$F[0],"0","-"),"\n";}' | sort >> orf_vs_genome.bed
grep -v '^#' tir_vs_genome| perl -ane 'if($F[8]<=$F[9]){print join("\t",$F[1],$F[8]-1,$F[9],$F[0],"0","+"),"\n";}else{print join("\t",$F[1],$F[9]-1,$F[8],$F[0],"0","-"),"\n";}' | sort >> tir_vs_genome.bed

# Update ORF name for downstream identification
while read bedline
do
	newname=$(echo "$bedline" | cut -d$'\t' -f1-4 | sed 's/\t/-/g')
	echo "$bedline" | sed "s/atf4_orf/$newname/g" >> updated_orf_vs_genome.bed
done <orf_vs_genome.bed

# Separate TIRs into start and end files
grep "TIR_start" tir_vs_genome.bed > start_tir.bed
grep "TIR_end" tir_vs_genome.bed > end_tir.bed


conda deactivate
# Identify neighbouring regions
bedtools closest -a updated_orf_vs_genome.bed -b start_tir.bed end_tir.bed -io -s -k 2 -mdb all > orf_near_itr.txt
# Identify if they are flanking and save
while read orfline
do
	orf=$(echo "$orfline" | cut -d$'\t' -f4)
	if [ $(grep $orf orf_near_itr.txt | wc -l) -eq 2 ] 
	then
		echo "$orfline" >> flanked_orfs.bed
	fi
done <updated_orf_vs_genome.bed

# 3) EXTRACT FLANKED ORFs SEQUENCE
bedtools getfasta -fi $GENOME -fo flanked_orfs.fasta -bed flanked_orfs.bed
prefix=$(echo $GENOME | sed 's/\.contigs\.fa//g')
sed -i "s/>/>$prefix:/g" flanked_orfs.fasta


