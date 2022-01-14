source activate blastenv
cd /data/pan_transposon

GENOME=ARAF001.contigs.fa
: <<'END'
# 1) PERFORM ALIGNMENTS
makeblastdb -in $GENOME -dbtype nucl -out genome_db
# Search the ORF gainst genome
blastn -db genome_db -query orf.fasta -out orf_vs_genome -outfmt "6 qseqid sseqid frames pident qcovhsp length qstart qend sstart send evalue"  
# Search TIR against genome
blastn -db genome_db -query tir.fasta -out tir_vs_genome -outfmt "6 qseqid sseqid frames pident qcovhsp length qstart qend sstart send evalue" 
END




# We could use bedtools closest to identify if they are flanking
# To do this, we need to first make a BED file out of the BLAST results
cat balst_tabular_output.txt | awk '{print($2"\t"$7-1"\t"$8)}' > blast_otput.bed



