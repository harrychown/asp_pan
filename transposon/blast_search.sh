source activate blastenv
cd /data/pan_transposon
#makeblastdb -in centroids.fasta -dbtype prot -out pan_blast_db
# Search the protein of ORF against pan-genome
#blastp -db pan_blast_db -query protein.fasta -out protein_vs_pangenome -outfmt "6 qseqid sseqid pident qcovhsp length qstart qend sstart send evalue" 


#makeblastdb -in ARAF001.pep.fa -dbtype prot -out araf_pep_db
# Search the protein of ORF against ARAF001 proteome
#blastp -db araf_pep_db -query protein.fasta -out protein_vs_proteome -outfmt "6 qseqid sseqid pident qcovhsp length qstart qend sstart send evalue" 

makeblastdb -in ARAF001.contigs.fa -dbtype nucl -out araf_db
# Search the ORF gainst ARAF001 genome
blastn -db araf_db -query orf.fasta -out orf_vs_araf -outfmt "6 qseqid sseqid frames pident qcovhsp length qstart qend sstart send evalue"  
# Search the locus against ARAF001 genome
#blastn -db araf_db -query locus.fasta -out locus_vs_araf -outfmt "6 qseqid sseqid pident qcovhsp length qstart qend sstart send evalue" 



# Search TIR against ARAF001 genome
blastn -db araf_db -query tir.fasta -out tir_vs_araf -outfmt "6 qseqid sseqid frames pident qcovhsp length qstart qend sstart send evalue" 


