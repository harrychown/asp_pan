conda activate gffread_env
cd ~/scratch/imperial/Afum/C30/

# Check start codons
grep -P "\tstart_codon\t" updated_masked_genemark.gtf > start_codon.gtf
sed 's/start_codon/gene/g' start_codon.gtf > fake_start.gtf
sed 's/start_codon/CDS/g' start_codon.gtf >> fake_start.gtf
gffread -w codons_C30_start.fa -g updated_masked_C30.contigs.fa fake_codon.gtf
# Remove FASTA header line
grep -v ">" codons_C30_start.fa | grep -v "ATG" > no_true_start.fa

# Check stop codons
grep -P "\tstop_codon\t" updated_masked_genemark.gtf > stop_codon.gtf
sed 's/stop_codon/gene/g' stop_codon.gtf > fake_stop.gtf
sed 's/stop_codon/CDS/g' stop_codon.gtf >> fake_stop.gtf
gffread -w codons_C30_stop.fa -g updated_masked_C30.contigs.fa fake_stop.gtf
# Remove FASTA header line
grep -v ">" codons_C30_stop.fa | grep -v "TGA" | grep -v "TAA" | grep -v "TAG" > no_true_stop.fa
