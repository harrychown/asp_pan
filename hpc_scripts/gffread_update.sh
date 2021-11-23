conda activate gffread_env
gffread -J -y ./recall_ep/true_C30.pep.fa  -w ./recall_ep/true_C30.masked.transcripts.fa -g ./updated_masked_C30.contigs.fa ./recall_ep/updated_masked_genemark.gtf
