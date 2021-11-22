# asp_pan
Generating large pangenome for Aspergillus fumigatus

**Obtaining assembly statistics**

Assembly statistics were obtained using BBTOOLS (version:):


A comparison of the number of contigs greater than 50KB was then performed to understand which of the assemblies were the most fragmented: 

```
for file in *.fa.stats
do
stat=$(grep "Number of scaffolds > 50 KB:" $file | cut -d$'\t' -f2)
prefix=$(echo $file | cut -d"_" -f3 | cut -d"." -f1)
echo "${stat},${prefix}" >> number_50kb.stats
done
```



**Annotating the pangenome**

The pangenome produced by PanOCT lacked true gene/protein names. I attempted to obtain names for these proteins by downloading all fungal proteins found on fungal OrthoDB (date accessed: 16/11/21, version: 9.1, url: https://www.orthodb.org/v9.1/index.html?page=filelist) alongside metadata tables which linked each accession to protein name and function.

FASTA files were then appended together to produce a master version ("odbfungi.fa"):

`cat *.fs > odbfungi.fa`

Furthermore the metadata was renamed to avoid confusion:
```
mv data metadata.tab
```


A BLAST (version database from the OrthoDB fungal proteins was produced and used to query the pangenome to find homologous proteins:
```
makeblastdb -in odbfungi.fa -dbtype prot -out odb_blast_db
blastp -db odb_blast_db -query centroids.fasta -outfmt 6 -max_hsps 1 -max_target_seqs 1 -out odb_search.tsv -num_threads 16
```

Results from this was then used to reannotate the pangenome FASTA and produce a centroid:accession key, enabling human-readable downstream analysis.
```
# Create a copy of the centroids FASTA
cp centroids.fasta updated_centroids.fasta
# Clear key file
rm centroid.key
while read line
do
# Extract the centroid and homologous gene accession
centroid=$(echo $line | cut -d' ' -f1)
gene_acc=$(echo $line | cut -d' ' -f2)
# Use this to search up appropriate metadata
key_line=$(grep $gene_acc metadata.tab)
# Create a new FASTA header
protein_id=$(echo $key_line | cut -d' ' -f4)
protein_name=$(echo $key_line | cut -d' ' -f7- | sed 's/ /_/g' | sed 's/\//_/g' | sed 's/\\/_/g')
new_fasta_header="${protein_id}_${protein_name}"
# Update the FASTA file and make a key
sed -i "s/$centroid /$new_fasta_header/g" updated_centroids.fasta
echo "${centroid},${new_fasta_header} " >> centroid.key
done <odb_search.tsv
```
