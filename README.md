# Generating a pan-genome for _Aspergillus fumigatus_

Some scripts found within this repository were used on HPC machines, therefore these may not run on local machines. The majority of software was downloaded/installed using Conda and a Conda environment was generated for each piece.

**Sequence data from 218 isolates was used to generate _A. fumigatus_ genomes**

Paired-end Illumina sequencing data was obtained alongside isolate metadata. (Links to these resources will be published once made publicly available). 
File "assembly_pipeline.job" assembles the genomes for the isolates in parallel. 

1. Raw data first passes through quality control (Trimmomatic)
2. Genomes are then assembled (Megahit)
3. Low-complexity regions are masked (RepeatMasker)
4. Assembly statistcs produced (BBtools)
5. Protein-coding regions predicted (ProtHint, GeneMark-EP+)
6. Descriptive files generated for PanOCT input

**Protein-protein homology was determined for all isolates**

Parallel BLASTP searches were performed by first building a database of all proteins from the study (script: "build_blast.job") and then aligning them to the database (script: "parallel_blast.job"). During the script for database construction, PanOCT inputs were produced.

**Pan-genome generation and analyses**

The pan-genome was constructed using PanOCT (script: panoct.job). Subsequent analysis was performed using "panoct_data.R" and association testing carried out using Scoary (script: "scoary.sh" followed by "scoaryout.r", naming dependent on DIAMOND).

**Annotation of pan-genome was performed using DIAMOND**

The scripts to annotate the pan-genome are also included here in "diamond.job"

**Additional scripts**
Two additional scripts are also included:
1. Grouping gene regions for future MSA analysis ("groupcentroids.job")
2. Generation of assembly statistics table ("generate_stats.job")



