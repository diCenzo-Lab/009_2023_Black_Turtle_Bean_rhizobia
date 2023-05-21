#!/usr/bin/env bash

# Help message
if [ "$1" == "-h" ]; then
  echo "
  Usage: `basename $0` nanopore.fastq.gz illumina-1.fastq.gz illumina-2.fastq.gz output threads unicycler_threads template.yaml template_submol.yaml

  Where:
    nanopore.fastq.gz: The gzipped basecalled nanopore data (uncorrected)
    illumina-1.fastq.gz: The first pair of Illumina reads from paired-end sequencing
    illumina-2.fastq.gz: The secpmd pair of Illumina reads from paired-end sequencing
    Output: The file name to use for the final assembly
    Threads: The number of threads to use
    unicycler_threads: The number of threads to use with Flye
    Template.yamp: The template.yaml file needed for PGAP annotation
    Template_submol.yaml: The template_submol.yaml file needed for PGAP annotation
    "
  exit 0
fi
if [ "$1" == "-help" ]; then
  echo "
  Usage: `basename $0` nanopore.fastq.gz illumina-1.fastq.gz illumina-2.fastq.gz output threads unicycler_threads

  Where:
    nanopore.fastq.gz: The gzipped basecalled nanopore data (uncorrected)
    illumina-1.fastq.gz: The first pair of Illumina reads from paired-end sequencing
    illumina-2.fastq.gz: The secpmd pair of Illumina reads from paired-end sequencing
    Output: The file name to use for the output directory and the final assembly
    Threads: The number of threads to use
    unicycler_threads: The number of threads to use with Flye
    Template.yamp: The template.yaml file needed for PGAP annotation
    Template_submol.yaml: The template_submol.yaml file needed for PGAP annotation
    "
  exit 0
fi
if [ "$1" == "--help" ]; then
  echo "
  Usage: `basename $0` nanopore.fastq.gz illumina-1.fastq.gz illumina-2.fastq.gz output threads unicycler_threads

  Where:
    nanopore.fastq.gz: The gzipped basecalled nanopore data (uncorrected)
    illumina-1.fastq.gz: The first pair of Illumina reads from paired-end sequencing
    illumina-2.fastq.gz: The secpmd pair of Illumina reads from paired-end sequencing
    Output: The file name to use for the final assembly
    Threads: The number of threads to use
    unicycler_threads: The number of threads to use with Flye
    Template.yamp: The template.yaml file needed for PGAP annotation
    Template_submol.yaml: The template_submol.yaml file needed for PGAP annotation
    "
  exit 0
fi

# Check for Output directory
if [ ! -d "${4}" ]; then
  mkdir "${4}"
fi

# Illumina QC
mkdir "${4}/illumina_qc"
bbduk.sh in=${2} in2=${3} ref=adapters,artifacts,phix,lambda out="${4}/illumina_qc/forward.bbduk.fastq.gz" out2="${4}/illumina_qc/reverse.bbduk.fastq.gz"
java -jar /home/Bioinformatics_programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads $5 "${4}/illumina_qc/forward.bbduk.fastq.gz" "${4}/illumina_qc/reverse.bbduk.fastq.gz" "${4}/illumina_qc/forward.bbduk.trimmed.fastq.gz" "${4}/illumina_qc/forward.bbduk.unpaired.fastq.gz" "${4}/illumina_qc/reverse.bbduk.trimmed.fastq.gz" "${4}/illumina_qc/reverse.bbduk.unpaired.fastq.gz" LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Perform the initial assembly
mkdir "${4}/unicycler_assembly"
unicycler-runner.py -1 "${4}/illumina_qc/forward.bbduk.trimmed.fastq.gz" -2 "${4}/illumina_qc/reverse.bbduk.trimmed.fastq.gz" -l $1 -o  "${4}/unicycler_assembly" -t $6

# Check for overlaps on ends and remove short contigs
mkdir "${4}/overlap_removal_1"
nucmer --maxmatch --nosimplify -p "${4}/overlap_removal_1/assembly_overlaps" "${4}/unicycler_assembly/assembly.fasta" "${4}/unicycler_assembly/assembly.fasta"
show-coords -lrcTH "${4}/overlap_removal_1/assembly_overlaps.delta" > "${4}/overlap_removal_1/assembly_overlaps.coord"
parseCoord-nanopore.pl "${4}/overlap_removal_1/assembly_overlaps.coord" > "${4}/overlap_removal_1/assembly_overlaps.ToTrim.coord"
removeEnds-nanopore.pl "${4}/overlap_removal_1/assembly_overlaps.ToTrim.coord" "${4}/unicycler_assembly/assembly.fasta" > "${4}/overlap_removal_1/trimmed_assembly.fasta"
mv "${4}/overlap_removal_1/trimmed_assembly.fasta" "${4}/overlap_removal_1/trimmed_assembly.full.fasta"
pullseq -i "${4}/overlap_removal_1/trimmed_assembly.full.fasta" -m 2000 > "${4}/overlap_removal_1/trimmed_assembly.fasta"

# Perform medaka polishing using the nanopore data
mkdir "${4}/medaka_correction"
medaka_consensus -i ${1} -d "${4}/overlap_removal_1/trimmed_assembly.fasta" -o "${4}/medaka_correction" -t $5 -m r941_min_sup_g507

# Perform polypolish polishing using the illumina data
mkdir "${4}/polypolish_polishing"
bwa index "${4}/medaka_correction/consensus.fasta"
bwa mem -t $5 -a "${4}/medaka_correction/consensus.fasta" "${4}/illumina_qc/forward.bbduk.trimmed.fastq.gz" > "${4}/polypolish_polishing/forward_mapped.sam"
bwa mem -t $5 -a "${4}/medaka_correction/consensus.fasta" "${4}/illumina_qc/reverse.bbduk.trimmed.fastq.gz" > "${4}/polypolish_polishing/reverse_mapped.sam"
polypolish_insert_filter.py --in1 "${4}/polypolish_polishing/forward_mapped.sam" --in2 "${4}/polypolish_polishing/reverse_mapped.sam" --out1 "${4}/polypolish_polishing/forward_mapped.filtered.sam" --out2 "${4}/polypolish_polishing/reverse_mapped.filtered.sam"
polypolish "${4}/medaka_correction/consensus.fasta" "${4}/polypolish_polishing/forward_mapped.filtered.sam" "${4}/polypolish_polishing/reverse_mapped.filtered.sam" > "${4}/polypolish_polishing/polished_assembly.fasta"

# Perform polca polishing using the illumina data
mkdir "${4}/polca_polishing"
polca.sh -t $5 -a "${4}/polypolish_polishing/polished_assembly.fasta" -r ""${4}/illumina_qc/forward.bbduk.trimmed.fastq.gz" "${4}/illumina_qc/reverse.bbduk.trimmed.fastq.gz""
mv polished_assembly* "${4}/polca_polishing"
mv bwa.err "${4}/polca_polishing"
mv samtools.err "${4}/polca_polishing"

# Check for overlaps on ends
mkdir "${4}/overlap_removal_2"
nucmer --maxmatch --nosimplify -p "${4}/overlap_removal_2/assembly_overlaps" "${4}/polca_polishing/polished_assembly.fasta.PolcaCorrected.fa" "${4}/polca_polishing/polished_assembly.fasta.PolcaCorrected.fa"
show-coords -lrcTH "${4}/overlap_removal_2/assembly_overlaps.delta" > "${4}/overlap_removal_2/assembly_overlaps.coord"
parseCoord-nanopore.pl "${4}/overlap_removal_2/assembly_overlaps.coord" > "${4}/overlap_removal_2/assembly_overlaps.ToTrim.coord"
removeEnds-nanopore.pl "${4}/overlap_removal_2/assembly_overlaps.ToTrim.coord" "${4}/polca_polishing/polished_assembly.fasta.PolcaCorrected.fa" > "${4}/overlap_removal_2/trimmed_assembly.fasta"

# Reorient the replicons
mkdir "${4}/circlator_output"
circlator fixstart "${4}/overlap_removal_2/trimmed_assembly.fasta" final_assembly
mv final_assembly* "${4}/circlator_output"

# Get final assembly and compress everything
cp "${4}/circlator_output/final_assembly.fasta" "${4}/${4}.fasta"
pigz -r -p $5 "${4}/unicycler_assembly"
pigz -r -p $5 "${4}/overlap_removal_1"
pigz -r -p $5 "${4}/medaka_correction"
pigz -r -p $5 "${4}/illumina_qc"
pigz -r -p $5 "${4}/polypolish_polishing"
pigz -r -p $5 "${4}/polca_polishing"
pigz -r -p $5 "${4}/overlap_removal_2"
pigz -r -p $5 "${4}/circlator_output"

# Annotate the assembly
mkdir "${4}/pgap_annotation"
cp $7 "${4}/pgap_annotation/template.yaml"
cp $8 "${4}/pgap_annotation/"
cp "${4}/${4}.fasta" "${4}/pgap_annotation/"
cd "${4}/pgap_annotation/"
pgap.py --no-self-update --report-usage-true -c $5 -m 20g template.yaml
rename "s/annot/${4}/" output/*
cd ../..

