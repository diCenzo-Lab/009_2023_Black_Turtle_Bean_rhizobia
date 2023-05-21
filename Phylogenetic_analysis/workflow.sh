# Make directories
mkdir Genome_files/
mkdir Genome_files_input/
mkdir Reannotate_genomes/
mkdir FastANI_output/
mkdir Phylogeny/

# Download genomes
perl Scripts/parseGenomeList.pl Input_files/Rhizobium_type_strains.txt # Parse the NCBI genome table to get info to download genomes
sed -i 's/_sp._/_sp_/' Input_files/genomeList.txt
perl Scripts/downloadGenomes.pl Input_files/genomeList.txt # Download the genomes of interest
cat Input_files/genomeList.txt Input_files/new_genomes.txt > temp.txt
mv temp.txt Input_files/genomeList.txt
sort -u Input_files/genomeList.txt > temp.txt
mv temp.txt Input_files/genomeList.txt
cp ../Rhizobium_*/*.fasta Genome_files
rename 's/fasta/fna/' Genome_files/*

# Calculate ANI
find Genome_files/*.fna > Input_files/genomePaths.txt # Get the genome paths
fastANI -q Genome_files/Rhizobium_leguminosarum_BT01.fna --rl Input_files/genomePaths.txt -o FastANI_output/Rhizobium_leguminosarum_BT01_fastani_output.txt
fastANI -q Genome_files/Rhizobium_sp_BT03.fna --rl Input_files/genomePaths.txt -o FastANI_output/Rhizobium_sp_BT03_fastani_output.txt
fastANI -q Genome_files/Rhizobium_sp_BT04.fna --rl Input_files/genomePaths.txt -o FastANI_output/Rhizobium_sp_BT04_fastani_output.txt

# Reannotate genomes
mv Genome_files/* Genome_files_input/
perl Scripts/runProkka.pl Input_files/genomeList.txt # Run prokka to annotate the genomes
perl Scripts/moveGenomes.pl Input_files/genomeList.txt # Collect important reannotated genome files

# Create phylogeny with IQ-TREE 2, with jackknife support values
roary -p 16 -f Roary_output -e -n -i 80 -g 300000 Genome_files/*.gff # Run roary
trimal -in Roary_output/core_gene_alignment.aln -out core_gene_alignment_trimmed.aln -fasta -automated1 # Trim the alignment made by Roary
mv core_gene_alignment_trimmed.aln Phylogeny/
cd Phylogeny/
iqtree2 -s core_gene_alignment_trimmed.aln -m MF -T 16 # Used results to identify best model for next step
iqtree2 -s core_gene_alignment_trimmed.aln -m GTR+F+I+I+R6 --alrt 1000 -J 1000 --jack-prop 0.4 -T 16 --prefix Rhizobium_phylogeny
cd ../
