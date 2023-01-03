#!/bin/bash
#
#SBATCH --job-name binner
#SBATCH --cpus-per-task=24
#SBATCH --mem=200GB
#SBATCH --output=binner-%j.out
#SBATCH --error=binner-%j.err
#SBATCH --partition=medium
#SBATCH --time=24:00:00 
#sbatch .sh <samplename>
#do sh #run in assembly folder with reads in ../*fastq.gz

#save current working directory
WD="$PWD"

module load anaconda3
source $ANACONDA3_ROOT/etc/profile.d/conda.sh 
conda activate /usr/users/benedict.yuen/conda_env/bbmap



#remove all contigs below 1000

reformat.sh in=scaffolds.fasta out=contigs-fixed.fa minlength=1000

#simplify names
sed -i 's/NODE_/c/g' contigs-fixed.fa
sed -i 's/_length.*//g' contigs-fixed.fa

rsync -v -L contigs-fixed.fa ../*fastq.gz $TMP_SCRATCH

#map reads to scaffolds and convert to bams, then sort and index the bam file

for f in $TMP_SCRATCH/*fastq.gz; do bbmap.sh -Xmx100g threads=$SLURM_CPUS_PER_TASK interleaved=auto unpigz=t ref=$TMP_SCRATCH/contigs-fixed.fa nodisk in1=$TMP_SCRATCH/${f##*/} out=$TMP_SCRATCH/${1}_${f##*/}.sam; done

conda deactivate

module load samtools/1.12

for f in $TMP_SCRATCH/*fastq.gz; do samtools view -b $TMP_SCRATCH/${1}_${f##*/}.sam > $TMP_SCRATCH/${1}_${f##*/}.bam; done 

for f in $TMP_SCRATCH/*fastq.gz; do samtools sort -o $TMP_SCRATCH/${1}_${f##*/}_sorted.bam $TMP_SCRATCH/${1}_${f##*/}.bam; done 

for f in $TMP_SCRATCH/*fastq.gz; do samtools index $TMP_SCRATCH/${1}_${f##*/}_sorted.bam; done 

echo 'bams made, sorted, and indexed'

rsync -v -L $TMP_SCRATCH/${1}_${f##*/}_sorted.bam* $WD

#run metabat wrapper script and metabat without coverage info
conda activate /scratch/projects/eei/software/metabat2

runMetaBat.sh --minContig 1500 -t $SLURM_CPUS_PER_TASK $TMP_SCRATCH/contigs-fixed.fa $TMP_SCRATCH/*_sorted.bam

metabat2 --minContig 1500 -t $SLURM_CPUS_PER_TASK -i $TMP_SCRATCH/contigs-fixed.fa -o $WD/metabat_nocov/${1}_nocov

conda deactivate

#run concoct

conda activate /scratch/projects/eei/software/concoct

cut_up_fasta.py $TMP_SCRATCH/contigs-fixed.fa -c 10000 -o 0 --merge_last -b $TMP_SCRATCH/contigs_10K.bed > $TMP_SCRATCH/contigs_10K.fa
concoct_coverage_table.py $TMP_SCRATCH/contigs_10K.bed $TMP_SCRATCH/*_sorted.bam >  $WD/coverage_table.tsv
concoct --threads $SLURM_CPUS_PER_TASK --composition_file $TMP_SCRATCH/contigs_10K.fa --coverage_file  $WD/coverage_table.tsv -b $WD/concoct_output/
merge_cutup_clustering.py $WD/concoct_output/clustering_gt1000.csv > $WD/concoct_output/clustering_merged.csv
mkdir $WD/concoct_output/fasta_bins
extract_fasta_bins.py $TMP_SCRATCH/contigs-fixed.fa $WD/concoct_output/clustering_merged.csv --output_path $WD/concoct_output/fasta_bins

conda deactivate


