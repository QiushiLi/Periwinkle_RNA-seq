#!/bin/bash
# ---------------------------------------------------------------------
# SLURM script for a pilon job on Helix.
# ---------------------------------------------------------------------
#SBATCH --partition=bigmem
#SBATCH --time=24:00:00
#SBATCH --mem=1000G
#SBATCH -n 64
#SBATCH --nodes=1
#SBATCH --job-name=peri_pilon_job1_p
#SBATCH --output=peri_pilon_job1_p.out
#SBATCH --error=peri_pilon_job1_p.err

##########
# Setup
##########

echo "Current working directory is `pwd`"

# Load modules

module load java

# activate conda environment

source ~/program/anaconda3/etc/profile.d/conda.sh
conda activate pilon

# Need date for record-keeping and unique identifier
DATE="`date +%b%d_%T`"
echo "Run start time: `date +%b%d_%T`"
# pilon before phasing with Hi-C from /home/qiushi.li/periwinkle/job8-ovl1/4-polish/cns-output
bwa index cns_p_ctg.fasta -p peri_job1_cns_p

bwa mem -t 64 peri_job1_cns_p Catharanthus_roseus_6678_4000ng_R1_clean.fastq Catharanthus_roseus_6678_4000ng_R2_clean.fastq > peri_job1_cns_p_vs_illumina6678.sam 2> bwa_mem_p.log

samtools view -bS -@ 64 peri_job1_cns_p_vs_illumina6678.sam > peri_job1_cns_p_vs_illumina6678.bam

samtools sort -@ 64 peri_job1_cns_p_vs_illumina6678.bam -o peri_job1_cns_p_vs_illumina6678_sorted.bam

samtools index peri_job1_cns_p_vs_illumina6678_sorted.bam peri_job1_cns_p_vs_illumina6678_sorted.bam.bai

java -Xmx800G -jar /home/qiushi.li/program/anaconda3/envs/pilon/share/pilon-1.23-2/pilon-1.23.jar --genome cns_p_ctg.fasta --threads 64 --frags peri_job1_cns_p_vs_illumina6678_sorted.bam --output peri_pilon_job1_p --outdir ./ --changes --tracks --diploid --fix bases --mindepth 10 --minmq 30

echo ""
echo "Run end time: `date +%b%d_%T`"

