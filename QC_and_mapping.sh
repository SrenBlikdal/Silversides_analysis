#!/bin/bash

#SBATCH --job-name=D1Gen5
#SBATCH --output=D1Gen5.out
#SBATCH --error=D1Gen5.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=8G    # memory per cpu-core
#SBATCH --time=30:00:00
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=soren.blikdal.hansen@sund.ku.dk
### Script below ###

module load anaconda3/2021.11
module load perl/5.24.4
module load pigz/2.6.0
module load openjdk/13.0.1
module load fastqc/0.12.1
module load trimgalore/0.6.6
module load samtools/1.18 
module load bowtie2/2.5.2
module load bismark/0.24.0

# change the file location accordingly
raw_data=/projects/mjolnir1/people/vrt191/Silversides_feb2024/1_raw_data
out_direc=/projects/mjolnir1/people/vrt191/Silversides_feb2024/
ref_genome=/projects/mjolnir1/people/vrt191/Silversides_feb2024/pop_ref/D1

for SAMPLE in D1Gen5_121  D1Gen5_237  D1Gen5_252  D1Gen5_266  D1Gen5_195  D1Gen5_239  D1Gen5_256  D1Gen5_339  D1Gen5_200  D1Gen5_240  D1Gen5_260  D1Gen5_342  D1Gen5_236  D1Gen5_245  D1Gen5_264

do

        mkdir -p $out_direc/2_fastqc_before/$SAMPLE
        mkdir -p $out_direc/3_trimgalore/$SAMPLE
        mkdir -p $out_direc/4_fastqc_after/$SAMPLE
        mkdir -p $out_direc/5_bismark/$SAMPLE

        # fastqc - before trimming
        cd $raw_data/$SAMPLE
        fastqc "$SAMPLE"_1.fq.gz "SAMPLE"_2.fq.gz -t 10 -o $out_direc/2_fastqc_before/$SAMPLE

        # Trimgalore (Adapterpair 2 seems to be best)
        cd $out_direc/3_trimgalore/$SAMPLE
        trim_galore -j 10  --path_to_cutadapt /projects/mjolnir1/apps/conda/cutadapt-3.7/bin/cutadapt --clip_R2 10 --three_prime_clip_R2 2 --three_prime_clip_R1 2 --paired $raw_data/$SAMPLE/"$SAMPLE"_1.fq.gz $raw_data/$SAMPLE/"$SAMPLE"_2.fq.gz 

        # fastqc - after trimming
        cd $out_direc/3_trimgalore/$SAMPLE
        fastqc "$SAMPLE"_1.fq.gz "SAMPLE"_2.fq.gz -t 10 -o $out_direc/4_fastqc_after/$SAMPLE

        # bismark
        #To increase the mapping rate
        #https://github.com/FelixKrueger/Bismark/blob/master/Docs/FAQ.md#issue-2-low-mapping-effiency-of-paired-end-bisulfite-seq-sample
        cd $out_direc/3_trimgalore/$SAMPLE
        bismark --parallel 6 --genome $ref_genome --output_dir $out_direc/5_bismark/$SAMPLE --score_min L,0,-0.6 -1 *_1.fq.gz -2 *_2.fq.gz

        # deduplication
        cd $out_direc/5_bismark/$SAMPLE
        deduplicate_bismark --paired --bam *bam

        # methylation extraction (excluded "--counts"; excluded "--multi")
        cd $out_direc/5_bismark/$SAMPLE
        bismark_methylation_extractor -p --gzip --scaffolds --buffer_size 85% --parallel 6 --ignore 0 --ignore_3prime 0  --ignore_3prime_r2 0 --cytosine_report -o $out_direc/5_bismark/$SAMPLE --genome_folder $ref_genome  *deduplicated.bam
done
