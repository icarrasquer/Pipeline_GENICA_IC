#!/bin/bash
#SBATCH --job-name="Test_cleaning"	
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu	
#SBATCH --account=paygo
#SBATCH --wckey=ips_genica				
#SBATCH --output=jobs/job_test_%j.out
#SBATCH --error=jobs/job_test_%j.err

# Script written by Ines Carrasquer

module load Anaconda3
source $(conda info --base)/etc/profile.d/conda.sh
module load cutadapt/3.4-GCCcore-10.3.0
module load FastQC/0.11.9-Java-11 
module load Perl
module load gzip



for i in D1_new_
    do
    for j in R1 R2
        do

        ##Merge all reads
        zcat ../GENICA/${i}/${i}_*${j}* ../GENICA_2/${i}/${i}_*${j}* > temp_${i}${j}.fastq
        nb_reads=$(grep -c "." temp_${i}${j}.fastq | awk '{print $1/4}')
        nb_Rd1_SP=$(grep -c TCCGATCT temp_${i}${j}.fastq)
        nb_Rd2_SP=$(grep -c AGATCGGA temp_${i}${j}.fastq)
        echo temp_${i}${j}.fastq ${nb_reads} ${nb_Rd1_SP} ${nb_Rd2_SP} >> summary_stats
        ##Stats
        conda activate seqkit
        seqkit stats temp_${i}${j}.fastq >> summary_lenght
        conda deactivate

        ##Adapter trimming
        cutadapt -j 4 -e 1 -a AGATCGGAAG temp_${i}${j}.fastq -o temp_cutadapt1_${i}${j}.fastq -m 30
        num=$(grep -c CTTCCGATCT temp_cutadapt1_${i}${j}.fastq)

        while (( $num > 100 )) 
            do
            cutadapt -j 4 -e 1 -g CTTCCGATCT temp_cutadapt1_${i}${j}.fastq -o temp_loop_${i}${j}.fastq -m 30
            mv temp_loop_${i}${j}.fastq temp_cutadapt1_${i}${j}.fastq
            num=$(grep -c CTTCCGATCT temp_cutadapt1_${i}${j}.fastq)
            done

        ##PolyG
        DropBpFastq_polyC.pl temp_cutadapt1_${i}${j}.fastq temp_cutadapt1_polyG_${i}${j}.fastq

        ##Quality & lenght 
        cutadapt -j 4 -q 30 temp_cutadapt1_polyG_${i}${j}.fastq -o clean_${i}${j}.fastq -m 30
        nb_reads=$(grep -c "." clean_${i}${j}.fastq | awk '{print $1/4}')
        nb_Rd1_SP=$(grep -c TCCGATCT clean_${i}${j}.fastq)
        nb_Rd2_SP=$(grep -c AGATCGGA clean_${i}${j}.fastq)
        echo clean_${i}${j}.fastq ${nb_reads} ${nb_Rd1_SP} ${nb_Rd2_SP} >> summary_stats
        
        ##Stats
        conda activate seqkit
        seqkit stats clean_${i}${j}.fastq >> summary_lenght
        conda deactivate

        fastqc clean_${i}${j}.fastq -o FASTQC/
        gzip clean_${i}${j}.fastq
        #rm temp*
        done
    conda activate seqkit
    seqkit pair \
    -1 clean_${i}R1.fastq.gz \
    -2 clean_${i}R2.fastq.gz \
    nb_reads=$(zgrep -c "." clean_${i}R1.paired.fastq.gz | awk '{print $1/4}')
    echo clean_${i}R1.paired.fastq.gz ${nb_reads} ${nb_Rd1_SP} ${nb_Rd2_SP} >> summary_stats


    # conda deactivate
    conda activate adapterremoval

    AdapterRemoval \
    --file1 clean_${i}R1.paired.fastq.gz \
    --file2 clean_${i}R2.paired.fastq.gz \
    --gzip \
    --threads 4 \
    --trimns \
    --trimqualities \
    --minquality 30 \
    --minlength 30 \
    --collapse \
    --outputcollapsed clean_${i}.collapsed.fastq.gz 
    rm your_output*

    nb_reads=$(zgrep -c "." clean_${i}.collapsed.fastq.gz | awk '{print $1/4}')
    echo clean_${i}.collapsed.fastq.gz ${nb_reads} ${nb_Rd1_SP} ${nb_Rd2_SP} >> summary_stats

    ##Stats
    conda activate seqkit
    seqkit stats clean_${i}.collapsed.fastq.gz >> summary_lenght
    conda deactivate

    done

    

