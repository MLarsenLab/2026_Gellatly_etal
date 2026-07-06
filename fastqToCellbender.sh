#!/bin/bash
#SBATCH --job-name=Gellatly2026	       		# Job name
#SBATCH --mail-type=END,FAIL           		# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --ntasks=1                      	# Run a single task
#SBATCH --mem=100gb                      	# Job memory request
#SBATCH --cpus-per-task=32               	# Number of CPU cores per task
#SBATCH --time=168:00:00                 	# Time limit hrs:min:sec
#SBATCH --output=Gellatly2026.log		# Standard output and error log

source ~/.bashrc
export PATH=../cellranger-8.0.0:$PATH
conda activate cellbender
mkdir cellranger
mkdir cellbender
for i in $(ls ./fastq | grep '_001.fastq.gz' | sed s/_S._L001_.._00..fastq.gz// | sort -u)
do
	cellranger count --id="$i" \
		--fastqs=./fastq \
		--sample="$i" \
   		--include-introns true \
		--create-bam true \
		--output-dir ./cellranger/$i \
   		--transcriptome=/network/rit/lab/larsenlab-rit/Next_Generation_Sequencing_Data/Reference_Genome/GRCm39-8-TdT-Cre
	echo $[j =  $(cat ./cellranger/"$i"/outs/metrics_summary.csv | head -2 | tail -1 | cut -d "," -f 1-2 | sed 's/[",]//g')]
	echo $[k = j + 10000]
	cellbender remove-background \
		--expected-cells "$j" \
		--total-droplets-included "$k" \
		--input "./cellranger/"$i"/outs/raw_feature_bc_matrix.h5" \
		--output "./cellbender/"$i"_cellbender.h5" \
		--estimator-multiple-cpu \
		--cpu-threads 32 \
		--learning-rate 0.0001 \
		--fpr 0.01 0.05 0.1 \
		--checkpoint-mins 30 \
		--debug
done

