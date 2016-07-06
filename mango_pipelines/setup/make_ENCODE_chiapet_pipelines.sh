ROOT_DIR=""

for i in `cut -f1 setup/EXPERIMENT_INFO.TXT | sort | uniq`
do

	PIPELINE_FILE="${ROOT_DIR}/pipeline_files/$i.pipeline.sh"

	printf "#!/bin/bash\n\n" > $PIPELINE_FILE	
	for j in `grep $i setup/EXPERIMENT_INFO.TXT | cut -f2`
	do
		echo "wget https://www.encodeproject.org/files/${j}/@@download/${j}.fastq.gz -O ${ROOT_DIR}/data/${j}.fastq.gz" >> $PIPELINE_FILE
	done

	for j in `grep $i setup/EXPERIMENT_INFO.TXT | awk '$3=="R1"' | cut -f2`
	do
		echo "zcat ${ROOT_DIR}/data/${j}.fastq.gz >> ${ROOT_DIR}/data/${i}.R1.fastq" >> $PIPELINE_FILE
		echo "rm ${ROOT_DIR}/data/${j}.fastq.gz" >> $PIPELINE_FILE
	done
	
	for j in `grep $i setup/EXPERIMENT_INFO.TXT | awk '$3=="R2"' | cut -f2`
	do
		echo "zcat ${ROOT_DIR}/data/${j}.fastq.gz >> ${ROOT_DIR}/data/${i}.R2.fastq" >> $PIPELINE_FILE
		echo "rm ${ROOT_DIR}/data/${j}.fastq.gz" >> $PIPELINE_FILE
	done
	
	mkdir ${ROOT_DIR}/mango_out/${i}_out

	echo "Rscript ${ROOT_DIR}/mango/mango.R --fastq1 ${ROOT_DIR}/data/${i}.R1.fastq --fastq2 ${ROOT_DIR}/data/$i.R2.fastq --prefix $i --chrominclude chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22 --bowtieref ${ROOT_DIR}/bowtie_index/hg19 --bedtoolsgenome ${ROOT_DIR}/ref/hg19.chrom.sizes --outdir ${ROOT_DIR}/mango_out/${i}_out --FDR 0.05 &> ${ROOT_DIR}/mango_out/${i}.stderr_stdout &" >> $PIPELINE_FILE
	
done
