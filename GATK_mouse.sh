#!/bin/bash

module load bwa
module load samtools
module load novocraft
module load GATK

ref=/GATK_resource_bundle/mm10/mm10.fa

read1=/data/01.RawData/FN2/FN2_CKDN240019070-1A_22HVVTLT4_L7_1.fq.gz
read2=/data/01.RawData/FN2/FN2_CKDN240019070-1A_22HVVTLT4_L7_2.fq.gz

# alignment
bwa mem -t 4 -B 4 -O 6 -E 1 -M -p -R  "@RG\\tID:FN2\\tSM:FN2\\tLB:FN2\\tPL:ILLUMINA" /igenomes/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/genome.fa $read1 $read2 \
      | samtools view -Sb -h - >  /data//FN2.unsorted.bam

# sorting
# submit using sbatch --cpus-per-task=16 --mem=400g  F2_novo.sh
novosort -c 16 -m 20G --markDuplicates --index /data//FN2.unsorted.bam \
       -o /data//FN2.sorted.bam

# GATK BaseRecalibrator - generates a table
# submit using sbatch --cpus-per-task=16 --mem=24g --gres=lscratch:10 --time 128:00:00 FN2_recal.sh
cd /lscratch/$SLURM_JOB_ID || exit 1
bam=/data//FN2.sorted.bam

gatk --java-options "-Xmx16g -XX:ParallelGCThreads=2 -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
    BaseRecalibrator \
    -R $ref \
    -I $bam \
    --known-sites /GATK_resource_bundle/mm10/dbsnp146_fixedNames.vcf.gz \
    -O FN2.table

# GATK Apply BaseRecalibrator
# submit using sbatch --cpus-per-task=16 --mem=24g --gres=lscratch:10 --time 128:00:00 FN2_recal_apply.sh
cd /lscratch/$SLURM_JOB_ID || exit 1
recal_file=/data//FN2.table

gatk --java-options "-Xmx16g -XX:ParallelGCThreads=2 -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
    ApplyBQSR \
    -R $ref \
    -I $bam \
    --bqsr-recal-file $recal_file \
    -O /data//FN2.sorted.recal.bam

# GATK Joint Genotyping (HaplotypeCaller) on each individual BAM files to generate GVCF
# contains a record for each position in the genome whether there is a variant call or not in there
bam=/data//FN2.sorted.recal.bam
cd /lscratch/$SLURM_JOB_ID || exit 1
gatk --java-options "-Xmx16g -XX:ParallelGCThreads=2 -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
    HaplotypeCaller \
    -R $ref \
    -I $bam \
    -O /data//FN2.g.vcf.gz \
    -ERC GVCF

# GATK GenomicsDBImport to merge GVCFs from multiple samples
# creates a workspace 'db'
cd /lscratch/$SLURM_JOB_ID || exit 1

gatk --java-options "-Xmx16g -XX:ParallelGCThreads=2 -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
    GenomicsDBImport \
    -R $ref \
    --sample-name-map /data//cohort.sample_map \
    --genomicsdb-workspace-path /data//db \
    -L intervals.list
# where intervals.list is:
# chr1
# chr2
# chrX
# chr3
# chr4
# chr5
# chr6
# chr7
# chr10
# chr8
# chr14
# chr9
# chr11
# chr13
# chr12
# chr15
# chr16
# chr17
# chrY
# chr18
# chr19

# cohort.sample_map is:
#
#FN1 /data//FN1.g.vcf.gz
#FN2 /data//FN2.g.vcf.gz
#FN3 /data//FN3.g.vcf.gz
#FN4 /data//FN4.g.vcf.gz
#FN5 /data//FN5.g.vcf.gz
#NN1 /data//NN1.g.vcf.gz
#NN2 /data//NN2.g.vcf.gz
#NN3 /data//NN3.g.vcf.gz
#NN4 /data//NN4.g.vcf.gz
#NN5 /data//NN5.g.vcf.gz

# GATK GenotypeGVCFs to merge GVCFs from multiple samples
# where joint genotyping is taking place
# run for each chromosome
ref=/GATK_resource_bundle/mm10/mm10.fa

cd /lscratch/$SLURM_JOB_ID || exit 1
db=/data//db

gatk --java-options "-Xmx16g -XX:ParallelGCThreads=2 -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
    GenotypeGVCFs \
    -R $ref \
    --use-new-qual-calculator true \
    --max-alternate-alleles 4 \
    -V gendb:////data//db2/chr1\$1\$195471971/ \
    -O /data//chr1.vcf.gz

# Separate Indels and SNP, INDEL and MIXED
ref=/GATK_resource_bundle/mm10/mm10.fa

cd /lscratch/$SLURM_JOB_ID || exit 1

gatk --java-options "-Xmx16g -XX:ParallelGCThreads=2 -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
    SelectVariants \
    -R $ref \
    -V /data//output.vcf.gz \
    --select-type-to-include SNP \
    -O /data//output.rawsnps.vcf.gz

# VariantFiltration
ref=/GATK_resource_bundle/mm10/mm10.fa

cd /lscratch/$SLURM_JOB_ID || exit 1

gatk --java-options "-Xmx16g -XX:ParallelGCThreads=2 -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
    VariantFiltration \
    -R $ref \
    -V /data//output.rawsnps.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O /data//output.filteredsnps.vcf.gz

gatk --java-options "-Xmx16g -XX:ParallelGCThreads=2 -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
    VariantFiltration \
    -R $ref \
    -V /data//output.rawindels.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-2" \
    -O /data//output.filteredindels.vcf.gz

gatk --java-options "-Xmx16g -XX:ParallelGCThreads=2 -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
    VariantFiltration \
    -R $ref \
    -V /data//output.rawmixed.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-2" \
    -O /data//output.filteredmixed.vcf.gz
    

# swarm -f vep_GRCm38.swarm -t 24 -g 100 --time 128:00:00 --module VEP/102
# VEP
module load VEP/102

vep --offline --cache --dir ${VEP_CACHEDIR} --input_file /data//output.filteredsnps.vcf.gz --species mouse --assembly GRCm38 --fasta ${VEP_CACHEDIR}/mouse.fa --vcf --output_file output.filteresnps.vepped.vcf.gz --everything --compress_output bgzip
vep --offline --cache --dir ${VEP_CACHEDIR} --input_file /data//output.filteredindels.vcf.gz --species mouse --assembly GRCm38 --fasta ${VEP_CACHEDIR}/mouse.fa --vcf --output_file output.filteredindels.vepped.vcf.gz --everything --compress_output bgzip
vep --offline --cache --dir ${VEP_CACHEDIR} --input_file /data//output.filteredmixed.vcf.gz --species mouse --assembly GRCm38 --fasta ${VEP_CACHEDIR}/mouse.fa --vcf --output_file output.filteredmixed.vepped.vcf.gz --everything --compress_output bgzip
