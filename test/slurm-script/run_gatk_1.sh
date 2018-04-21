### GATK pipeline created by huskeR
### Thu May 04 15:50:38 2017

### Generate a SAM file containing aligned reads
bwa mem -M -R '@RG\tID:g1\tSM:s1\tPL:illumina\tLB:lib1\tPU:unit1' -T 5 -t 1 ~/dbcenter/Ecoli/reference/Ecoli_k12_MG1655.fasta fq_1.fq f1_2.fq > mysample.aln.sam

java -Xmx9g -jar $HOME/bin/picard-tools-2.1.1/picard.jar SortSam \
    INPUT=mysample.aln.sam \
    OUTPUT=mysample.sorted.bam \
    SORT_ORDER=coordinate
#rm mysample.aln.sam

### Mark duplicates
java -Xmx9g -jar $HOME/bin/picard-tools-2.1.1/picard.jar MarkDuplicates \
INPUT=mysample.sorted.bam \
OUTPUT=mysample.sorted.dedup.bam \
METRICS_FILE=mysample_metrics.txt

java -Xmx9g -jar $HOME/bin/picard-tools-2.1.1/picard.jar BuildBamIndex \
INPUT=mysample.sorted.dedup.bam

samtools flagstat mysample.sorted.bam > mysample.sorted.bam.log

### Call variants in your sequence data
java -Xmx9g -jar $HOME/bin/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R ~/dbcenter/Ecoli/reference/Ecoli_k12_MG1655.fasta \
-I mysample.sorted.dedup.bam \
-ERC GVCF \
-o mysample.g.vcf
