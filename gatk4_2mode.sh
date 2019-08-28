#! /bin/bash
doc=data/         #input directory
out=multilane/
ref=/data/all_data/ref/hg38/   #reference directory
prefix=sample1     #sample name
gatkdoc=/data/all_data/GATK/gatkdoc/
annovar=/home/XDye/annovar/

#set mode
MULTILANE=0

if [ $MULTILANE -eq 1 ]
then
#multilanes mode
echo "MULTILANE mode processing"
echo "input files from :${doc}"
echo "output files in directory : ${out}"
echo "output prefix is : ${prefix}"
mkdir -p ${out}/clean ${out}/fastqc ${out}/bwa ${out}/annotation

#Quality contral
# trimmomatic
for i in ${doc}/V300014518_*_*_1.fq;do fq1=${i} fq2=${i::-4}2.fq fq1trim=${fq1::-2}clean.fq.gz fq2trim=${fq2::-2}clean.fq.gz fq1unpaired=${fq1::-2}unpaired.fq.gz fq2unpaired=${fq2::-2}unpaired.fq.gz;trimmomatic PE -threads 6 -phred33 $fq1 $fq2 $fq1trim $fq1unpaired $fq2trim $fq2unpaired HEADCROP:10 ;done && rm ${doc}/*.unpaired.* && echo "**trimming done**"

#quality check
time fastqc -o ${out}/fastqc/ ${doc}/*.clean.fq.gz -t 6
time multiqc ${out}/fastqc/*.zip -o ${out}/fastqc/

# alignment
for i in ${doc}/V300014518_*_*_1.clean.fq.gz;do fq1=${i} fq2=${i::-13}2.clean.fq.gz fq1fn=${fq1##*/} fq1bam=${fq1fn::-14}.bam laneid=${fq1bam#*_} bwaid=${laneid%.*} label="@RG\tID:${bwaid}\tSM:V300021733\tLB:WES\tPL:COMPLETE"; time bwa mem -t 6 -R $label ${ref}/hg38.fa $fq1 $fq2 | samtools view -Sb - > ${out}/bwa/${fq1bam};done && echo "** bwa mapping done **"

#merge bam files
time samtools merge -@ 6  ${out}${prefix}.bam ${out}/bwa/*.bam && rm -f ${out}/bwa/*.bam && mv ${out}${prefix}.bam ${out}/bwa/ && echo "** bam merge done**"

#single lane mode
else
echo "SINGAL LANE mode processing"
echo "input files from :${doc}"
echo "output files in directory : ${out}"
echo "output prefix is : ${prefix}"
mkdir -p ${out}/clean ${out}/fastqc ${out}/bwa ${out}/annotation

#Quality contral
# SOAPnuke or trimmomatic
time SOAPnuke2.0 filter -1 ${doc}sample1_1.fq -2  ${doc}sample1_2.fq -T 6 -l 5 -q 0.5 -n 0.1 -f AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA -rAAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG -Q 2 -G 2 --seqType 0 -o ${out}/clean/ -C ${prefix}_1.clean.fq.gz -D ${prefix}_2.clean.fq.gz  && echo "**SOAPnuke done **"
trimmomatic PE -threads 6 -phred33 ${doc}/sample1_1.fq ${doc}/sample1_2.fq ${out}/clean/${prefix}_1.single.clean.trim.fq.gz ${out}/clean/${prefix}_1.clean.unpaired.fq.gz ${out}/clean/${prefix}_2.single.clean.trim.fq.gz ${out}/clean/${prefix}_2.clean.unpaired.fq.gz HEADCROP:10 && rm ${out}/clean/*.unpaired.fq.gz  && echo "**trimmomatic done **"

#quality check
time fastqc -o ${out}/fastqc/ ${out}/clean/*.fq.gz

# alignment
time bwa mem -t 6 -R "@RG\tID:sample1\tSM:sample1\tLB:WES\tPL:illumina" ${ref}/hg38.fa ${out}/clean/${prefix}_1.single.clean.trim.fq.gz ${out}/clean/${prefix}_2.single.clean.trim.fq.gz | samtools view -Sb - > ${out}/bwa/${prefix}.bam && echo "**bwa mapping done **"

fi

# sam/bam file processing  
#if error： open too much file $ ulimit -n 10240
#if error: GC-overhead limit exceeded ; use option -Xmx64g，64g for 64G memory，
time samtools sort -@ 6 -o ${out}${prefix}.sorted.bam ${out}/bwa/${prefix}.bam && echo "** bam sort done **"
time picard MarkDuplicates -Xmx64g I=${out}/${prefix}.sorted.bam O=${out}${prefix}.sorted.markdup.bam M=${out}${prefix}.sorted.markdup.txt REMOVE_DUPLICATES=true && rm -f ${out}${prefix}.sorted.bam && echo "** picard Markdup done **"
time picard BuildBamIndex -Xmx64g I=${out}${prefix}.sorted.markdup.bam && echo "** picard bam index done **"

#BQSR
time gatk BaseRecalibrator -R ${ref}/hg38.fa -I ${out}${prefix}.sorted.markdup.bam \
--known-sites ${gatkdoc}/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
--known-sites ${gatkdoc}/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf \
--known-sites ${gatkdoc}/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf \
-O ${out}/sample1.BQSR.table && echo "** bqsr table done **" 
time gatk ApplyBQSR -R ${ref}/hg38.fa -I ${out}${prefix}.sorted.markdup.bam -bqsr ${out}/${prefix}.BQSR.table -O ${out}/${prefix}.sorted.markdup.BQSR.bam && echo "** BQSR done **"

# candidate variant calling
time gatk HaplotypeCaller -R ${ref}/hg38.fa -I ${out}/${prefix}.sorted.markdup.BQSR.bam -O ${out}/${prefix}.HC.vcf && echo "** Haplotype variant calling done **"

# VCF annotation by ANNOVAR
time ${annovar}table_annovar.pl ${out}/${prefix}.HC.vcf ${annovar}humandb/ -buildver hg38 -out ${out}/annotation/${prefix} -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a -operation gx,r,f,f,f -nastring . -vcfinput -polish && echo "** annotation done **"
