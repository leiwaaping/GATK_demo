# GATK_demo
a gatk4 pipeline with NA12878 demo dateset ,conda environment used(no GPU)

### datasat  
PE100 rawdata from NA12878

### mode   
single sample and single lane, with only 2 fastq file (you can also modify the script to run with 1 fastq file)  
single sample and multi lane,with more than 2 fastq file.

### software  
fastq  
Soapnuke 2.0  
trimmomatic  
multiqc  
samtools  
bwa  
gatk4  
annovar

### related documents 
**REFERENCE**  
hg38：http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/  

```
#pre-procesing  
$ bwa index -a bwtsw /path/to/hg38.fa  
$ samtools faidx hg38.fa  
$ picard CreateSequenceDictionary R=/path/to/hg38.fa O=hg38.dict  
```
**GATK**  
gatkdoc：https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0  
```
# download files and index first
$ gatk IndexFeatureFile -F *_Mills_and_1000G_gold_standard.indels.hg38.vcf
```

**ANNOTATION**  
```
# download annotated documents and save in annovar/humandb/  
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/  
annotate_variation.pl -buildver hg38 -downdb cytoBand humandb/  
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar exac03 humandb/  
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar avsnp147 humandb/  
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbnsfp30a humandb/  

# pre-pracessing  
annotate_variation.pl -geneanno -dbtype refGene -buildver hg38 example/ex1.avinput humandb/  
annotate_variation.pl -regionanno -dbtype cytoBand -buildver hg38 example/ex1.avinput humandb/  
annotate_variation.pl -filter -dbtype exac03 -buildver hg38 example/ex1.avinput humandb/  
```

### LOAD IN CONDA ENVIRONMANT AND RUN SCRIPTS  
```
conda env create -n gatk4 -f condaenv_gatk.yaml
source activate gatk4
bash gatk4_2mode.sh |tee ./out.log
```

output info will saved in out.log ,and you should change MULTILANE = 0 if you are running single sample single lane mode.  

