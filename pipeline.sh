#!/bin/sh
#this script is a shell script to run through the whole mutation calling pipeline 
# input will be
# -s as sample name  ;
# -i as raw reads directory  --absolute directory here 
# -o as output_directory   --absolute directory here
# -d design file --a bed file indicating which regions are captured --majorly for QC purpose
# -n sample number(S1,S2,...depend on the name of sample file)
# -g gene list file for filtering   --absolute directory here
# -c SeqCNV control file  --absolute directory here
# -w run WES pipeline
# -f compressed format(miseq(default),nextseq,bz)
# -t frequency cutoff(default:0.005)
#shellPath=$(cd "$(dirname "$0")"; pwd)
shellPath=...
export PATH=${shellPath}/bin/jdk1.8.0_121/bin/:$PATH
fastpPath=${shellPath}/bin/fastp
bwaPath=${shellPath}/bin/bwa
samtoolsPath=${shellPath}/bin/samtools-1.2/samtools
gatkPath=${shellPath}/bin/gatk4/gatk
pythonPath=${shellPath}/python/python
geneFile="ig"
compressedFormat="miseq"
cutoff=0.005
while getopts ":s:i:o:d:n:g:c:w:f:t:p:v:x:" opt
do
    case $opt in
        s)
            sName=$OPTARG
        ;;
        i)
            iPath=$OPTARG
        ;;
        o)
            oPath=$OPTARG
        ;;
        d)
            designFile=$OPTARG
        ;;
        n)
            sNum=$OPTARG
        ;;
        g)
            geneFile=$OPTARG
        ;;
        c)
            CNVFile=$OPTARG
        ;;
        w)
            ifWES=$OPTARG
        ;;
        f)
            compressedFormat=$OPTARG
        ;;
        t)
            cutoff=$OPTARG
        ;;
        p)
            disease=$OPTARG
        ;;
        v)
             wgs=$OPTARG
        ;;
        x)
             svtyper_ok=$OPTARG
        ;;
        ?)
            echo "unknown argument"
            exit 1
        ;;
    esac
done
if [ ! ${sName} ]; then
    echo "-s invalid"
    exit 1
fi
if [ ! ${iPath} ]; then
    echo "-i invalid"
    exit 1
fi
if [ ! ${oPath} ]; then
    echo "-o invalid"
    exit 1
fi
if [ ! ${designFile} ]; then
    echo "-d invalid"
    exit 1
fi
if [ ! ${disease} ]; then
    echo "-p invalid"
    exit 1
fi

mkdir -p ${oPath}/${sName}
exec 1>>${oPath}/${sName}/${sName}_log.txt
exec 2>>${oPath}/${sName}/${sName}_log.txt

echo -e "${iPath}/${sName} startTime:::::::\c" ; date

echo "## Pipeline Starts to Run##"

touch ${oPath}/${sName}/${sName}_status.txt
date_start=$(date +%s)

echo "## Uncompress Fastq File ##" > ${oPath}/${sName}/${sName}_status.txt
################# 1. Unzipping fastq files #########################################################################
if [ ${compressedFormat} == "bz" ]; then
    bunzip2 -k ${iPath}/${sName}.read1.bz2
    bunzip2 -k ${iPath}/${sName}.read2.bz2
    mv ${iPath}/${sName}.read1 ${oPath}/${sName}/${sName}.read1.fq
    mv ${iPath}/${sName}.read2 ${oPath}/${sName}/${sName}.read2.fq
elif [ ${compressedFormat} == "nextseq" ];then
    gunzip -c ${iPath}/${sName}_${sNum}_L001_R1_001.fastq.gz > ${oPath}/${sName}/${sName}_${sNum}_L001_R1_001.fastq
    gunzip -c ${iPath}/${sName}_${sNum}_L002_R1_001.fastq.gz > ${oPath}/${sName}/${sName}_${sNum}_L002_R1_001.fastq
    gunzip -c ${iPath}/${sName}_${sNum}_L003_R1_001.fastq.gz > ${oPath}/${sName}/${sName}_${sNum}_L003_R1_001.fastq
    gunzip -c ${iPath}/${sName}_${sNum}_L004_R1_001.fastq.gz > ${oPath}/${sName}/${sName}_${sNum}_L004_R1_001.fastq
    gunzip -c ${iPath}/${sName}_${sNum}_L001_R2_001.fastq.gz > ${oPath}/${sName}/${sName}_${sNum}_L001_R2_001.fastq
    gunzip -c ${iPath}/${sName}_${sNum}_L002_R2_001.fastq.gz > ${oPath}/${sName}/${sName}_${sNum}_L002_R2_001.fastq
    gunzip -c ${iPath}/${sName}_${sNum}_L003_R2_001.fastq.gz > ${oPath}/${sName}/${sName}_${sNum}_L003_R2_001.fastq
    gunzip -c ${iPath}/${sName}_${sNum}_L004_R2_001.fastq.gz > ${oPath}/${sName}_${sNum}_L004_R2_001.fastq

    touch ${iPath}/${sName}.read1.fq
    cat ${iPath}/${sName}_${sNum}_L001_R1_001.fastq >> ${oPath}/${sName}/${sName}.read1.fq
    rm ${iPath}/${sName}_${sNum}_L001_R1_001.fastq
    cat ${iPath}/${sName}_${sNum}_L002_R1_001.fastq >> ${oPath}/${sName}/${sName}.read1.fq
    rm ${iPath}/${sName}_${sNum}_L002_R1_001.fastq
    cat ${iPath}/${sName}_${sNum}_L003_R1_001.fastq >> ${oPath}/${sName}/${sName}.read1.fq
    rm ${iPath}/${sName}_${sNum}_L003_R1_001.fastq
    cat ${iPath}/${sName}_${sNum}_L004_R1_001.fastq >> ${oPath}/${sName}/${sName}.read1.fq
    rm ${iPath}/${sName}_${sNum}_L004_R1_001.fastq

    touch ${iPath}/${sName}.read2.fq
    cat ${iPath}/${sName}_${sNum}_L001_R2_001.fastq >> ${oPath}/${sName}/${sName}.read2.fq
    rm ${iPath}/${sName}_${sNum}_L001_R2_001.fastq
    cat ${iPath}/${sName}_${sNum}_L002_R2_001.fastq >> ${oPath}/${sName}/${sName}.read2.fq
    rm ${iPath}/${sName}_${sNum}_L002_R2_001.fastq
    cat ${iPath}/${sName}_${sNum}_L003_R2_001.fastq >> ${oPath}/${sName}/${sName}.read2.fq
    rm ${iPath}/${sName}_${sNum}_L003_R2_001.fastq
    cat ${iPath}/${sName}_${sNum}_L004_R2_001.fastq >> ${oPath}/${sName}/${sName}.read2.fq
    rm ${iPath}/${sName}_${sNum}_L004_R2_001.fastq
elif [ -n "${compressedFormat}" ] && ([ ${compressedFormat} == "fastq" ] || [ ${compressedFormat} == "bwa" ] || [ ${compressedFormat} == "markup" ]|| [ ${compressedFormat} == "bam" ] || [ ${compressedFormat} == "recal" ]|| [ ${compressedFormat} == "realign" ]|| [ ${compressedFormat} == "gatk" ]|| [ ${compressedFormat} == "annov" ]|| [ ${compressedFormat} == "filter" ]|| [ ${compressedFormat} == "clean" ] || [ ${compressedFormat} == "qc" ] || [ ${compressedFormat} == "gatk" ] ||[ ${compressedFormat} == "gatkmark" ] || [ ${compressedFormat} == "cleanall" ]|| [ ${compressedFormat} == "force" ]||[ ${compressedFormat} == "test" ]||[ ${compressedFormat} == "exomiser" ]||[ ${compressedFormat} == "sv" ] && [ ${compressedFormat} == "mysql" ] && [ ${compressedFormat} == "spliceai" ]); then 
    echo "fastq ready"
else
    gunzip -c ${iPath}/${sName}_R1_001.fastq.gz > ${oPath}/${sName}/${sName}_R1_001.fastq
    gunzip -c ${iPath}/${sName}_R2_001.fastq.gz > ${oPath}/${sName}/${sName}_R2_001.fastq
    mv ${oPath}/${sName}/${sName}_R1_001.fastq ${oPath}/${sName}/${sName}.read1.fq
    mv ${oPath}/${sName}/${sName}_R2_001.fastq ${oPath}/${sName}/${sName}.read2.fq
fi

################# 2. trimming fastq files #########################################################################
cd ${shellPath}
echo -e '${shellPath} is :::::::::\c';echo "${shellPath}"
if [ -n "${compressedFormat}" ] && ([ ${compressedFormat} == "bwa" ] || [ ${compressedFormat} == "bam" ] ||  [ ${compressedFormat} == "markup" ]||[ ${compressedFormat} == "recal" ]|| [ ${compressedFormat} == "realign" ]|| [ ${compressedFormat} == "gatk" ]|| [ ${compressedFormat} == "annov" ]|| [ ${compressedFormat} == "filter" ]|| [ ${compressedFormat} == "clean" ] || [ ${compressedFormat} == "qc" ] || [ ${compressedFormat} == "gatk" ] ||[ ${compressedFormat} == "gatkmark" ] || [ ${compressedFormat} == "cleanall" ]|| [ ${compressedFormat} == "force" ]||[ ${compressedFormat} == "test" ]||[ ${compressedFormat} == "exomiser" ]||[ ${compressedFormat} == "sv" ] && [ ${compressedFormat} == "mysql" ] && [ ${compressedFormat} == "spliceai" ]); then
   echo "fastq ready"
elif [ -n "${compressedFormat}" ] && [ ${compressedFormat} == "bz" ]; then
   java -Xmx4g -classpath ${shellPath}/bin/SNPfilter/ FastqShaver  ${oPath}/${sName}/${sName}.read1.fq  ${shellPath}/reference/adapters/cap_RP/adapters.txt >${oPath}/${sName}/${sName}_1_sequence.sheared.txt
   java -Xmx4g -classpath ${shellPath}/bin/SNPfilter/ FastqShaver  ${oPath}/${sName}/${sName}.read2.fq  ${shellPath}/reference/adapters/cap_RP/adapters.txt > ${oPath}/${sName}/${sName}_2_sequence.sheared.txt
else
echo "## Run fastp ##" > ${oPath}/${sName}/${sName}_status.txt
${fastpPath} -i ${oPath}/${sName}/${sName}.read1.fq -I ${oPath}/${sName}/${sName}.read2.fq -o ${oPath}/${sName}/${sName}_1_sequence.sheared.txt -O ${oPath}/${sName}/${sName}_2_sequence.sheared.txt -h ${oPath}/${sName}/${sName}.fastp.html -j ${oPath}/${sName}/${sName}.fastp.json
fi
# remove fastq files after fastq trimming
rm ${iPath}/${sName}.read1.fq
rm ${iPath}/${sName}.read2.fq

################# 3. mapping fastq to bam files #########################################################################
if [ -n "${compressedFormat}" ] && ([ ${compressedFormat} == "bam" ] || [ ${compressedFormat} == "recal" ]|| [ ${compressedFormat} == "markup" ]|| [ ${compressedFormat} == "realign" ]|| [ ${compressedFormat} == "gatk" ]|| [ ${compressedFormat} == "annov" ]|| [ ${compressedFormat} == "filter" ]|| [ ${compressedFormat} == "clean" ] || [ ${compressedFormat} == "qc" ]) ||[ ${compressedFormat} == "gatk" ] || [ ${compressedFormat} == "gatkmark" ] || ([ ${compressedFormat} == "cleanall" ]|| [ ${compressedFormat} == "force" ]||[ ${compressedFormat} == "test" ]||[ ${compressedFormat} == "exomiser" ]||[ ${compressedFormat} == "sv" ] && [ ${compressedFormat} == "mysql" ] && [ ${compressedFormat} == "spliceai" ]); then
   echo "bwa ready"
else
echo "## Run BWA ##" > ${oPath}/${sName}/${sName}_status.txt
${bwaPath} mem -R '@RG\tID:foo\tSM:'${sName}'\tLB:bar\tPL:illumina\tPU:run_std' -t 4 ${shellPath}/reference/hg19/hg19.fasta  ${oPath}/${sName}/${sName}_1_sequence.sheared.txt  ${oPath}/${sName}/${sName}_2_sequence.sheared.txt|${samtoolsPath} view -bS -o ${oPath}/${sName}/${sName}.bam -
fi
# remove sheared files after bwa
rm ${oPath}/${sName}/*sheared.* 

################# 4.preprocessing bam files #########################################################################
mkdir -p ${oPath}/${sName}/BAM
if [ -n "${compressedFormat}" ] && ([ ${compressedFormat} == "recal" ]|| [ ${compressedFormat} == "realign" ]||  [ ${compressedFormat} == "markup" ]||[ ${compressedFormat} == "gatk" ]|| [ ${compressedFormat} == "annov" ]|| [ ${compressedFormat} == "filter" ]|| [ ${compressedFormat} == "clean" ] || [ ${compressedFormat} == "qc" ] ||[ ${compressedFormat} == "gatk" ] || [ ${compressedFormat} == "gatkmark" ] || [ ${compressedFormat} == "cleanall" ]|| [ ${compressedFormat} == "force" ]||[ ${compressedFormat} == "test" ]||[ ${compressedFormat} == "exomiser" ]||[ ${compressedFormat} == "sv" ] && [ ${compressedFormat} == "mysql" ] && [ ${compressedFormat} == "spliceai" ]); then
   echo "bam ready"
else
echo "## Run samtools ##"
${samtoolsPath} sort -n  ${oPath}/${sName}/${sName}.bam  ${oPath}/${sName}/BAM/${sName}.sorted
${samtoolsPath} fixmate ${oPath}/${sName}/BAM/${sName}.sorted.bam  ${oPath}/${sName}/BAM/${sName}.matefixed.bam
${samtoolsPath} sort  ${oPath}/${sName}/BAM/${sName}.matefixed.bam  ${oPath}/${sName}/BAM/${sName}.matefixed.sorted
fi
# remove tmp bam after bam preprocess
rm ${oPath}/${sName}/*.bam
rm ${oPath}/${sName}/BAM/${sName}.sorted.bam
rm ${oPath}/${sName}/BAM/${sName}.matefixed.bam

################# 5. dedupping bam files #########################################################################
if [ -n "${compressedFormat}" ] && ([ ${compressedFormat} == "recal" ]|| [ ${compressedFormat} == "realign" ]||[ ${compressedFormat} == "gatk" ]|| [ ${compressedFormat} == "annov" ]|| [ ${compressedFormat} == "filter" ]|| [ ${compressedFormat} == "clean" ] || [ ${compressedFormat} == "qc" ] ||[ ${compressedFormat} == "gatk" ] || [ ${compressedFormat} == "gatkmark" ] || [ ${compressedFormat} == "cleanall" ]|| [ ${compressedFormat} == "force" ]||[ ${compressedFormat} == "test" ]||[ ${compressedFormat} == "exomiser" ]||[ ${compressedFormat} == "sv" ] && [ ${compressedFormat} == "mysql" ] && [ ${compressedFormat} == "spliceai" ]); then
   echo "markdup ready"
else
echo "## Run samtools ##"
${gatkPath} --java-options "-Xmx8g" MarkDuplicates -I ${oPath}/${sName}/BAM/${sName}.matefixed.sorted.bam -M ${oPath}/${sName}/BAM/${sName}.metrics_file -O ${oPath}/${sName}/BAM/${sName}.rmdup.bam --VALIDATION_STRINGENCY LENIENT --VERBOSITY INFO --COMPRESSION_LEVEL 5 --CREATE_INDEX true --CREATE_MD5_FILE false --ASSUME_SORTED true --REMOVE_DUPLICATES true --TMP_DIR=/${shellPath}/novaseq/tmp
${samtoolsPath} index ${oPath}/${sName}/BAM/${sName}.rmdup.bam
fi
# remove tmp bam after bam preprocess
rm ${oPath}/${sName}/BAM/${sName}.matefixed.sorted.bam

mkdir -p ${oPath}/${sName}/RECALIBRATION
echo "## Realignment ##" > ${oPath}/${sName}/${sName}_status.txt
echo "## Run BaseRecalibrator ##"
if [ -n "${compressedFormat}" ] && ([ ${compressedFormat} == "realign" ]|| [ ${compressedFormat} == "gatk" ]|| [ ${compressedFormat} == "annov" ]|| [ ${compressedFormat} == "filter" ]|| [ ${compressedFormat} == "clean" ] || [ ${compressedFormat} == "gatk" ] || [ ${compressedFormat} == "gatkmark" ] || [ ${compressedFormat} == "qc" ]|| [ ${compressedFormat} == "cleanall" ]|| [ ${compressedFormat} == "force" ]||[ ${compressedFormat} == "test" ]||[ ${compressedFormat} == "exomiser" ]||[ ${compressedFormat} == "sv" ] && [ ${compressedFormat} == "mysql" ] && [ ${compressedFormat} == "spliceai" ]); then
   echo "recal ready"
else
${gatkPath} --java-options "-Xmx8g" BaseRecalibrator -R ${shellPath}/reference/hg19/hg19.fasta --known-sites ${shellPath}/reference/dbsnp/dbsnp_132.hg19.vcf -I ${oPath}/${sName}/BAM/${sName}.rmdup.bam -O ${oPath}/${sName}/RECALIBRATION/${sName}.recal_data.cvs
echo "## Run ApplyBQSR ##"
${gatkPath} --java-options "-Xmx8g" ApplyBQSR -R ${shellPath}/reference/hg19/hg19.fasta -I ${oPath}/${sName}/BAM/${sName}.rmdup.bam -bqsr ${oPath}/${sName}/RECALIBRATION/${sName}.recal_data.cvs -O ${oPath}/${sName}/RECALIBRATION/${sName}.recal.bam
${samtoolsPath} index ${oPath}/${sName}/RECALIBRATION/${sName}.recal.bam
fi
# remove tmp bam after bam preprocess
rm ${oPath}/${sName}/BAM/${sName}.rmdup.bam

################# 6. realign bam files #########################################################################
mkdir -p ${oPath}/${sName}/REALIGNMENT
#Indel realignment is no longer necessary for variant discovery if you plan to use a variant caller that performs a haplotype assembly step, such as HaplotypeCaller or MuTect2.
if [ -n "${compressedFormat}" ] && ([ ${compressedFormat} == "gatk" ]|| [ ${compressedFormat} == "annov" ]|| [ ${compressedFormat} == "filter" ]|| [ ${compressedFormat} == "clean" ] || [ ${compressedFormat} == "qc" ] || [ ${compressedFormat} == "gatk" ] || [ ${compressedFormat} == "gatkmark" ] || [ ${compressedFormat} == "cleanall" ]|| [ ${compressedFormat} == "force" ]||[ ${compressedFormat} == "test" ]||[ ${compressedFormat} == "exomiser" ]||[ ${compressedFormat} == "sv" ] && [ ${compressedFormat} == "mysql" ] && [ ${compressedFormat} == "spliceai" ]); then
   echo "realign ready"
else
${samtoolsPath} sort  ${oPath}/${sName}/RECALIBRATION/${sName}.recal.bam  ${oPath}/${sName}/REALIGNMENT/${sName}.local_realigned.sorted
${samtoolsPath} index ${oPath}/${sName}/REALIGNMENT/${sName}.local_realigned.sorted.bam
# remove tmp bam after bam preprocess
rm ${oPath}/${sName}/RECALIBRATION/${sName}.recal.bam

#generate chrX reads ratio for gender QC
${samtoolsPath} view -o ${oPath}/${sName}/REALIGNMENT/${sName}.X.sam ${oPath}/${sName}/REALIGNMENT/${sName}.local_realigned.sorted.bam chrX
${samtoolsPath} view -o ${oPath}/${sName}/REALIGNMENT/${sName}.sam ${oPath}/${sName}/REALIGNMENT/${sName}.local_realigned.sorted.bam
wc -l ${oPath}/${sName}/REALIGNMENT/${sName}.X.sam > ${oPath}/${sName}/REALIGNMENT/${sName}.X.ratio
wc -l ${oPath}/${sName}/REALIGNMENT/${sName}.sam >> ${oPath}/${sName}/REALIGNMENT/${sName}.X.ratio
rm ${oPath}/${sName}/REALIGNMENT/${sName}.X.sam
rm ${oPath}/${sName}/REALIGNMENT/${sName}.sam
fi

################# 7. GATK snv calling #########################################################################
echo "## Run GATK HaplotypeCaller ##" > ${oPath}/${sName}/${sName}_status.txt
if [ -n "${compressedFormat}" ] && ([ ${compressedFormat} == "annov" ]|| [ ${compressedFormat} == "filter" ]|| [ ${compressedFormat} == "clean" ] || [ ${compressedFormat} == "qc" ] || [ ${compressedFormat} == "gatkmark" ] || [ ${compressedFormat} == "cleanall" ]|| [ ${compressedFormat} == "force" ]||[ ${compressedFormat} == "test" ]||[ ${compressedFormat} == "exomiser" ]||[ ${compressedFormat} == "sv" ] && [ ${compressedFormat} == "mysql" ] && [ ${compressedFormat} == "spliceai" ]); then
   echo "gatk ready"
else
${gatkPath} --java-options "-Xmx8g" HaplotypeCaller -R ${shellPath}/reference/hg19/hg19.fasta -I ${oPath}/${sName}/REALIGNMENT/${sName}.local_realigned.sorted.bam -O ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.g.vcf -ERC GVCF 
#${gatkPath} --java-options "-Xmx8g" HaplotypeCaller -R ${shellPath}/reference/hg19/hg19.fasta -I ${oPath}/${sName}/REALIGNMENT/${sName}.local_realigned.sorted.bam -O ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.g.vcf -ERC GVCF  -bamout ${oPath}/${sName}/REALIGNMENT/${sName}.HaplotypeCaller.bam
${gatkPath} --java-options "-Xmx4g" GenotypeGVCFs -R ${shellPath}/reference/hg19/hg19.fasta --variant ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.g.vcf -O ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.vcf
  if [ $? -eq 0 ]; then
    message=${sName}" GATK Done" 
else
    message=${sName}" GATK FAIL"
  fi
python ${shellPath}/bin/mail/job_mail.py ${message}
fi

################# 8. GATK QC filtering #########################################################################
if [ -n "${compressedFormat}" ] && ([ ${compressedFormat} == "annov" ]|| [ ${compressedFormat} == "filter" ]|| [ ${compressedFormat} == "clean" ] || [ ${compressedFormat} == "qc" ] || [ ${compressedFormat} == "cleanall" ]|| [ ${compressedFormat} == "force" ]||[ ${compressedFormat} == "test" ]||[ ${compressedFormat} == "exomiser" ]||[ ${compressedFormat} == "sv" ] && [ ${compressedFormat} == "mysql" ] && [ ${compressedFormat} == "spliceai" ]); then
   echo "gatk mark ready"
else
${pythonPath} ${shellPath}/bin/GATKScript/gatk_cis_trans.py -i ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.vcf -o ${oPath}/${sName}/${sName}.GATK.phase.txt
${gatkPath} SelectVariants -R ${shellPath}/reference/hg19/hg19.fasta -V ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.vcf -select-type SNP -O ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.snp.vcf
${gatkPath} VariantFiltration -R ${shellPath}/reference/hg19/hg19.fasta -V ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.snp.vcf --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "snp_filter" -O ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.snp.mark.vcf
${gatkPath} SelectVariants -R ${shellPath}/reference/hg19/hg19.fasta -V ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.vcf -select-type INDEL -O ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.indel.vcf
${gatkPath} VariantFiltration -R ${shellPath}/reference/hg19/hg19.fasta -V ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.indel.vcf --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filter-name "indel_filter" -O ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.indel.mark.vcf
${gatkPath} MergeVcfs -I ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.snp.mark.vcf -I ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.indel.mark.vcf -O ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.mark.vcf
fi

################# x. force GATK mutect2 calling #########################################################################
if [ -n "${compressedFormat}" ] && [ ${compressedFormat} == "mutect2" ]; then
${samtoolsPath} view -H ${oPath}/${sName}/REALIGNMENT/${sName}.local_realigned.sorted.bam > ${oPath}/${sName}/REALIGNMENT/${sName}.head
sed 's/$1/${x}/g' ${oPath}/${sName}/REALIGNMENT/${sName}.head > ${oPath}/${sName}/REALIGNMENT/${sName}.header
${samtoolsPath} reheader -P ${oPath}/${sName}/REALIGNMENT/${sName}.header ${oPath}/${sName}/REALIGNMENT/${sName}.local_realigned.sorted.bam > ${oPath}/${sName}/REALIGNMENT/${sName}.local_realigned.sorted.hd.bam
${samtoolsPath} index  ${oPath}/${sName}/REALIGNMENT/${sName}.local_realigned.sorted.hd.bam
${gatkPath} Mutect2 -R ${shellPath}/reference/hg19/hg19.fasta -I  ${oPath}/${sName}/REALIGNMENT/${sName}.local_realigned.sorted.hd.bam -tumor $x -O  ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.mutect2.vcf
${pythonPath} ${shellPath}/bin/GATKScript/gatk_cis_trans.py -i ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.vcf -o ${oPath}/${sName}/${sName}.GATK.phase.txt
${gatkPath} SelectVariants -R ${shellPath}/reference/hg19/hg19.fasta -V ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.mutect2.vcf -select-type SNP -O ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.snp.vcf
${gatkPath} VariantFiltration -R ${shellPath}/reference/hg19/hg19.fasta -V ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.snp.vcf --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "snp_filter" -O ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.snp.mark.vcf
${gatkPath} SelectVariants -R ${shellPath}/reference/hg19/hg19.fasta -V ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.mutect2.vcf -select-type INDEL -O ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.indel.vcf
${gatkPath} VariantFiltration -R ${shellPath}/reference/hg19/hg19.fasta -V ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.indel.vcf --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filter-name "indel_filter" -O ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.indel.mark.vcf
fi
# remove tmp bam after bam preprocess
rm ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.snp.*
rm ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.indel.*
rm ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.g.*

################# 9. annovar annotating #########################################################################
if [ -n "${compressedFormat}" ] && ([ ${compressedFormat} == "filter" ] || [ ${compressedFormat} == "clean" ] || [ ${compressedFormat} == "qc" ] || [ ${compressedFormat} == "cleanall" ]|| [ ${compressedFormat} == "force" ]||[ ${compressedFormat} == "test" ]||[ ${compressedFormat} == "exomiser" ]||[ ${compressedFormat} == "sv" ] && [ ${compressedFormat} == "mysql" ] && [ ${compressedFormat} == "spliceai" ]); then
   echo "annov ready"
else
  mkdir -p ${oPath}/${sName}/SNP_CALL
  ${pythonPath} ${shellPath}/bin/GATKScript/format_GATK_SNP.py -i ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.mark.vcf -o ${oPath}/${sName}/SNP_CALL/${sName}.snp.vcf
  ${pythonPath} ${shellPath}/bin/annotate_filter/pipeline_filter_and_annotate.WGS.py -i ${oPath}/${sName}/SNP_CALL/${sName}.snp.vcf -o ${oPath}/${sName}/SNP_CALL/${sName}.snp.filt+anno.vcf -d ${sName} -t SNP -c ${cutoff} -w ${oPath}/${sName}/SNP_CALL/${sName}.snp.filt+anno.vcf.all 

  mkdir -p ${oPath}/${sName}/INDEL_CALL
  ${pythonPath} ${shellPath}/bin/GATKScript/format_GATK_INDEL.py -i ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.mark.vcf -o ${oPath}/${sName}/INDEL_CALL/${sName}.indel.vcf
  ${pythonPath} ${shellPath}/bin/annotate_filter/pipeline_filter_and_annotate.WGS.py -i ${oPath}/${sName}/INDEL_CALL/${sName}.indel.vcf -o ${oPath}/${sName}/INDEL_CALL/${sName}.indel.filt+anno.vcf -d ${sName} -t INDEL -c ${cutoff} -w ${oPath}/${sName}/INDEL_CALL/${sName}.indel.filt+anno.vcf.all
fi
  
    ${pythonPath} ${shellPath}/bin/annotate_filter/convert_to_readable_form.py -s ${oPath}/${sName}/SNP_CALL/${sName}.snp.filt+anno.vcf -i ${oPath}/${sName}/INDEL_CALL/${sName}.indel.filt+anno.vcf -o ${oPath}/${sName}/${sName}.format.sorted.analysis


    ${pythonPath} ${shellPath}/bin/annotate_filter/convert_to_readable_form.py -s ${oPath}/${sName}/SNP_CALL/${sName}.snp.filt+anno.vcf.all -i ${oPath}/${sName}/INDEL_CALL/${sName}.indel.filt+anno.vcf.all -o ${oPath}/${sName}/${sName}.format.sorted.analysis.all

if [ -n "${ifWES}" ] && [ ${ifWES} == "n" ]; then
    ${pythonPath} ${shellPath}/bin/annotate_filter/convert_to_readable_form.py -s ${oPath}/${sName}/SNP_CALL/${sName}.snp.filt+anno.vcf.all -i ${oPath}/${sName}/INDEL_CALL/${sName}.indel.filt+anno.vcf.all -o ${oPath}/${sName}/${sName}.format.sorted.analysis.all
    cp ${oPath}/${sName}/${sName}.format.sorted.analysis ${oPath}/${sName}/${sName}.format.sorted.analysis.wes
    cp ${oPath}/${sName}/${sName}.format.sorted.analysis.all ${oPath}/${sName}/${sName}.format.sorted.analysis
fi

################# x. force QC checking #########################################################################
if [ -n "${compressedFormat}" ] && [ ${compressedFormat} == "qc" ]; then
   echo "## Run QC ##" > ${oPath}/${sName}/${sName}_status.txt
   sh ${shellPath}/qc.sh  -s ${sName} -d ${designFile} -o ${oPath}
   exit
fi
################# x. force file cleaning #########################################################################
if [ -n "${compressedFormat}" ] && ([ ${compressedFormat} == "clean" ]); then
   sh ${shellPath}/bam.cleanup.sh ${sName} ${oPath}
   exit
fi

################# noncoding annotation ########################################################################
sh ${shellPath}/bin/otherScript/noncoding_annotation.sh ${sName} ${oPath} ${cutoff}

################# 10. filter skipping #########################################################################

if [ -n "${compressedFormat}" ] && ([ ${compressedFormat} == "test" ] && [ ${compressedFormat} == "exomiser" ] && [ ${compressedFormat} == "mysql" ] && [ ${compressedFormat} == "spliceai" ] && [ ${compressedFormat} == "sv" ]); then
   echo "test ready"
else
echo "## Filter and Annotate ##" > ${oPath}/${sName}/${sName}_status.txt

#use gnomad to filter variants(max allele frequency > 0.005)
${pythonPath} ${shellPath}/bin/otherScript/gnomad_wes_wgs_filter.py ${oPath}/${sName}/${sName}.format.sorted.analysis ${oPath}/${sName}/${sName}.origin.analysis ${cutoff}
#${pythonPath} ${shellPath}/bin/otherScript/gnomad_wgs_filter.py ${oPath}/${sName}/${sName}.format.sorted.analysis > ${oPath}/${sName}/${sName}.origin.analysis

#modify hgmd annotation
${pythonPath} ${shellPath}/bin/otherScript/hgmd_patch.py ${oPath}/${sName}/${sName}.origin.analysis > ${oPath}/${sName}/${sName}.hgmd.analysis

#modify Annovar wrong annotation and add DXDB annotation
${pythonPath} ${shellPath}/bin/otherScript/anno_modify.py ${oPath}/${sName}/${sName}.hgmd.analysis > ${oPath}/${sName}/${sName}.anno.modify.analysis

#add ClinVar annotation
${pythonPath} ${shellPath}/bin/otherScript/clinvar_anno.py ${oPath}/${sName}/${sName}.anno.modify.analysis > ${oPath}/${sName}/${sName}.anno.clinvar.analysis
    
${pythonPath} ${shellPath}/bin/otherScript/variant_filter.py ${geneFile} ${shellPath}/reference/variant_filtered.txt ${oPath}/${sName}/${sName}.anno.clinvar.analysis > ${oPath}/${sName}/${sName}.filter.analysis
 
#internal filter(Baylor)
${pythonPath} ${shellPath}/bin/otherScript/wes_internal_filter.py ${oPath}/${sName}/${sName}.filter.analysis > ${oPath}/${sName}/${sName}.inter.filter.analysis

#ATAC binding site annotation
${pythonPath} ${shellPath}/bin/otherScript/ATAC_label.py ${oPath}/${sName}/${sName}.inter.filter.analysis ${oPath}/${sName}/${sName}.analysis.all
#${pythonPath} ${shellPath}/bin/otherScript/ATAC_filter.py ${oPath}/${sName}/${sName}.analysis.all ${oPath}/${sName}/${sName}.analysis.tmp
 

#${pythonPath} ${shellPath}/bin/otherScript/wes_internal_annotate.py ${oPath}/${sName}/${sName}.inter.filter.analysis > ${oPath}/${sName}/${sName}.inter.filter.anno.analysis

#rf percentage annotation
${pythonPath} ${shellPath}/bin/otherScript/rf_anno.py ${oPath}/${sName}/${sName}.analysis.all > ${oPath}/${sName}/${sName}.atac.rf.analysis
   
#ExAC LOF annotation
${pythonPath} ${shellPath}/bin/otherScript/ExAC_annotate.py ${oPath}/${sName}/${sName}.atac.rf.analysis > ${oPath}/${sName}/${sName}.atac.exac.analysis

#GATK phase annotation
${pythonPath} ${shellPath}/bin/GATKScript/GATK_phase_annotate.py ${oPath}/${sName}/${sName}.atac.exac.analysis  ${oPath}/${sName}/${sName}.GATK.phase.txt > ${oPath}/${sName}/${sName}.atac.exac.phase.analysis

cp ${oPath}/${sName}/${sName}.atac.exac.phase.analysis ${oPath}/${sName}/${sName}.unclassified.analysis

#classify analysis file into WES format
${pythonPath} ${shellPath}/bin/otherScript/divide_WGS.py ${oPath}/${sName}/${sName}.unclassified.analysis > ${oPath}/${sName}/${sName}.classified.analysis 
   
#change column order
#${pythonPath} ${shellPath}/bin/otherScript/column_change_cap.py ${oPath}/${sName}/${sName}.classified.analysis > ${oPath}/${sName}/${sName}.analysis.tmp2

${pythonPath} ${shellPath}/bin/annotate_filter/filter_wgs_final.py -i ${oPath}/${sName}/${sName}.classified.analysis -o ${oPath}/${sName}/${sName}.analysis -c 0.005
fi 

# exomiser annotation
if [ -n "${compressedFormat}" ] && ( [ ${compressedFormat} == "mysql" ] && [ ${compressedFormat} == "spliceai" ] && [ ${compressedFormat} == "sv" ]);then
echo "exomiser ready"
else
if [ -s ${oPath}/${sName}/${sName}.GATK.HaplotypeCaller.mark.vcf ]; then
  sh ${shellPath}/exomiser/run_exomiser_pipeline.sh -d $disease -p ${oPath}/${sName} -s ${sName}
  if [ $? -eq 0 ]; then
    message=${sName}" exomiser Done" 
else
    message=${sName}" exomiser FAIL"
  fi
python ${shellPath}/bin/mail/job_mail.py ${message}  
else
  sh ${shellPath}/exomiser/run_exomiser_pipeline.sim.sh -d $disease -p ${oPath}/${sName} -s ${sName}
  if [ $? -eq 0 ]; then
    message=${sName}" exomiser Done" 
else
    message=${sName}" exomiser FAIL"
  fi
python ${shellPath}/bin/mail/job_mail.py ${message}
fi

# MGI/ZFN database annotation
${pythonPath} ${shellPath}/bin/otherScript/dbSNP_MGI_ZFIN_anno.py ${oPath}/${sName}/${sName}.exomiser.analysis > ${oPath}/${sName}/${sName}.atac.exo.mgi.analysis
fi



if [ -n "${compressedFormat}" ] && ([ ${compressedFormat} == "spliceai" ] && [ ${compressedFormat} == "sv" ]);then
echo "mysql dxdb annotation ready"
else
#internal annotation(Clinbytes)
disease1=$(awk -v d="$disease" '($1==d){print $2}' ${shellPath}/disease.table)
if [ -n ${disease1} ]; then
echo ${disease1}
else
disease1="RP"
fi
${pythonPath} ${shellPath}/bin/otherScript/wes_internal_annotate.py ${oPath}/${sName}/${sName}.atac.exo.mgi.analysis  ${disease1} > ${oPath}/${sName}/${sName}.final.snv.analysis
  if [ $? -eq 0 ]; then
    message=${sName}" mysql Done" 
else
    message=${sName}" mysql FAIL"
  fi
python ${shellPath}/bin/mail/job_mail.py ${message}
fi

if [ -n "${compressedFormat}" ] && [ ${compressedFormat} == "sv" ];then
echo "spliceai ready"
else
#revel score annotation [dbSNP3.5 already build-in]
${pythonPath} ${shellPath}/bin/otherScript/revel_anno.py ${oPath}/${sName}/${sName}.exo.mgi.analysis > ${oPath}/${sName}/${sName}.atac.exo.mgi.rev.analysis
sbatch ${shellPath}/Software/spliceai-1.1.1/run.spliceai.sh ${oPath}/${sName}/${sName}.format.sorted.analysis
sbatch ${shellPath}/Software/spliceai-1.1.1/run.spliceai.wgs.sh ${oPath}/${sName}/${sName}.format.sorted.analysis.all.filtered
  if [ $? -eq 0 ]; then
    message=${sName}" spliceai Done" 
else
    message=${sName}" spliceai FAIL"
  fi
python ${shellPath}/Pipeline/pipeline/pipeline_restructure/bin/mail/job_mail.py ${message}
fi

#change column order
${pythonPath} ${shellPath}/bin/otherScript/column_change_wgs.py ${oPath}/${sName}/${sName}.atac.exo.mgi.repeat.phaseCon.gnomadwgs.tfbs.sim.analysis > ${oPath}/${sName}/${sName}.atac.exo.mgi.repeat.phaseCon.gnomadwgs.tfbs.sim.cc.analysis
rm ${oPath}/${sName}/${sName}.atac.exo.mgi.svcnv.analysis

if [ -n "${compressedFormat}" ] && [ ${compressedFormat} == "final" ];then
echo "sv ready"
else
# SVCNV annotation combine
   cd ${oPath}/${sName}/REALIGNMENT
   sh ${shellPath}/svcnv_sim.sh ${sName} ${oPath}/${sName}/REALIGNMENT ${oPath}/${sName}/${sName}.final.snv.analysis ${cutoff} ${wgs} ${svtyper_ok}
${pythonPath} ${shellPath}/bin/otherScript/wes_internal_annotate.py ${oPath}/${sName}/${sName}.atac.exo.mgi.analysis  ${disease1} > ${oPath}/${sName}/${sName}.final.snv.analysis
  if [ $? -eq 0 ]; then
    message=${sName}" SV Done" 
else
    message=${sName}" SV FAIL"
  fi
python ${shellPath}/bin/mail/job_mail.py ${message}
fi
   
cat ${oPath}/${sName}/${sName}.final.snv.analysis ${oPath}/${sName}/REALIGNMENT/{sName}.svtyper.QCDBGVpass.merge.anno.vcf > ${oPath}/${sName}/${sName}.final.analysis


echo "## Run QC ##" > ${oPath}/${sName}/${sName}_status.txt
sh ${shellPath}/qc.sh  -s ${sName} -d ${designFile}  -i ${iPath} -o ${oPath}
${shellPath}/qualimap_v2.2.1/qualimap --java-mem-size=3600m bamqc -bam ${oPath}/${sName}/REALIGNMENT/${sName}.local_realigned.sorted.bam -outdir  ${oPath}/${sName}/REALIGNMENT

echo -e "${iPath}/${sName} endTime:::::::\c" ; date

date_end=$(date +%s)

echo "## Completed(without batch CNV) ##" > ${oPath}/${sName}/${sName}_status.txt
echo "total time: $((date_end-date_start)) s"


