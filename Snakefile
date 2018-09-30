import os
from os.path import join, basename, dirname
import glob
######WES############

configfile: 'config.yml'


#prepare
#global variables
outputdir = config["output_dir"]
sample_info = config["sample_info"]
WORKDIR = join(outputdir,"Workspace")
FASTQDIR = WORKDIR + "/01.Data"
QCDIR = WORKDIR + "/02.Qc"
BAMDIR = WORKDIR + "/03.Alin"
BAMQCDIR = WORKDIR + "/04.Bamqc"
CALLDIR = WORKDIR + "/05.Call"
VCFDIR = WORKDIR + "/06.Result"

#Read in sample information
sample_dict = json.load(open(sample_info,"r"))
contigs = config["contigs"]

os.system("mkdir -p %s %s %s %s %s %s" % (FASTQDIR,QCDIR,BAMDIR,BAMQCDIR,CALLDIR,VCFDIR))

rule all:
    input:
        expand(join(QCDIR,"{sample}.R1_fastqc.zip"),sample=sorted(sample_dict.keys())),
        expand(join(QCDIR,"{sample}.R1_fastqc.html"),sample=sorted(sample_dict.keys())),
        expand(join(QCDIR,"{sample}.R2_fastqc.zip"),sample=sorted(sample_dict.keys())),
        expand(join(QCDIR,"{sample}.R1_fastqc.html"),sample=sorted(sample_dict.keys())),
        expand(join(BAMQCDIR,"{sample}.qc.tsv.gz"),sample=sorted(sample_dict.keys())),
        expand(join(VCFDIR,"{sample}.snp.filter.vcf.hg19_multianno.txt"),sample=sorted(sample_dict.keys())),
        expand(join(VCFDIR,"{sample}.indel.filter.vcf.hg19_multianno.txt"),sample=sorted(sample_dict.keys()))

rule findfastq:
    input:
        r1 = lambda wildcards:sample_dict[wildcards.sample][0],
        r2 = lambda wildcards:sample_dict[wildcards.sample][1]
    output:
        f1 = join(FASTQDIR,"{sample}.R1.fq.gz"),
        f2 = join(FASTQDIR,"{sample}.R2.fq.gz")
    shell:
        "ln -s {input.r1} {output.f1};ln -s {input.r2} {output.f2}"

#Quality control for fastq
rule fastqc:
    input:
        f1 = join(FASTQDIR,"{sample}.R1.fq.gz"),
        f2 = join(FASTQDIR,"{sample}.R2.fq.gz")
    output:
        join(QCDIR,"{sample}.R1_fastqc.zip"),
        join(QCDIR,"{sample}.R1_fastqc.html"),
        join(QCDIR,"{sample}.R2_fastqc.zip"),
        join(QCDIR,"{sample}.R2_fastqc.html")
    shell:
        'software/FastQC/fastqc -f fastq -o %s {input.f1} {input.f2}'% (QCDIR)

#Mapping 
rule map:
    input:
        f1 = join(FASTQDIR,"{sample}.R1.fq.gz"),
        f2 = join(FASTQDIR,"{sample}.R2.fq.gz")
    output:
        join(BAMDIR,"{sample}.raw.bam")
    params:
        prefix="\"@RG\\tID:{sample}\\tLB:{sample}\\tPL:ILLUMINA\\tSM:{sample}\\tPU:ILLUMINA\""
    shell:
        "bwa/bwa-0.7.13/bwa/bwa mem -t 8 %s {input.f1} {input.f2} -R {params.prefix}| /share/software/VariantCalling/samtools/samtools-1.3/samtools view -bS - -o {output}" % config["ref_hg19"]


#-R \"\@RG\\tID:{sample}\\tLB:{sample}\\tPL:ILLUMINA\\tSM:{sample}\\tPU:ILLUMINA\"

#sort 
rule sortbam:
    input:
        join(BAMDIR,"{sample}.raw.bam")
    output:
        join(BAMDIR,"{sample}.sort.bam")
    shell:
        "samtools-1.3/samtools sort -@ 10 -m 2G {input} -o {output}"



#rmdup
rule rmdup:
    input:
        join(BAMDIR,"{sample}.sort.bam")
    output:
        bam = join(BAMDIR,"{sample}.rmdup.bam"),
        index = join(BAMDIR,"{sample}.rmdup.bam.bai"),
        metrics = join(BAMDIR,"{sample}.rmdup.bam.metrics")
    shell:
        "java -Xms5g -Xmx5g -Djava.io.tmpdir=%s -jar picard-tools-1.140/picard.jar MarkDuplicates I={input} O={output.bam} M={output.metrics};/share/software/VariantCalling/samtools/samtools-1.3/samtools index {output.bam}" % BAMDIR
#Quality control for BAM
rule bamqc:
    input:
        join(BAMDIR,"{sample}.rmdup.bam")
    output:
        join(BAMQCDIR,"{sample}.qc.tsv.gz")
    shell:
        "software/bamqc/alfred/alfred_v0.1.9_linux_x86_64bit qc -r %s -b %s -o {output} {input}" % (config["ref_hg19"],config["first_bed"])


#calling using GATK
#
#

rule gatk_hc:
    input:
        bam = join(BAMDIR,"{sample}.rmdup.bam"),
        bed = join(config["beddir"],"target.{contig}.bed")
    output:
        join(CALLDIR,"{sample}.raw.{contig}.vcf")
    shell:
        "java -Xms5g -Xmx5g -Djava.io.tmpdir=%s -jar /GATK/v3.8/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T HaplotypeCaller -nct 8 -I {input.bam} -o {output} -L {input.bed} -R %s -G StandardAnnotation -G StandardHCAnnotation -G RankSumTest -G RMSAnnotation -A QualByDepth -A FisherStrand -A AlleleBalance -A Coverage -A MappingQualityZero -A TandemRepeatAnnotator -A VariantType -A DepthPerAlleleBySample -stand_call_conf 30.0 -rf BadCigar --dontUseSoftClippedBases" % (CALLDIR,config["ref_hg19"])

rule merge_vcf:
    input:
        expand(join(CALLDIR,"{{sample}}.raw.{contig}.vcf"),contig=contigs)
    output:
        join(CALLDIR,"{sample}.raw.vcf")
    run:
        import os
        seq = ""
        for one in input:
            seq = seq + " -V "+ one
       
        os.system("java -cp /GATK/v3.8/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants -R %s %s -out %s" % (config["ref_hg19"],seq,output))


rule extract_snp:
    input:
        join(CALLDIR,"{sample}.raw.vcf")
    output:
        join(CALLDIR,"{sample}.snp.vcf")
    shell:
        "java -Xms5g -Xmx5g -Djava.io.tmpdir=%s -jar GATK/v3.8/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T SelectVariants -R %s -V {input} --selectTypeToExclude INDEL -o {output}" % (CALLDIR,config["ref_hg19"])

rule extract_indel:
    input:
        join(CALLDIR,"{sample}.raw.vcf")
    output:
        join(CALLDIR,"{sample}.indel.vcf")
    shell:
        "java -Xms5g -Xmx5g -Djava.io.tmpdir=%s -jar GATK/v3.8/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T SelectVariants -R %s -V {input} --selectTypeToInclude INDEL -o {output}" % (CALLDIR,config["ref_hg19"])

rule filter_snp:
    input:
        join(CALLDIR,"{sample}.snp.vcf")
    output:
        join(CALLDIR,"{sample}.snp.filter.vcf")
    params:
        snp_f = '--filterExpression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0" --filterName "Standard" '
    shell:
        "java -Xms5g -Xmx5g -Djava.io.tmpdir=%s -jar GATK/v3.8/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T VariantFiltration -R %s -V {input} -o {output} {params.snp_f}" % (CALLDIR,config["ref_hg19"])


rule filter_indel:
    input:
        join(CALLDIR,"{sample}.indel.vcf")
    output:
        join(CALLDIR,"{sample}.indel.filter.vcf")
    params:
        indel_f = ' --filterExpression "QD < 2.0 || ReadPosRankSum < -20.0 || FS > 200.0 || SOR > 10.0" --filterName "Standard" '
    shell:
        "java -Xms5g -Xmx5g -Djava.io.tmpdir=%s -jar GATK/v3.8/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T VariantFiltration -R %s -V {input} -o {output} {params.indel_f}" % (CALLDIR,config["ref_hg19"])

rule annotater_snp:
    input:
        join(CALLDIR,"{sample}.snp.filter.vcf")
    output:
        join(VCFDIR,"{sample}.snp.filter.vcf.hg19_multianno.txt"),
        join(VCFDIR,"{sample}.snp.filter.vcf.hg19_multianno.vcf")
    params:
        prefix = join(VCFDIR,"{sample}.snp.filter.vcf")
    shell:
        "perl /share/software/Annotation/annovar/table_annovar.pl  --buildver hg19 --otherinfo --nastring . {input} /share/software/Annotation/annovar/humandb/ -protocol %s -operation g,f,f,f,f,f,f,f,f,f --vcfinput --remove --outfile {params.prefix}" % config["annovar_database"]

rule annotater_indel:
   input:
       join(CALLDIR,"{sample}.indel.filter.vcf")
   output:
       join(VCFDIR,"{sample}.indel.filter.vcf.hg19_multianno.txt"),
       join(VCFDIR,"{sample}.indel.filter.vcf.hg19_multianno.vcf")
   params:
       prefix = join(VCFDIR,"{sample}.indel.filter.vcf")
   shell:
       "perl /share/software/Annotation/annovar/table_annovar.pl  --buildver hg19 --otherinfo --nastring . {input} /share/software/Annotation/annovar/humandb/ -protocol %s -operation g,f,f,f,f,f,f,f,f,f --vcfinput --remove --outfile {params.prefix}" % config["annovar_database"]


