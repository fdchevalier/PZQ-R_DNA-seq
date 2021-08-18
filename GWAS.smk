from __future__ import print_function
import os
import fnmatch
import pandas as pd

SAMPLES = os.listdir("data/libraries/1-GWAS")
CONTIGS = pd.read_table("data/genome/schistosoma_mansoni.PRJEA36577.WBPS14.genomic.fa.fai", header=None, usecols=[0], squeeze=True, dtype=str)
GENOME = "data/genome/schistosoma_mansoni.PRJEA36577.WBPS14.genomic.fa"
VCF_PFX = "PZQ_GWAS"

rule all:
    input:
        expand("data/libraries/1-GWAS/{sample}/{sample}_sorted.bam", sample=SAMPLES),
        expand("data/libraries/1-GWAS/{sample}/{sample}_sorted.bam.bai", sample=SAMPLES),
        expand("data/libraries/1-GWAS/{sample}/{sample}_sorted_MD.bam", sample=SAMPLES),
        expand("data/libraries/1-GWAS/{sample}/{sample}_sorted_MD.log", sample=SAMPLES),
        expand("data/libraries/1-GWAS/{sample}/{sample}_sorted_MD.bam.bai", sample=SAMPLES),
        expand("data/libraries/1-GWAS/{sample}/{sample}_sorted_MD.grp", sample=SAMPLES),
        expand("data/libraries/1-GWAS/{sample}/{sample}_sorted_MD_recal.bam", sample=SAMPLES),
        expand("data/libraries/1-GWAS/{sample}/{sample}_sorted_MD_recal.bam.bai", sample=SAMPLES),
        expand("data/libraries/1-GWAS/{sample}/{sample}_sorted_MD_recal.flagstat", sample=SAMPLES),
        expand("data/libraries/1-GWAS/{sample}/{sample}.gvcf.gz", sample=SAMPLES),
        expand("data/calling/{vcf_pfx}.gvcf.gz", vcf_pfx=VCF_PFX),
        expand("data/calling/{vcf_pfx}.{contig}.vcf.gz", vcf_pfx=VCF_PFX, contig=CONTIGS),
        expand("data/calling/{vcf_pfx}.vcf.gz", vcf_pfx=VCF_PFX)

rule alignment:
    input:
        read1="data/libraries/1-GWAS/{sample}/{sample}_R1.fastq.gz",
        read2="data/libraries/1-GWAS/{sample}/{sample}_R2.fastq.gz",
        genome=GENOME
    output:
        temp("data/libraries/1-GWAS/{sample}/{sample}_sorted.bam")
    params:
        rg=r"@RG\tID:{sample}\tPL:illumina\tLB:{sample}\tSM:{sample}"
    shell:
        'bwa mem -t $(nproc) -M -R "{params.rg}" "{input.genome}" "{input.read1}" "{input.read2}" | samtools sort -@8 -o "{output}" -'

rule indexing1:
    input:
        "data/libraries/1-GWAS/{sample}/{sample}_sorted.bam"
    output:
        temp("data/libraries/1-GWAS/{sample}/{sample}_sorted.bam.bai")
    shell:
        'samtools index "{input}"'

rule mark_duplicates:
    input:
        bam="data/libraries/1-GWAS/{sample}/{sample}_sorted.bam",
        bai="data/libraries/1-GWAS/{sample}/{sample}_sorted.bam.bai"
    output:
        bam=temp("data/libraries/1-GWAS/{sample}/{sample}_sorted_MD.bam"),
        log=protected("data/libraries/1-GWAS/{sample}/{sample}_sorted_MD.log")
    shell:
        'gatk --java-options "-Xmx2g" MarkDuplicates -I "{input.bam}" -O "{output.bam}" -M "{output.log}" --VALIDATION_STRINGENCY LENIENT --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP $(ulimit -n)'

rule indexing2:
    input:
        "data/libraries/1-GWAS/{sample}/{sample}_sorted_MD.bam"
    output:
        temp("data/libraries/1-GWAS/{sample}/{sample}_sorted_MD.bam.bai")
    shell:
        'samtools index "{input}"'

rule bsqr:
    input:
        bam="data/libraries/1-GWAS/{sample}/{sample}_sorted_MD.bam",
        bai="data/libraries/1-GWAS/{sample}/{sample}_sorted_MD.bam.bai",
        genome=GENOME,
        sites="data/genome/sm_dbSNP_v7.vcf"
    output:
        table=temp("data/libraries/1-GWAS/{sample}/{sample}_sorted_MD.grp"),
        bam=protected("data/libraries/1-GWAS/{sample}/{sample}_sorted_MD_recal.bam")
    params:
        table=r"data/libraries/1-GWAS/{sample}/{sample}_sorted_MD.grp"
    shell:
        'gatk --java-options "-Xmx2g" BaseRecalibrator -R "{input.genome}" -I "{input.bam}" --known-sites "{input.sites}" -O "{output.table}" &&\
         gatk --java-options "-Xmx2g" ApplyBQSR -R "{input.genome}" -I "{input.bam}" --bqsr-recal-file "{params.table}" -O "{output.bam}"'

rule indexing3:
    input:
        "data/libraries/1-GWAS/{sample}/{sample}_sorted_MD_recal.bam"
    output:
        protected("data/libraries/1-GWAS/{sample}/{sample}_sorted_MD_recal.bam.bai")
    shell:
        'samtools index "{input}"'

rule stats:
    input:
        bam="data/libraries/1-GWAS/{sample}/{sample}_sorted_MD_recal.bam",
        bai="data/libraries/1-GWAS/{sample}/{sample}_sorted_MD_recal.bam.bai"
    output:
        protected("data/libraries/1-GWAS/{sample}/{sample}_sorted_MD_recal.flagstat")
    params:
        spl=r"data/libraries/1-GWAS/{sample}/{sample}_sorted_MD_recal.bam"
    shell:
        'echo -e "File: {params.spl}" > "{output}" ; \
         samtools flagstat "{input.bam}" >> "{output}" ; \
         echo  "\n" >> "{output}"'

rule calling:
    input:
        bam="data/libraries/1-GWAS/{sample}/{sample}_sorted_MD_recal.bam",
        bai="data/libraries/1-GWAS/{sample}/{sample}_sorted_MD_recal.bam.bai",
        genome=GENOME,
        sites="data/genome/sm_dbSNP_v7.vcf"
    output:
        "data/libraries/1-GWAS/{sample}/{sample}.gvcf.gz"
    shell:
        'gatk --java-options "-Xmx2g" HaplotypeCaller -R "{input.genome}" -I "{input.bam}" -D "{input.sites}" --output-mode EMIT_ALL_ACTIVE_SITES -ERC GVCF -O "{output}"'

rule combining:
    input:
        gvcfs=expand("data/libraries/1-GWAS/{sample}/{sample}.gvcf.gz", sample=SAMPLES),
        genome=GENOME,
        sites="data/genome/sm_dbSNP_v7.vcf"
    output:
        vcf=temp("data/calling/{vcf_pfx}.gvcf.gz"),
    run:
        gvcfs=" --variant ".join(input.gvcfs)
        shell('gatk --java-options "-Xmx2g" CombineGVCFs -R "{input.genome}" --variant {gvcfs} -D "{input.sites}" -O "{output.vcf}"')

rule indexing_gvcf:
    input:
        gvcf="data/calling/{vcf_pfx}.gvcf.gz"
    output:
        tbi=temp("data/calling/{vcf_pfx}.gvcf.gz.tbi")
    shell:
        'gatk --java-options "-Xmx2g" IndexFeatureFile -I "{input.gvcf}"' 

rule genotype_variants:
    input:
        gvcf="data/calling/{vcf_pfx}.gvcf.gz",
        tbi="data/calling/{vcf_pfx}.gvcf.gz.tbi",
        genome=GENOME,
        sites="data/genome/sm_dbSNP_v7.vcf"
    output:
        vcf=temp("data/calling/{vcf_pfx}.{contig}.vcf.gz"),
        tbi=temp("data/calling/{vcf_pfx}.{contig}.vcf.gz.tbi")
    params:
        contig=r"{contig}"
    shell:
        'gatk --java-options "-Xmx2g"  GenotypeGVCFs -R "{input.genome}" -V "{input.gvcf}" -D "{input.sites}" -L {params.contig} -O "{output.vcf}"'

rule merge_variants:
    input:
        vcfs=expand("data/calling/{vcf_pfx}.{contig}.vcf.gz", vcf_pfx=VCF_PFX, contig=CONTIGS)
    output:
        protected("data/calling/{vcf_pfx}.vcf.gz")
    run:
        vcfs=" -I ".join('"{0}"'.format(w) for w in input.vcfs)
        shell('gatk --java-options "-Xmx2g"  MergeVcfs -I {vcfs} -O "{output}"')
