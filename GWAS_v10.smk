from __future__ import print_function
import os
import fnmatch
import pandas as pd

SAMPLES = os.listdir("data/libraries/1-GWAS")
CONTIGS = pd.read_table("data/genome/schistosoma_mansoni.PRJEA36577.WBPS18.genomic.fa.fai", header=None, usecols=[0], squeeze=True, dtype=str)
GENOME = "data/genome/schistosoma_mansoni.PRJEA36577.WBPS18.genomic.fa"
SITES = "data/genome/dbSNP_SM_V10.vcf"
VCF_PFX = "PZQ_GWAS_v10"

rule all:
    input:
        expand("data/libraries/1-GWAS/{sample}/{sample}_v10_sorted.bam", sample=SAMPLES),
        expand("data/libraries/1-GWAS/{sample}/{sample}_v10_sorted.bam.bai", sample=SAMPLES),
        expand("data/libraries/1-GWAS/{sample}/{sample}_v10_sorted_MD.bam", sample=SAMPLES),
        expand("data/libraries/1-GWAS/{sample}/{sample}_v10_sorted_MD.log", sample=SAMPLES),
        expand("data/libraries/1-GWAS/{sample}/{sample}_v10_sorted_MD.bam.bai", sample=SAMPLES),
        expand("data/libraries/1-GWAS/{sample}/{sample}_v10_sorted_MD.grp", sample=SAMPLES),
        expand("data/libraries/1-GWAS/{sample}/{sample}_v10_sorted_MD_recal.bam", sample=SAMPLES),
        expand("data/libraries/1-GWAS/{sample}/{sample}_v10_sorted_MD_recal.bam.bai", sample=SAMPLES),
        expand("data/libraries/1-GWAS/{sample}/{sample}_v10_sorted_MD_recal.flagstat", sample=SAMPLES),
        expand("data/libraries/1-GWAS/{sample}/{sample}_v10.gvcf.gz", sample=SAMPLES),
        expand("data/calling/{vcf_pfx}.gvcf.gz", vcf_pfx=VCF_PFX),
        expand("data/calling/{vcf_pfx}.{contig}.vcf.gz", vcf_pfx=VCF_PFX, contig=CONTIGS),
        expand("data/calling/{vcf_pfx}.vcf.gz", vcf_pfx=VCF_PFX)

rule alignment:
    input:
        read1="data/libraries/1-GWAS/{sample}/{sample}_R1.fastq.gz",
        read2="data/libraries/1-GWAS/{sample}/{sample}_R2.fastq.gz",
        genome=GENOME
    output:
        temp("data/libraries/1-GWAS/{sample}/{sample}_v10_sorted.bam")
    params:
        rg=r"@RG\tID:{sample}\tPL:illumina\tLB:{sample}\tSM:{sample}"
    resources:
        hmem=5 * 10 ** 9,
        cores=48
    shell:
        'bwa mem -t {resources.cores} -M -R "{params.rg}" "{input.genome}" "{input.read1}" "{input.read2}" | samtools sort -@{resources.cores} -o "{output}" -'

rule indexing1:
    input:
        "data/libraries/1-GWAS/{sample}/{sample}_v10_sorted.bam"
    output:
        temp("data/libraries/1-GWAS/{sample}/{sample}_v10_sorted.bam.bai")
    shell:
        'samtools index "{input}"'

rule mark_duplicates:
    input:
        bam="data/libraries/1-GWAS/{sample}/{sample}_v10_sorted.bam",
        bai="data/libraries/1-GWAS/{sample}/{sample}_v10_sorted.bam.bai"
    output:
        bam=temp("data/libraries/1-GWAS/{sample}/{sample}_v10_sorted_MD.bam"),
        log=protected("data/libraries/1-GWAS/{sample}/{sample}_v10_sorted_MD.log")
    resources:
        hmem=10 * 10 ** 9
    shell:
        'gatk --java-options "-Xmx{resources.hmem}" MarkDuplicates -I "{input.bam}" -O "{output.bam}" -M "{output.log}" --VALIDATION_STRINGENCY LENIENT --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP $(ulimit -n)'

rule indexing2:
    input:
        "data/libraries/1-GWAS/{sample}/{sample}_v10_sorted_MD.bam"
    output:
        temp("data/libraries/1-GWAS/{sample}/{sample}_v10_sorted_MD.bam.bai")
    shell:
        'samtools index "{input}"'

rule bsqr:
    input:
        bam="data/libraries/1-GWAS/{sample}/{sample}_v10_sorted_MD.bam",
        bai="data/libraries/1-GWAS/{sample}/{sample}_v10_sorted_MD.bam.bai",
        genome=GENOME,
        sites=SITES
    output:
        table=temp("data/libraries/1-GWAS/{sample}/{sample}_v10_sorted_MD.grp"),
        bam=protected("data/libraries/1-GWAS/{sample}/{sample}_v10_sorted_MD_recal.bam")
    params:
        table=r"data/libraries/1-GWAS/{sample}/{sample}_v10_sorted_MD.grp"
    shell:
        'gatk --java-options "-Xmx{resources.hmem}" BaseRecalibrator -R "{input.genome}" -I "{input.bam}" --known-sites "{input.sites}" -O "{output.table}" &&\
         gatk --java-options "-Xmx{resources.hmem}" ApplyBQSR -R "{input.genome}" -I "{input.bam}" --bqsr-recal-file "{params.table}" -O "{output.bam}"'

rule indexing3:
    input:
        "data/libraries/1-GWAS/{sample}/{sample}_v10_sorted_MD_recal.bam"
    output:
        protected("data/libraries/1-GWAS/{sample}/{sample}_v10_sorted_MD_recal.bam.bai")
    shell:
        'samtools index "{input}"'

rule stats:
    input:
        bam="data/libraries/1-GWAS/{sample}/{sample}_v10_sorted_MD_recal.bam",
        bai="data/libraries/1-GWAS/{sample}/{sample}_v10_sorted_MD_recal.bam.bai"
    output:
        protected("data/libraries/1-GWAS/{sample}/{sample}_v10_sorted_MD_recal.flagstat")
    params:
        spl=r"data/libraries/1-GWAS/{sample}/{sample}_v10_sorted_MD_recal.bam"
    resources:
        hmem=10 * 10 ** 9,
        cores=8
    shell:
        'echo -e "File: {params.spl}" > "{output}" ; \
         samtools flagstat -@{resources.cores} "{input.bam}" >> "{output}" ; \
         echo  "\n" >> "{output}"'

rule calling:
    input:
        bam="data/libraries/1-GWAS/{sample}/{sample}_v10_sorted_MD_recal.bam",
        bai="data/libraries/1-GWAS/{sample}/{sample}_v10_sorted_MD_recal.bam.bai",
        genome=GENOME,
        sites=SITES
    output:
        gvcf="data/libraries/1-GWAS/{sample}/{sample}_v10.gvcf.gz"
    resources:
        hmem=10 * 10 ** 9
    shell:
        'gatk --java-options "-Xmx{resources.hmem}" HaplotypeCaller -R "{input.genome}" -I "{input.bam}" -D "{input.sites}" --output-mode EMIT_ALL_ACTIVE_SITES -ERC GVCF -O "{output.gvcf}"'

rule combining:
    input:
        gvcfs=expand("data/libraries/1-GWAS/{sample}/{sample}_v10.gvcf.gz", sample=SAMPLES),
        genome=GENOME,
        sites=SITES
    output:
        vcf=temp("data/calling/{vcf_pfx}.gvcf.gz"),
    resources:
        hmem=10 * 10 ** 9
    run:
        gvcfs=" --variant ".join(input.gvcfs)
        shell('gatk --java-options "-Xmx{resources.hmem}" CombineGVCFs -R "{input.genome}" -V {gvcfs} -D "{input.sites}" -O "{output}"')

rule indexing_gvcf:
    input:
        gvcf="data/calling/{vcf_pfx}.gvcf.gz"
    output:
        tbi=temp("data/calling/{vcf_pfx}.gvcf.gz.tbi")
    resources:
        hmem=2 * 10 ** 9
    shell:
        'gatk --java-options "-Xmx{resources.hmem}" IndexFeatureFile -I "{input.gvcf}"'

rule genotype_variants:
    input:
        gvcf="data/calling/{vcf_pfx}.gvcf.gz",
        tbi="data/calling/{vcf_pfx}.gvcf.gz.tbi",
        genome=GENOME,
        sites=SITES
    output:
        vcf=temp("data/calling/{vcf_pfx}.{contig}.vcf.gz"),
        tbi=temp("data/calling/{vcf_pfx}.{contig}.vcf.gz.tbi")
    params:
        contig=r"{contig}"
    shell:
        'gatk --java-options "-Xmx{resources.hmem}"  GenotypeGVCFs -R "{input.genome}" -V "{input.gvcf}" -D "{input.sites}" -L {params.contig} -O "{output.vcf}"'

rule merge_variants:
    input:
        vcfs=expand("data/calling/{vcf_pfx}.{contig}.vcf.gz", vcf_pfx=VCF_PFX, contig=CONTIGS)
    output:
        protected("data/calling/{vcf_pfx}.vcf.gz")
    resources:
        hmem=200 * 10 ** 9
    run:
        vcfs=" -I ".join('"{0}"'.format(w) for w in input.vcfs)
        shell('gatk --java-options "-Xmx{resources.hmem}"  MergeVcfs -I {vcfs} -O "{output}"')
