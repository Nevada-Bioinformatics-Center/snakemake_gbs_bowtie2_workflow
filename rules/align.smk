#def get_trim_fastq1(wildcards):
#    fq1 = expand("trimmed/{trimmer}/{sample}.1.fastq.gz", **wildcards)
#    return fq1
#
#def get_trim_fastq2(wildcards):
#    fq2 = expand("trimmed/{trimmer}/{sample}.2.fastq.gz", **wildcards)
#    return fq2

#def get_trim_fastq1(wildcards):
#    fq1 = "trimmed/{trimmer}/{sample}.1.fastq.gz"
#    return fq1
#
#def get_trim_fastq2(wildcards):
#    fq2 = "trimmed/{trimmer}/{sample}.2.fastq.gz"
#    return fq2

##Bowtie alignment
rule bowtie2_build:
    input:
        reference=config["ref"]["genomefa"]
    output:
        multiext(
            config["ref"]["index"],
            ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2",
        ),
    log:
        "logs/bowtie2_build/build.log"
    params:
        extra=""  # optional parameters
    resources: time_min=480, mem_mb=40000, cpus=16
    threads: 16
    wrapper:
        "0.84.0/bio/bowtie2/build"

rule bowtie2:
    input:
        #sample=[get_trim_fastq1, get_trim_fastq2]
        sample=["trimmed/{trimmer}/{sample}.1.fastq.gz", "trimmed/{trimmer}/{sample}.2.fastq.gz"],
	ref=config["ref"]["index"] + ".1.bt2"
    output:
        "bowtie2/{trimmer}/{sample}.bam"
    log:
        "logs/bowtie2/{trimmer}/{sample}.log"
    params:
        index=config["ref"]["index"],  # prefix of reference genome index (built with bowtie2-build)
        extra="--rg-id {sample} --rg SM:{sample} --rg LB:{sample} --rg PL:Illumina"  # optional parameters
    resources: time_min=480, mem_mb=30000, cpus=32
    threads: 32  # Use at least two threads
    conda:
        "../envs/bowtie2.yaml"
    shell:
        "(bowtie2 --threads {threads} {params.extra} -x {params.index} -1 {input.sample[0]} -2 {input.sample[1]} | samtools view -Sbh -o {output[0]} -) 2> {log}"
    #wrapper:
    #    "0.84.0/bio/bowtie2/build"


rule sambamba_sort:
    input:
        "bowtie2/{trimmer}/{sample}.bam"
    output:
        "bowtie2/{trimmer}/{sample}.sorted.bam"
    log:
        "logs/bowtie2/{trimmer}/sambamba-sort/{sample}.log"
    params: ""
    threads: 8
    resources: time_min=480, mem_mb=30000, cpus=8
    wrapper:
        "0.73.0/bio/sambamba/sort"

rule samtools_index:
    input:
        "bowtie2/{trimmer}/{sample}.sorted.bam"
    output:
        "bowtie2/{trimmer}/{sample}.sorted.bam.bai"
    params:
        "" # optional params string
    resources: time_min=480, mem_mb=2000, cpus=1
    wrapper:
        "0.73.0/bio/samtools/index"

rule sambamba_markdup:
    input:
        "bowtie2/{trimmer}/{sample}.sorted.bam"
    output:
        "bowtie2/{trimmer}/{sample}.sorted.rmdup.bam"
    params:
        extra=""  # optional parameters
    log:
        "logs/sambamba-markdup/{trimmer}/{sample}_rmdup.log"
    threads: 48
    resources: time_min=480, mem_mb=60000, cpus=48
    wrapper:
        "0.84.0/bio/sambamba/markdup"

rule sambamba_merge_fastp:
    input:
        expand("bowtie2/fastp/{sample}.sorted.rmdup.bam", sample=samples)
    output:
        "bowtie2/fastp/all_merged.sorted.rmdup.bam"
    params:
        extra=""  # optional parameters
    log:
        "logs/sambamba-merge/fastp/all_merged.log"
    threads: 48
    resources: time_min=480, mem_mb=60000, cpus=48
    wrapper:
        "0.84.0/bio/sambamba/merge"


rule sambamba_merge_trimmomatic:
    input:
        expand("bowtie2/trimmomatic/{sample}.sorted.rmdup.bam", sample=samples)
    output:
        "bowtie2/trimmomatic/all_merged.sorted.rmdup.bam"
    params:
        extra=""  # optional parameters
    log:
        "logs/sambamba-merge/trimmomatic/all_merged.log"
    threads: 48
    resources: time_min=480, mem_mb=20000, cpus=48
    wrapper:
        "0.84.0/bio/sambamba/merge"

