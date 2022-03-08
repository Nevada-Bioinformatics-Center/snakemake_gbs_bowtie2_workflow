def get_fastq(wildcards):
    return units.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()

def get_fastq1(wildcards):
    fq1 = units.loc[(wildcards.sample), ["fq1"]].dropna()
    return fq1

def get_fastq2(wildcards):
    fq2 = units.loc[(wildcards.sample), ["fq2"]].dropna()
    return fq2



rule fastqc_pretrim_r1:
    input:
       get_fastq1
    output:
        html="qc/fastqc_pretrim/{sample}_r1.html",
        zip="qc/fastqc_pretrim/{sample}_r1_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    log:
        "logs/fastqc_pretrim/{sample}_r1.log"
    resources: time_min=320, mem_mb=8000, cpus=1
    threads: 1
    wrapper:
        "v0.75.0/bio/fastqc"

rule fastqc_pretrim_r2:
    input:
       get_fastq2
    output:
        html="qc/fastqc_pretrim/{sample}_r2.html",
        zip="qc/fastqc_pretrim/{sample}_r2_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    log:
        "logs/fastqc_pretrim/{sample}_r2.log"
    resources: time_min=320, mem_mb=8000, cpus=1
    threads: 1
    wrapper:
        "v0.75.0/bio/fastqc"

rule fastqc_posttrim_r1:
    input:
        "trimmed/{trimmer}/{sample}.1.fastq.gz"
    output:
        html="qc/fastqc_posttrim/{trimmer}/{sample}_r1.html",
        zip="qc/fastqc_posttrim/{trimmer}/{sample}_r1_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    log:
        "logs/fastqc_posttrim/{trimmer}/{sample}_r1.log"
    resources: time_min=320, mem_mb=8000, cpus=1
    threads: 1
    wrapper:
        "v0.75.0/bio/fastqc"

rule fastqc_posttrim_r2:
    input:
        "trimmed/{trimmer}/{sample}.2.fastq.gz"
    output:
        html="qc/fastqc_posttrim/{trimmer}/{sample}_r2.html",
        zip="qc/fastqc_posttrim/{trimmer}/{sample}_r2_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    log:
        "logs/fastqc_posttrim/{trimmer}/{sample}_r2.log"
    resources: time_min=320, mem_mb=8000, cpus=1
    threads: 1
    wrapper:
        "v0.75.0/bio/fastqc"

rule multiqc_pre:
    input:
        expand("qc/fastqc_pretrim/{sample}_r1_fastqc.zip", sample=samples),
        expand("qc/fastqc_pretrim/{sample}_r2_fastqc.zip", sample=samples)
    output:
        "qc/multiqc_report_pretrim.html"
    log:
        "logs/multiqc_pre.log"
    resources: time_min=320, mem_mb=8000, cpus=1
    wrapper:
        "v0.75.0/bio/multiqc"

rule multiqc_post_trimmomatic:
    input:
        expand("logs/trimmomatic/{sample}.log", sample=samples),
        expand("qc/fastqc_posttrim/trimmomatic/{sample}_r1_fastqc.zip", sample=samples),
        expand("qc/fastqc_posttrim/trimmomatic/{sample}_r2_fastqc.zip", sample=samples)
    output:
        "qc/multiqc_report_posttrim_trimmomatic.html"
    log:
        "logs/multiqc_posttrim_trimmomatic.log"
    resources: time_min=320, mem_mb=8000, cpus=1
    wrapper:
        "v0.75.0/bio/multiqc"

rule multiqc_post_fastp:
    input:
        expand("report/fastp/{sample}.fastp.json", sample=samples),
        expand("qc/fastqc_posttrim/fastp/{sample}_r1_fastqc.zip", sample=samples),
        expand("qc/fastqc_posttrim/fastp/{sample}_r2_fastqc.zip", sample=samples)
    output:
        "qc/multiqc_report_posttrim_fastp.html"
    log:
        "logs/multiqc_posttrim_fastp.log"
    resources: time_min=320, mem_mb=8000, cpus=1
    wrapper:
        "v0.75.0/bio/multiqc"

rule multiqc_bowtie2_trimmomatic:
    input:
        expand("bowtie2/trimmomatic/{sample}.sorted.rmdup.bam", sample=samples),
        expand("logs/trimmomatic/{sample}.log", sample=samples),
        expand("logs/bowtie2/trimmomatic/{sample}.log", sample=samples),
        expand("qc/fastqc_posttrim/trimmomatic/{sample}_r1_fastqc.zip", sample=samples),
        expand("qc/fastqc_posttrim/trimmomatic/{sample}_r2_fastqc.zip", sample=samples)
    output:
        "qc/multiqc_report_bowtie2_trimmomatic.html"
    log:
        "logs/multiqc_bowtie2_trimmomatic.log"
    resources: time_min=320, mem_mb=8000, cpus=1
    wrapper:
        "v0.75.0/bio/multiqc"

rule multiqc_bowtie2_fastp:
    input:
        expand("bowtie2/fastp/{sample}.sorted.rmdup.bam", sample=samples),
        expand("logs/bowtie2/fastp/{sample}.log", sample=samples),
        expand("report/fastp/{sample}.fastp.json", sample=samples),
        expand("qc/fastqc_posttrim/fastp/{sample}_r1_fastqc.zip", sample=samples),
        expand("qc/fastqc_posttrim/fastp/{sample}_r2_fastqc.zip", sample=samples)
    output:
        "qc/multiqc_report_bowtie2_fastp.html"
    log:
        "logs/multiqc_bowtie2_fastp.log"
    resources: time_min=320, mem_mb=8000, cpus=1
    wrapper:
        "v0.75.0/bio/multiqc"
