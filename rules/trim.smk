rule trimmomatic_pe:
    input:
        r1=get_fastq1,
        r2=get_fastq2
    output:
        r1="trimmed/trimmomatic/{sample}.1.fastq.gz",
        r2="trimmed/trimmomatic/{sample}.2.fastq.gz",
        # reads where trimming entirely removed the mate
        r1_unpaired="trimmed/trimmomatic/{sample}.1.unpaired.fastq.gz",
        r2_unpaired="trimmed/trimmomatic/{sample}.2.unpaired.fastq.gz"
    log:
        "logs/trimmomatic/{sample}.log"
    params:
        # list of trimmers (see manual)
        trimmer = [f"ILLUMINACLIP:{config['ref']['adapter']}:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36"],
        #trimmer = [f"ILLUMINACLIP:{config['ref']['adapter']}:2:30:10 SLIDINGWINDOW:4:15 LEADING:30 MINLEN:36"],
        # optional parameters
        extra="",
        compression_level="-9"
    threads: 16
    resources: time_min=480, mem_mb=20000, cpus=16
    wrapper:
        "0.73.0/bio/trimmomatic/pe"

rule fastp_pe:
    input:
        sample=get_fastq
    output:
        trimmed=["trimmed/fastp/{sample}.1.fastq.gz", "trimmed/fastp/{sample}.2.fastq.gz"],
        html="report/fastp/{sample}.html",
        json="report/fastp/{sample}.fastp.json"
    log:
        "logs/fastp/{sample}.log"
    params:
        #adapters="--adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        adapters="--detect_adapter_for_pe",
        extra=""
    threads: 16
    resources: time_min=480, mem_mb=20000, cpus=16
    wrapper:
        "0.73.0/bio/fastp"
