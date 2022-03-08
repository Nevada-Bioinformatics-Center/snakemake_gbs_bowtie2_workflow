import pandas as pd
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.1.2")


##### load config and sample sheets #####

configfile: "config.yaml"
wildcard_constraints:
    sample="\w+"

#units = pd.read_table(config["units"], dtype=str).set_index(["sample", "unit"], drop=False)
units = pd.read_table(config["units"], dtype=str).set_index("sample", drop=False)
#units = pd.read_table(config["units"], dtype=str)
trimmers=config["params"]["trimmers"].split(",")
samples=units["sample"].tolist()
print("Trimmers:", trimmers)
print("Samples:", samples)

def strip_suffix(pattern, suffix):
    return pattern[: -len(suffix)]


##### target rules #####
rule all:
    input:
        "qc/multiqc_report_pretrim.html",
        expand("qc/multiqc_report_posttrim_{trimmer}.html", trimmer=trimmers),
        expand("bowtie2/{trimmer}/all_merged.sorted.rmdup.bam", trimmer=trimmers),
        expand("qc/multiqc_report_bowtie2_{trimmer}.html", trimmer=trimmers),


include: "rules/qc.smk"
include: "rules/trim.smk"
include: "rules/align.smk"
