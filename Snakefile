workdir: config["workdir"]

rule all:
    input:
        fastqc_forward = ["results/" + sample + "/fastqc/" + sample + "_1_fastqc.html" for sample in config["samples"]],
        fastqc_reverse = ["results/" + sample + "/fastqc/" + sample + "_2_fastqc.html" for sample in config["samples"]],
        edena = ["results/" + sample + "/edena/" + sample + "_contigs.fasta" for sample in config["samples"]],
        arforward = ["results/" + sample + "/adapterremoval/" + sample + "_1.fastq" for sample in config["samples"]]

# TODO: remember to remove files extracted at the end of the pipeline

rule fastqc:
    input:
        forward = lambda wildcards: os.path.abspath(config["samples"][wildcards.sample]["forward"]),
        reverse = lambda wildcards: os.path.abspath(config["samples"][wildcards.sample]["reverse"])
    params:
        outdir = "results/{sample}/fastqc"
    output:
        forward = "results/{sample}/fastqc/{sample}_1_fastqc.html",
        reverse = "results/{sample}/fastqc/{sample}_2_fastqc.html"
    log:
        stdout = "results/{sample}/fastqc/log-stdout.txt",
        stderr = "results/{sample}/fastqc/log-stderr.txt"
    conda:
        "envs/fastqc.yaml"
    benchmark:
        "results/{sample}/fastqc/benchmark.txt"
    threads:
        config["threads"]
    shell:
        "fastqc --threads {threads} --outdir {params.outdir} {input.forward} {input.reverse} > {log.stdout} 2> {log.stderr}"

rule extract:
    input:
        forward = lambda wildcards: os.path.abspath(config["samples"][wildcards.sample]["forward"]),
        reverse = lambda wildcards: os.path.abspath(config["samples"][wildcards.sample]["reverse"])
    output:
        forward = "results/{sample}/extract-file/{sample}_1.fastq",
        reverse = "results/{sample}/extract-file/{sample}_2.fastq"
    conda:
        "envs/extract-file.yaml"
    benchmark:
        "results/{sample}/extract-file/benchmark.txt"
    shell:
        """
        ./scripts/extract-file.sh {input.forward} {output.forward}
        ./scripts/extract-file.sh {input.reverse} {output.reverse}
        """

rule edena:
    input:
        forward = "results/{sample}/extract-file/{sample}_1.fastq",
        reverse = "results/{sample}/extract-file/{sample}_2.fastq"
    params:
        prefix = "results/{sample}/edena/{sample}"
    output:
        "results/{sample}/edena/{sample}_contigs.fasta"
    log:
        stdout = "results/{sample}/edena/log-stdout.txt",
        stderr = "results/{sample}/edena/log-stderr.txt"
    conda:
        "envs/edena.yaml"
    benchmark:
        "results/{sample}/edena/benchmark.txt"
    threads:
        config["threads"]
    shell:
        """
        edena -nThreads {threads} -paired {input.forward} {input.reverse} -prefix {params.prefix} > {log.stdout} 2> {log.stderr}
        edena -edenaFile {params.prefix}.ovl -prefix {params.prefix} >> {log.stdout} 2>> {log.stderr}
        rm {params.prefix}.ovl
        """

rule adapterremoval:
    input:
        forward = "results/{sample}/extract-file/{sample}_1.fastq",
        reverse = "results/{sample}/extract-file/{sample}_2.fastq"
    params:
        minquality = config['ar-minquality'],
        minlength = config['ar-minlength'],
        optional = config['ar-optional'],
        mm = config['ar-mm'],
        minalignmentlength = config['ar-minalignmentlength'],
        collapsed = "results/{sample}/adapterremoval/{sample}_collapsed.fastq",
        collapsed_truncated = "results/{sample}/adapterremoval/{sample}_collapsed_truncated.fastq"
    output:
        forward = "results/{sample}/adapterremoval/{sample}_1.fastq",
        reverse = "results/{sample}/adapterremoval/{sample}_2.fastq",
        singleton = "results/{sample}/adapterremoval/{sample}_singleton.fastq",
        discarded = "results/{sample}/adapterremoval/{sample}_discarded.fastq",
        settings = "results/{sample}/adapterremoval/{sample}_settings"
    log:
        stdout = "results/{sample}/adapterremoval/log-stdout.txt",
        stderr = "results/{sample}/adapterremoval/log-stderr.txt"
    conda:
        "envs/adapterremoval.yaml"
    benchmark:
        "results/{sample}/adapterremoval/benchmark.txt"
    threads:
        config["threads"]
    shell:
        """
        AdapterRemoval \
--file1 {input.forward} \
--file2 {input.reverse} \
--threads {threads} \
--output1 {output.forward} \
--output2 {output.reverse} \
--singleton {output.singleton} \
--outputcollapsed {params.collapsed} \
--outputcollapsedtruncated {params.collapsed_truncated} \
--discarded {output.discarded} \
{params.optional} \
--minquality {params.minquality} \
--minlength {params.minlength} \
--minalignmentlength {params.minalignmentlength} \
--mm {params.mm} \
--settings {output.settings} \
> {log.stdout} \
2> {log.stderr}
        """
