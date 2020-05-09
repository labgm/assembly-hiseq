workdir: config["workdir"]

rule all:
    input:
        fastqc_forward = ["results/" + sample + "/fastqc/" + sample + "_1_fastqc.html" for sample in config["samples"]],
        fastqc_reverse = ["results/" + sample + "/fastqc/" + sample + "_2_fastqc.html" for sample in config["samples"]],
        edena = ["results/" + sample + "/edena/" + sample + "_contigs.fasta" for sample in config["samples"]]

# TODO: remember to remove files extracted by zcat at the end of the pipeline

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

rule zcat:
    input:
        forward = lambda wildcards: os.path.abspath(config["samples"][wildcards.sample]["forward"]),
        reverse = lambda wildcards: os.path.abspath(config["samples"][wildcards.sample]["reverse"])
    output:
        forward = "results/{sample}/zcat/{sample}_1_.fastq",
        reverse = "results/{sample}/zcat/{sample}_2_fastq"
    benchmark:
        "results/{sample}/zcat/benchmark.txt"
    shell:
        """
        zcat {input.forward} > {output.forward}
        zcat {input.reverse} > {output.reverse}
        """

rule edena:
    input:
        forward = "results/{sample}/zcat/{sample}_1_.fastq",
        reverse = "results/{sample}/zcat/{sample}_2_fastq"
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


