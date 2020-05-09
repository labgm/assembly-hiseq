workdir: config["workdir"]

rule all:
    input:
        fastqc = ["results/" + sample + "/fastqc/*.html" for sample in config["samples"]]

rule fastqc:
    input:
        forward = lambda wildcards: os.path.abspath(config["samples"][wildcards.sample]["forward"]),
        reverse = lambda wildcards: os.path.abspath(config["samples"][wildcards.sample]["reverse"])
    params:
        output = "results/{sample}/fastqc"
    output:
        "results/{sample}/fastqc/*.html"
    conda:
        "envs/fastqc.yaml"
    benchmark:
        "results/{sample}/fastqc/benchmark.txt"
    threads:
        config["threads"]
    shell:
        "fastqc --threads {threads} --outdir {params.output} {input.forward} {input.reverse}"

