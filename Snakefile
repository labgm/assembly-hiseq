workdir: config["workdir"]

rule all:
    input:
        forward = ["results/" + sample + "/fastqc/" + sample + "_1_fastqc.html" for sample in config["samples"]],
        reverse = ["results/" + sample + "/fastqc/" + sample + "_2_fastqc.html" for sample in config["samples"]]

rule fastqc:
    input:
        forward = lambda wildcards: os.path.abspath(config["samples"][wildcards.sample]["forward"]),
        reverse = lambda wildcards: os.path.abspath(config["samples"][wildcards.sample]["reverse"])
    params:
        output = "results/{sample}/fastqc"
    output:
        forward = "results/{sample}/fastqc/{sample}_1_fastqc.html",
        reverse = "results/{sample}/fastqc/{sample}_2_fastqc.html"
    conda:
        "envs/fastqc.yaml"
    benchmark:
        "results/{sample}/fastqc/benchmark.txt"
    threads:
        config["threads"]
    shell:
        "fastqc --threads {threads} --outdir {params.output} {input.forward} {input.reverse}"

