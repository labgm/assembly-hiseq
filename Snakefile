workdir: config["workdir"]

rule all:
    input:
        fastqc_forward = ["results/" + sample + "/fastqc/" + sample + "_1_fastqc.html" for sample in config["samples"]],
        fastqc_reverse = ["results/" + sample + "/fastqc/" + sample + "_2_fastqc.html" for sample in config["samples"]],
        cdhit = ["results/" + sample + "/cdhit/contigs.fasta" for sample in config["samples"]]

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
    resources:
        config["mem_mb"]
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

rule kmerstream:
    input:
        forward = "results/{sample}/adapterremoval/{sample}_1.fastq",
        reverse = "results/{sample}/adapterremoval/{sample}_2.fastq",
        singleton = "results/{sample}/adapterremoval/{sample}_singleton.fastq"
    params:
        collapsed = "results/{sample}/adapterremoval/{sample}_collapsed.fastq",
        collapsed_truncated = "results/{sample}/adapterremoval/{sample}_collapsed_truncated.fastq"
    output:
        "results/{sample}/kmerstream/ar-{sample}.tsv"
    log:
        stdout = "results/{sample}/kmerstream/log-stdout.txt",
        stderr = "results/{sample}/kmerstream/log-stderr.txt"
    conda:
        "envs/kmerstream.yaml"
    benchmark:
        "results/{sample}/kmerstream/benchmark.txt"
    threads:
        config["threads"]
    shell:
        """
params=()
if [[ -f {params.collapsed} && -f {params.collapsed_truncated} ]]; then
    params+=({params.collapsed} {params.collapsed_truncated})
fi
KmerStream \
--kmer-size=7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,\
51,53,55,57,59,61,63,65,67,69,71,73,75,77,79,81,83,85,87,89,91,93,95,97,99,\
101,103,105,107,109,111,113,115,117,119,121,123,125,127 \
--output={output} \
--threads={threads} \
--tsv \
{input.forward} \
{input.reverse} \
{input.singleton} "${{params[@]}}"
        """

rule spades:
    input:
        kmerstream = "results/{sample}/kmerstream/ar-{sample}.tsv",
        forward = "results/{sample}/adapterremoval/{sample}_1.fastq",
        reverse = "results/{sample}/adapterremoval/{sample}_2.fastq",
        singleton = "results/{sample}/adapterremoval/{sample}_singleton.fastq"
    params:
        prefix = "results/{sample}/spades",
        collapsed = "results/{sample}/adapterremoval/{sample}_collapsed.fastq",
        collapsed_truncated = "results/{sample}/adapterremoval/{sample}_collapsed_truncated.fastq"
    output:
        "results/{sample}/spades/scaffolds.fasta"
    log:
        stdout = "results/{sample}/spades/log-stdout.txt",
        stderr = "results/{sample}/spades/log-stderr.txt"
    conda:
        "envs/spades.yaml"
    benchmark:
        "results/{sample}/spades/benchmark.txt"
    threads:
        config["threads"]
    resources:
        config["mem_mb"]
    shell:
        """
kmers=$(tail -n +2 {input.kmerstream} \
| sort -g -r -k3 \
| cut -f2 \
| head -n 5 \
| tr '\\n' ',' \
| rev \
| cut -c 2- \
| rev)
params=()
if [[ -f {params.collapsed} && -f {params.collapsed_truncated} ]]; then
    params+=(--merged {params.collapsed} --merged {params.collapsed_truncated})
fi
spades.py \
--memory {resources.mem_mb} \
-1 {input.forward} \
-2 {input.reverse} \
-s {input.singleton} "${{params[@]}}" \
--threads {threads} \
-k $kmers \
-o {params.prefix} \
> {log.stdout} \
2> {log.stderr}
rm -rf {params.prefix}/corrected
        """

rule unicycler:
    input:
        forward = "results/{sample}/extract-file/{sample}_1.fastq",
        reverse = "results/{sample}/extract-file/{sample}_2.fastq"
    params:
        prefix = "results/{sample}/unicycler"
    output:
        "results/{sample}/unicycler/assembly.fasta"
    log:
        stdout = "results/{sample}/unicycler/log-stdout.txt",
        stderr = "results/{sample}/unicycler/log-stderr.txt"
    conda:
        "envs/unicycler.yaml"
    benchmark:
        "results/{sample}/unicycler/benchmark.txt"
    threads:
        config["threads"]
    shell:
        """
        unicycler -1 {input.forward} -2 {input.reverse} -o {params.prefix} > {log.stdout} 2> {log.stderr}
        """

rule cdhit:
    input:
        edena = "results/{sample}/edena/{sample}_contigs.fasta",
        spades = "results/{sample}/spades/scaffolds.fasta",
        unicycler = "results/{sample}/unicycler/assembly.fasta"
    params:
        version = config['ch-version'],
        identity = config['ch-identity']
        program = config['ch-program']
        circle = config['ch-circle']
    output:
        "results/{sample}/cdhit/contigs.fasta"
    log:
        stdout = "results/{sample}/cdhit/log-stdout.txt",
        stderr = "results/{sample}/cdhit/log-stderr.txt"
    conda:
        "envs/cdhit.yaml"
    benchmark:
        "results/{sample}/cdhit/benchmark.txt"
    threads:
        config["threads"]
    resources:
        config["mem_mb"]
    run:
        if params.version == 'est':
            shell:
                """
                cat {input.edena} {input.spades} {input.unicycler} > results/{sample}/cdhit/concat.fasta
                cd-hit-est -i results/{sample}/cdhit/concat.fasta -o {output} -T {threads} -mask N -c {params.identity} -M {resources.mem_mb}
                """
        else:
            shell:
                """
                cat {input.edena} {input.spades} {input.unicycler} > results/{sample}/cdhit/concat.fasta
                psi-cd-hit.pl -i results/{sample}/cdhit/concat.fasta -o {output} -c {params.identity} -prog {params.program} -circle {params.circle}
                """
