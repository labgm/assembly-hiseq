configfile: "config.yaml"

workdir: config["workdir"]

rule all:
    input:
        fastqc_forward = ["results/" + sample + "/fastqc/" + sample + "_1_fastqc.html" for sample in config["samples"]],
        fastqc_reverseR = ["results/" + sample + "/fastqc/" + sample + "_2_fastqc.html" for sample in config["samples"]],
        filter = ["results/" + sample + "/mob_recon/chromosome_filtered.fasta" for sample in config["samples"]],
        prokka = ["results/" + sample + "/prokka/" + sample + ".gbk" for sample in config["samples"]],
        quast = ["results/" + sample + "/quast/report.tsv" for sample in config["samples"]]

# TODO: remember to remove files extracted at the end of the pipeline

rule fastqc:
    input:
        forward = lambda wildcards: os.path.abspath(config["samples"][wildcards.sample]["forward"]),
        reverseR = lambda wildcards: os.path.abspath(config["samples"][wildcards.sample]["reverseR"])
    params:
        outdir = "results/{sample}/fastqc"
    output:
        forward = "results/{sample}/fastqc/{sample}_1_fastqc.html",
        reverseR = "results/{sample}/fastqc/{sample}_2_fastqc.html"
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
        "fastqc --threads {threads} --outdir {params.outdir} {input.forward} {input.reverseR} > {log.stdout} 2> {log.stderr}"

rule extract:
    priority: 10
    input:
        forward = lambda wildcards: os.path.abspath(config["samples"][wildcards.sample]["forward"]),
        reverseR = lambda wildcards: os.path.abspath(config["samples"][wildcards.sample]["reverseR"])
    output:
        forward = "results/{sample}/extract-file/{sample}_1.fastq",
        reverseR = "results/{sample}/extract-file/{sample}_2.fastq"
    conda:
        "envs/extract-file.yaml"
    benchmark:
        "results/{sample}/extract-file/benchmark.txt"
    shell:
        """
        ./scripts/extract-file.sh {input.forward} {output.forward}
        ./scripts/extract-file.sh {input.reverseR} {output.reverseR}
        """

rule adapterremoval:
    priority: 9
    input:
        forward = "results/{sample}/extract-file/{sample}_1.fastq",
        reverseR = "results/{sample}/extract-file/{sample}_2.fastq"
    params:
        minquality = config['adapterremoval']['minquality'],
        minlength = config['adapterremoval']['minlength'],
        optional = config['adapterremoval']['optional'],
        mm = config['adapterremoval']['mm'],
        minalignmentlength = config['adapterremoval']['minalignmentlength'],
        collapsed = "results/{sample}/adapterremoval/{sample}_collapsed.fastq",
        collapsed_truncated = "results/{sample}/adapterremoval/{sample}_collapsed_truncated.fastq"
    output:
        forward = "results/{sample}/adapterremoval/{sample}_1.fastq",
        reverseR = "results/{sample}/adapterremoval/{sample}_2.fastq",
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
        AdapterRemoval --file1 {input.forward} --file2 {input.reverseR} --threads {threads} --output1 {output.forward} --output2 {output.reverseR} --singleton {output.singleton} --outputcollapsed {params.collapsed} --outputcollapsedtruncated {params.collapsed_truncated} --discarded {output.discarded} {params.optional} --minquality {params.minquality} --minlength {params.minlength} --minalignmentlength {params.minalignmentlength} --mm {params.mm} --settings {output.settings} > {log.stdout} 2> {log.stderr}
        """

rule edena:
    input:
        forward = "results/{sample}/extract-file/{sample}_1.fastq",
        reverseR = "results/{sample}/extract-file/{sample}_2.fastq"
    params:
        prefix = "results/{sample}/edena/{sample}",
        edena_r1 = "results/{sample}/adapterremoval/{sample}_1_edena.fastq",
        edena_r2 = "results/{sample}/adapterremoval/{sample}_2_edena.fastq"
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
        mem_mb = config["mem_mb"]
    shell:
        """
        ./scripts/filter_reads.py 150 {input.forward} {params.edena_r1}  > {log.stdout} 2> {log.stderr}
        ./scripts/filter_reads.py 150 {input.reverseR} {params.edena_r2}  > {log.stdout} 2> {log.stderr}
        edena -nThreads {threads} -paired {params.edena_r1} {params.edena_r2} -prefix {params.prefix} > {log.stdout} 2> {log.stderr}
        edena -edenaFile {params.prefix}.ovl -prefix {params.prefix} >> {log.stdout} 2>> {log.stderr}
        rm {params.prefix}.ovl
        """

rule kmerstream:
    input:
        forward = "results/{sample}/adapterremoval/{sample}_1.fastq",
        reverseR = "results/{sample}/adapterremoval/{sample}_2.fastq",
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
        KmerStream --kmer-size=7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57,59,61,63,65,67,69,71,73,75,77,79,81,83,85,87,89,91,93,95,97,99,101,103,105,107,109,111,113,115,117,119,121,123,125,127 --output={output} --threads={threads} --tsv {input.forward} {input.reverseR} {input.singleton} "${{params[@]}}"
        """

rule spades:
    input:
        kmerstream = "results/{sample}/kmerstream/ar-{sample}.tsv",
        forward = "results/{sample}/adapterremoval/{sample}_1.fastq",
        reverseR = "results/{sample}/adapterremoval/{sample}_2.fastq",
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
        mem_gb = 100
    shell:
        """
        kmers=$(tail -n +2 {input.kmerstream} | sort -g -r -k3 | cut -f2 | head -n 5 | tr '\\n' ',' | rev | cut -c 2- | rev)
        params=()
        if [[ -f {params.collapsed} && -f {params.collapsed_truncated} ]]; then
            params+=(--merged {params.collapsed} --merged {params.collapsed_truncated})
        fi
        spades.py --memory {resources.mem_gb} -1 {input.forward} -2 {input.reverseR} -s {input.singleton} "${{params[@]}}" --threads {threads} -k $kmers -o {params.prefix} > {log.stdout} 2> {log.stderr}
        rm -rf {params.prefix}/corrected
        """

rule unicycler:
    input:
        forward = "results/{sample}/adapterremoval/{sample}_1.fastq",
        reverseR = "results/{sample}/adapterremoval/{sample}_2.fastq"
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
    resources:
        mem_mb = config["mem_mb"]
    shell:
        """
        unicycler -1 {input.forward} -2 {input.reverseR} -o {params.prefix} > {log.stdout} 2> {log.stderr}
        """

# Cloning the cdhit repository because psi-cd-hit.pl is not available in bioconda

rule install_cdhit:
    output:
        "results/bin/cdhit/psi-cd-hit/psi-cd-hit.pl"
    conda:
        "envs/cdhit.yaml"
    shell:
        """
        rm -rf results/bin/cdhit
        git clone https://github.com/weizhongli/cdhit.git results/bin/cdhit > /dev/null 2> /dev/null
        """

rule cdhit:
    input:
        cdhit = "results/bin/cdhit/psi-cd-hit/psi-cd-hit.pl",
        edena = "results/{sample}/edena/{sample}_contigs.fasta",
        spades = "results/{sample}/spades/scaffolds.fasta",
        unicycler = "results/{sample}/unicycler/assembly.fasta"
    params:
        identity = config['cdhit']['identity'],
        program = config['cdhit']['program'],
        circle = config['cdhit']['circle'],
        concat = "results/{sample}/cdhit/concat.fasta"
    output:
        "results/{sample}/cdhit/contigs.fasta"
    log:
        stdout = "results/{sample}/cdhit/log-stdout.txt",
        stderr = "results/{sample}/cdhit/log-stderr.txt"
    benchmark:
        "results/{sample}/cdhit/benchmark.txt"
    resources:
        mem_mb = config["mem_mb"]
    conda:
        "envs/cdhit.yaml"
    shell:
        """
        cat {input.edena} {input.spades} {input.unicycler} > {params.concat}
        {input.cdhit} -i {params.concat} -o {output} -c {params.identity} -prog {params.program} -circle {params.circle} > {log.stdout} 2> {log.stderr}
        """

rule quast:
    input:
        edena = "results/{sample}/edena/{sample}_contigs.fasta",
        spades = "results/{sample}/spades/scaffolds.fasta",
        unicycler = "results/{sample}/unicycler/assembly.fasta",
        cdhit = "results/{sample}/cdhit/contigs.fasta",
    params:
        prefix = "results/{sample}/quast",
        reference = config['quast']['reference'],
        genes = config['quast']['genes'],
	cminlength= config['contigs']['minlength']
    output:
        "results/{sample}/quast/report.tsv"
    log:
        stdout = "results/{sample}/quast/log-stdout.txt",
        stderr = "results/{sample}/quast/log-stderr.txt"
    conda:
        "envs/quast.yaml"
    benchmark:
        "results/{sample}/quast/benchmark.txt"
    threads:
        config["threads"]
    resources:
        mem_mb = config["mem_mb"]
    shell:
        """
        params=()
        if [ -f {params.reference} ]; then
            params+=(-r {params.reference})
        fi
        if [ -f {params.genes} ]; then
            params+=(-g {params.genes})
        fi
        quast.py "${{params[@]}}" -L -t {threads} -m {params.cminlength} -o {params.prefix} {input.edena} {input.spades} {input.unicycler} {input.cdhit} > {log.stdout} 2> {log.stderr}
        """

rule mob_recon:
    input:
        "results/{sample}/cdhit/contigs.fasta"
    params:
        output = "results/{sample}/mob_recon"
    output:
        "results/{sample}/mob_recon/chromosome.fasta"
    log:
        stdout = "results/{sample}/mob_recon/log-stdout.txt",
        stderr = "results/{sample}/mob_recon/log-stderr.txt"
    conda:
        "envs/mobsuite.yaml"
    benchmark:
        "results/{sample}/mob_recon/benchmark.txt"
    threads:
        config["threads"]
    shell:
        """
        mob_recon --force -u -c -t -n {threads} -i {input} -o {params.output} > {log.stdout} 2> {log.stderr}
        """

rule filter_contigs:
    input:
        "results/{sample}/mob_recon/chromosome.fasta"
    output:
        "results/{sample}/mob_recon/chromosome_filtered.fasta"
    conda:
        "envs/filter_contigs.yaml"
    params:
        length = config['contigs']['minlength']
    shell:
        """
        ./scripts/filter_contigs.py {params.length} {input} {output}
        """

# Cloning the prokka repository because tbl2asn in bioconda is old and throws errors

rule install_prokka:
    output:
        "results/bin/prokka/binaries/linux/tbl2asn"
    conda:
        "envs/prokka.yaml"
    shell:
        """
        rm -rf results/bin/prokka
        git clone https://github.com/tseemann/prokka.git results/bin/prokka > /dev/null 2> /dev/null
        """

rule prokka:
    input:
        prokka = "results/bin/prokka/binaries/linux/tbl2asn",
        chromosome = "results/{sample}/mob_recon/chromosome.fasta"
    params:
        outdir = "results/{sample}/prokka",
        prefix = "{sample}",
        prokka = "results/bin/prokka/binaries/linux"
    output:
        "results/{sample}/prokka/{sample}.gbk"
    log:
        stdout = "results/{sample}/prokka/log-stdout.txt",
        stderr = "results/{sample}/prokka/log-stderr.txt"
    conda:
        "envs/prokka.yaml"
    benchmark:
        "results/{sample}/prokka/benchmark.txt"
    threads:
        config["threads"]
    shell:
        """
        export PATH={params.prokka}:$PATH
        prokka --force --cpus {threads} --outdir {params.outdir} --prefix {params.prefix} {input.chromosome} --centre X --compliant > {log.stdout} 2> {log.stderr}
        """

rule cleardir:
    input:
         "results/{sample}/prokka/{sample}.gbk"
    output:
         "results/{sample}/clear/{sample}.txt"
    conda:
         "envs/cleardir.yaml"
    shell:
        """
        echo "Removing fastq files..."
        rm -rf results/{sample}/extract-file/*.fastq
        rm -rf results/{sample}/adapterremoval/*_[0-9].fastq
        rm -rf results/{sample}/adapterremoval/*_[0-9]_edena.fastq
        echo "Fastq files removed..."
        """
