# Assembly HiSeq

A snakemake pipeline to assemble sequencing data produced by Illumina HiSeq

## Installation

The only requirement to run this pipeline is **Anaconda 3**. Please, install the latest version of **Anaconda 3** on your machine beforehand ([https://www.anaconda.com/products/individual](https://www.anaconda.com/products/individual)).

With Anaconda installed, we can create a new environment with snakemake by running the following command:

```
conda create -c conda-forge -c bioconda -n assembly-hiseq snakemake -y
```

For some reason, Snakemake cannot activate conda environments automatically, hence you need to copy the activate binary from your conda main folder to the environment you created to use this pipeline. For example:

```
cp /home/fabio/anaconda3/bin/activate /home/fabio/anaconda3/envs/assembly-hiseq/bin/activate
```

To copy the contents of this repository, run the following command:

```
git clone https://github.com/softwarecomputationalbiologyufpa/assembly-hiseq.git
```

After cloning the repository, it is possible to enter the newly created directory and modify the config file with the appropriate parameters and samples:

```
cd assembly-hiseq
pwd
```

Copy the output of the pwd command and update the workdir on the config.yaml file.

```
nano config.yaml
```

It is a good idea to run the pipeline with the test sample provided for the first time to test your installation. However, to add your samples, simply append them to the samples section of the config file.

```
samples:
    'NG-19302_B208_lib327116_6377_1':
        forward: 'data/NG-19302_B208_lib327116_6377/NG-19302_B208_lib327116_6377_1_1.fastq.gz'
        reverse: 'data/NG-19302_B208_lib327116_6377/NG-19302_B208_lib327116_6377_1_2.fastq.gz'
    'another_sample':
        forward: 'data/another_sample/another_sample_1.fastq.gz'
        reverse: 'data/another_sample/another_sample_2.fastq.gz'
```

**It is important that you write the identifier of the sample with the common prefix between the paired FASTQ files, as shown above, otherwise FastQC will throw an error.**

## Running

After a successful installation, we can activate the newly created environment and run the pipeline (please don't forget to modify the config file with your workdir, samples and parameters in advance).

```
conda activate assembly-hiseq
snakemake -kpr --use-conda --cores 12
```

The parameter --use-conda is necessary to indicate that conda will be used to manage the software dependencies of the pipeline, while the parameter --cores tells Snakemake how many cpus can be used to assemble the samples (the more cpus you can spare, the faster the assemblies will be completed).

The assemblies for each sample are saved in the results folder.
