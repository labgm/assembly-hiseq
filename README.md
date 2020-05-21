# Assembly HiSeq

A snakemake pipeline to assemble sequencing data produced by Illumina HiSeq

## Installation

The only requirement to run this pipeline is conda. Please, install conda on your machine beforehand ([https://www.anaconda.com/products/individual](https://www.anaconda.com/products/individual)).

With conda installed, we can create a new environment with snakemake by running the following command:

```
conda create -c conda-forge -c bioconda -n assembly-hiseq snakemake -y
```

To copy the contents of this repository, run the following command:

```
git clone https://github.com/softwarecomputationalbiologyufpa/assembly-hiseq.git
```

After cloning the repository, it is possible to enter the newly created directory and modify the config file with the appropriate parameters and samples:

```
cd assembly-hiseq
nano config.yaml
```

## Running

After a successful installation, we can activate the newly created environment and run the pipeline (please don't forget to modify the config file with your samples and parameters in advance).

```
conda activate assembly-hiseq
snakemake -s Snakefile --configfile config.yaml --use-conda --cores 12
```

The -s parameter is used to select the script containing the workflow rules. Similarly, the --configfile parameter is used to select the config file with parameters and samples to be processed, while the --use-conda is necessary to indicate that conda will be used to manage the software dependencies of the pipeline. Finally, the --cores parameters is used to indicate how many cpus can be used by Snakemake to assemble the samples (the more cpus you can use, the faster the assemblies will be completed).

The assemblies for each sample are be saved in the results folder.

## Troubleshooting

If Snakemake cannot activate conda environments, you can just copy the activate binary from your conda main folder to the environment you created to use this pipeline. For example:

```
cp /home/fabio/anaconda3/bin/activate /home/fabio/anaconda3/envs/assembly-hiseq/bin/activate
```
