# Assembly HiSeq
A snakemake pipeline to assemble sequencing data produced by Illumina HiSeq

## Installation
```
conda create -c conda-forge -c bioconda -n assembly-hiseq snakemake
```

## Running
```
conda activate assembly-hiseq
snakemake -s Snakefile --configfile config.yaml --use-conda -kpr --cores 12
```

## Troubleshooting

```
cp /home/fabio/anaconda3/bin/activate /home/fabio/anaconda3/envs/assembly-hiseq/bin/activate
```
