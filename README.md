# nf_trim_fastq

The workflow trims fastq files for adapter content etc. using `fastp`. 

https://github.com/OpenGene/fastp

It then runs `fastqc` on the trimmed files to generate quality control reports, 
and gathers both the `fastqc` and `fastp` reports into a single report using 
`multiqc`.

https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

https://multiqc.info/

## Usage

```
nextflow run crukci-bioinformatics/nf_trim_fastq
```

## Conda environment

The workflow uses a conda environment to manage dependencies.
The environment file is `conda.yaml`.

## Input

The workflow expects a directory containing fastq files. The directory should
be specified using the `--inputDir` parameter. By default, the workflow expects
the fastq files to be in a directory called `fastq` in the current working
directory.

Additional  parameters are:

* `--outputDir` - the directory to write the trimmed fastq files to. By default,
  the trimmed files and reports are written to a directory called
  `fastq_trimmed` in the current working directory.
* `--fastqPattern` - a glob pattern to match the fastq files. By default, the
  workflow will match files ending in `r_{1,2}.fq.gz`. The `r_{1,2}` part of
  the pattern detects the read 1 and read 2 files of paired end data.
* `--fastpOptions` - additional options to pass to `fastp`. By default, the
  workflow uses:
```
    --detect_adapter_for_pe \
    -g \
    -x
```


    
