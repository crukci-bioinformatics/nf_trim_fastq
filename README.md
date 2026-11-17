# nf_trim_fastq

The workflow trims fastq files for adapter content etc. using `fastp`. 

https://github.com/OpenGene/fastp

It then runs `fastqc` on the trimmed files to generate quality control reports, 
and gathers both the `fastqc` and `fastp` reports into a single report using 
`multiqc`.

https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

https://multiqc.info/

## Requirements

Nextflow 25.10 or later. The pipeline uses the strict Nextflow syntax and
workflow output definitions.

## Usage

```
nextflow run crukci-bioinformatics/nf_trim_fastq
```

Profiles:

* `standard` (default) - run locally with up to 4 CPUs and 8 GB of memory.
* `bigserver` - run locally with up to 50 CPUs and 128 GB of memory.
* `cluster` - submit jobs to SLURM.
* `epyc` - submit jobs to the SLURM `epyc` queue.
* `test` - a stub run on placeholder inputs; see [Testing](#testing).

## Conda environment

The workflow uses a conda environment to manage dependencies.
The environment file is `conda.yaml`:

* fastp 1.4.0
* FastQC 0.13.0
* MultiQC 1.35

## Input

The workflow expects a directory containing fastq files. The directory should
be specified using the `--inputDir` parameter. By default, the workflow expects
the fastq files to be in a directory called `fastq` in the current working
directory.

Additional  parameters are:

* `--outputDir` - the directory to write the trimmed fastq files to. By default,
  the trimmed files and reports are written to a directory called
  `fastq_trimmed` in the current working directory. Published files are hard
  links to the files in the Nextflow work directory.
* `--fastqPattern` - a glob pattern to match the fastq files. By default, the
  workflow will match files ending in `r_{1,2}.fq.gz`. The `r_{1,2}` part of
  the pattern detects the read 1 and read 2 files of paired end data.
* `--fastpOptions` - additional options to pass to `fastp`. These are added to
  the options the workflow always uses:
```
    --detect_adapter_for_pe \
    -g \
    -x
```

## Output

Written to `--outputDir`:

* `<id>.r_1.trimmed.fastq.gz`, `<id>.r_2.trimmed.fastq.gz` - trimmed reads
* `<id>.fastp.json`, `<id>.fastp.html` - fastp reports
* `<id>.r_{1,2}.trimmed_fastqc.{zip,html}` - FastQC reports
* `Trimming_Report.html`, `Trimming_Report_data.zip` - MultiQC report

## Testing

A stub run checks the wiring of the whole workflow in a few seconds, without
running any of the tools or building the conda environment:

```
nextflow run main.nf -stub-run -profile test
```

It uses the empty placeholder files in `test/fastq/` and writes to
`test_outputs/`. See [test/README.md](test/README.md).
