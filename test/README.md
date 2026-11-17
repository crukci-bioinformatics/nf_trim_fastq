# Test fixtures

Empty placeholder FASTQ files for a stub run of the pipeline:

```
nextflow run main.nf -stub-run -profile test
```

Two paired-end samples (`sampleA`, `sampleB`) match the default `*r_{1,2}.fq.gz`
glob. The stub run exercises the whole workflow graph (fastp → FastQC → MultiQC)
and publishing without running any of the tools or building the conda
environment. Output is written to `test_outputs/` in the launch directory.

The stub outputs are empty files: a stub run checks the wiring of the pipeline,
not the trimming itself.
