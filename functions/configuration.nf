/*
 * Write a log message summarising how the pipeline is configured and the
 * locations of reference files that will be used.
 */

def displayParameters(params)
{
    params.with
    {
        log.info "Input fastq direcotry: ${inputDir}"
        log.info "Timmed fastq directory: ${outputDir}"
        log.info "Fastq glob: ${fastqPattern}"
        log.info "Additional options for fastp: ${fastpOptions}"
    }
}