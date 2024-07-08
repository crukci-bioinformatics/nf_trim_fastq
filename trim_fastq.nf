#!/usr/bin/env nextflow

// enable DSL 2 syntax
nextflow.enable.dsl = 2

process runFastp {
    publishDir "${launchDir}/${params.outputDir}", mode: 'link'

    cpus 8
    memory '32.G'
    time '2h'

    input:
        tuple val(id), path(fastqR1In), path(fastqR2In)

    output:
        tuple val(id), path(fastqR1Out), path(fastqR2Out), emit: trimmedFastq
        path fastpJson, emit: fastpJson
        path fastpHtml, emit: fastpHtml

    script:
        fastqR1Out = "${id}.r_1.trimmed.fastq.gz"
        fastqR2Out = "${id}.r_2.trimmed.fastq.gz"
        fastpHtml = "${id}.fastp.html"
        fastpJson = "${id}.fastp.json"
        """
        fastp \
            -i ${fastqR1In} \
            -I ${fastqR2In} \
            -o ${fastqR1Out} \
            -O ${fastqR2Out} \
            --detect_adapter_for_pe \
            -g \
            -x \
            -j ${fastpJson} \
            -h ${fastpHtml} \
            --thread 8 \
            ${params.fastpOptions}
        """
}

process runFastqc {
    publishDir "${launchDir}/${params.outputDir}", mode: 'link'

    cpus 8
    memory '32.G'
    time '2h'

    input:
        tuple val(id), path(fastqR1In), path(fastqR2In)

    output:
        path '*_fastqc.zip', arity: '2', emit: fastqcZip 
        path '*_fastqc.html', arity: '2'

    script:
        """
        fastqc -t 8 ${fastqR1In} ${fastqR2In}
        """

}

process runMultiqc {
    publishDir "${launchDir}/${params.outputDir}", mode: 'link'

    cpus 1
    memory '8.G'
    time '1h'

    input:
        path reportFiles

    output:
        path 'Trimming_Report.html'
        path 'Trimming_Report_data.zip'

    script:
        """
        multiqc -n Trimming_Report.html -z .
        """
}

/*
 * Write a log message summarising how the pipeline is configured and the
 * locations of reference files that will be used.
 */
params.with
{
    log.info "Input fastq direcotry: ${inputDir}"
    log.info "Timmed fastq directory: ${outputDir}"
    log.info "Fastq glob: ${fastqPattern}"
    log.info "Additional options for fastp: ${fastpOptions}"
}

workflow {
    fastqChannel = Channel.fromFilePairs("${params.inputDir}/${params.fastqPattern}", flat: true)

    runFastp(fastqChannel)
    runFastqc(runFastp.out.trimmedFastq)

    reportsChannel = runFastp.out.fastpJson
        .concat(runFastqc.out.fastqcZip)
        .flatten()
        .toList()

    runMultiqc(reportsChannel)
}
