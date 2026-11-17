#!/usr/bin/env nextflow

// Parameter names and types are declared here. Defaults are set in
// nextflow.config.
params {
    // Directory containing the untrimmed fastq files.
    inputDir : String

    // Directory, relative to the launch directory, to which to write the
    // trimmed fastq files and reports.
    outputDir : String

    // Glob pattern to match the read 1/2 fastq files.
    fastqPattern : String

    // Additional options to pass to fastp.
    fastpOptions : String = ""
}

// Trim adapters, poly-G and poly-X tails from paired-end reads.
process runFastp {
    input:
        tuple val(id), path(fastqR1In), path(fastqR2In)
        val fastpOptions

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
            --thread ${task.cpus} \
            ${fastpOptions}
        """

    stub:
        fastqR1Out = "${id}.r_1.trimmed.fastq.gz"
        fastqR2Out = "${id}.r_2.trimmed.fastq.gz"
        fastpHtml = "${id}.fastp.html"
        fastpJson = "${id}.fastp.json"
        """
        touch ${fastqR1Out} ${fastqR2Out} ${fastpHtml} ${fastpJson}
        """
}

// Quality control reports on the trimmed reads.
process runFastqc {
    input:
        tuple val(id), path(fastqR1In), path(fastqR2In)

    output:
        path '*_fastqc.zip', arity: '2', emit: fastqcZip
        path '*_fastqc.html', arity: '2', emit: fastqcHtml

    script:
        """
        fastqc -t ${task.cpus} ${fastqR1In} ${fastqR2In}
        """

    stub:
        // FastQC names its reports after the input with .fastq.gz/.fq.gz removed
        def stems = [fastqR1In, fastqR2In].collect { f -> f.name.replaceAll(/\.f(ast)?q(\.gz)?$/, '') }
        """
        touch ${stems.collect { st -> "${st}_fastqc.zip ${st}_fastqc.html" }.join(' ')}
        """
}

// Gather the fastp and FastQC reports into a single report.
process runMultiqc {
    input:
        path reportFiles

    output:
        path 'Trimming_Report.html', emit: report
        path 'Trimming_Report_data.zip', emit: data

    script:
        """
        multiqc -n Trimming_Report.html -z .
        """

    stub:
        """
        touch Trimming_Report.html Trimming_Report_data.zip
        """
}

workflow {
    main:
    // Summarise how the pipeline is configured
    log.info "Input fastq directory: ${params.inputDir}"
    log.info "Trimmed fastq directory: ${params.outputDir}"
    log.info "Fastq glob: ${params.fastqPattern}"
    log.info "Additional options for fastp: ${params.fastpOptions}"

    def fastqChannel = channel.fromFilePairs("${params.inputDir}/${params.fastqPattern}", flat: true, checkIfExists: true)

    runFastp(fastqChannel, params.fastpOptions)
    runFastqc(runFastp.out.trimmedFastq)

    def reportsChannel = runFastp.out.fastpJson
        .concat(runFastqc.out.fastqcZip)
        .flatten()
        .toList()

    runMultiqc(reportsChannel)

    publish:
    trimmed_fastq = runFastp.out.trimmedFastq.flatMap { _id, r1, r2 -> [r1, r2] }
    fastp_reports = runFastp.out.fastpJson.mix(runFastp.out.fastpHtml)
    fastqc_reports = runFastqc.out.fastqcZip.mix(runFastqc.out.fastqcHtml).flatten()
    multiqc_report = runMultiqc.out.report.mix(runMultiqc.out.data)
}

output {
    trimmed_fastq {
        path '.'
        mode 'link'
    }

    fastp_reports {
        path '.'
        mode 'link'
    }

    fastqc_reports {
        path '.'
        mode 'link'
    }

    multiqc_report {
        path '.'
        mode 'link'
    }
}
