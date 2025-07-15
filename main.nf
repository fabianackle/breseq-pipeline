#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/* Modules */
process TRIMMOMATIC {
    conda "bioconda::trimmomatic=0.39"
    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_R1_trim_paired.fastq.gz"), path("${sample_id}_R2_trim_paired.fastq.gz"), emit: reads
    path("${sample_id}_trim_out.log"), emit: log

    script:
    def illumina_adapters = file(params.illumina_adapters)
    """
    trimmomatic PE \\
        -threads ${task.cpus} \\
        -phred33 \\
        ${reads[0]} ${reads[1]} \\
        ${sample_id}_R1_trim_paired.fastq.gz ${sample_id}_R1_trim_unpaired.fastq.gz \\
        ${sample_id}_R2_trim_paired.fastq.gz ${sample_id}_R2_trim_unpaired.fastq.gz \\
        ILLUMINACLIP:${illumina_adapters}:2:30:10:2 \\
        SLIDINGWINDOW:4:12 \\
        MINLEN:100 \\
        2> ${sample_id}_trim_out.log
    """
}

process BRESEQ {
    conda "bioconda::breseq=0.39.0"
    tag "$sample_id"

    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    path("$sample_id"), emit: results

    script:
    def reference_seq = file(params.reference_seq)
    """
    breseq \\
        -j ${task.cpus} \\
        -r ${reference_seq} \\
        -o ${sample_id} \\
        ${read1} ${read2}
    """
}

process FASTQC {
    conda "bioconda::fastqc=0.12.1"
    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    path("*.{html,zip}"), emit: stats

    script:
    """
    fastqc --threads $task.cpus \
        ${reads[0]} ${reads[1]}
    """
}

process MULTIQC {
    conda "bioconda::multiqc=1.28"
    tag "MultiQC"

    publishDir params.outdir, mode: 'copy'

    input:
    path("*")

    output:
    path("multiqc_report.html"), emit: report
    path("multiqc_data"), emit: data

    script:
    def multiqc_config = file("multiqc_config.yml")
    """
    multiqc --config ${multiqc_config} .
    """
}

/* Workflow */
workflow {
    reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)

    multiqc_files_ch = Channel.empty()

    TRIMMOMATIC(reads_ch)
    FASTQC(reads_ch)
    BRESEQ(TRIMMOMATIC.out.reads)
    
    multiqc_files_ch = multiqc_files_ch.mix(FASTQC.out.stats)
    multiqc_files_ch = multiqc_files_ch.mix(TRIMMOMATIC.out.log)
    multiqc_files_ch = multiqc_files_ch.collect()
    MULTIQC(multiqc_files_ch) 
}