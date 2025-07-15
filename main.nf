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
    path("${sample_id}.quality_read_trimm_info.txt"), emit: log

    script:
    """
    trimmomatic PE \\
        -threads ${task.cpus} \\
        -phred33 \\
        ${reads[0]} ${reads[1]} \\
        ${sample_id}_R1_trim_paired.fastq.gz ${sample_id}_R1_trim_unpaired.fastq.gz \\
        ${sample_id}_R2_trim_paired.fastq.gz ${sample_id}_R2_trim_unpaired.fastq.gz \\
        ILLUMINACLIP:${params.illumina_adapters}:2:30:10:2 \\
        SLIDINGWINDOW:4:12 \\
        MINLEN:100 \\
        2> ${sample_id}.quality_read_trimm_info.txt
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
    path("${sample_id}.breseq.log"), emit: log

    script:
    """
    breseq \\
        -r ${params.reference} \\
        -o ${sample_id} \\
        ${read1} ${read2} \\
        &> ${sample_id}.breseq.log
    """
}

/* Workflow */
workflow {
    reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)

    TRIMMOMATIC(reads_ch)
    BRESEQ(TRIMMOMATIC.out.reads)
}