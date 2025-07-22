/*
 * TAPS methylation conversion subworkflow
 *
 * Uses Rastair to assess C->T conversion as a readout for methylation in a genome-wide basis
 */

include { RASTAIR_CALL              } from '../../../modules/nf-core/rastair/call/main'
include { RASTAIR_MBIAS             } from '../../../modules/nf-core/rastair/mbias/main'

workflow TAPS_CONVERSION {

    take:
    ch_bam                 // channel: [ val(meta), [ bam ] ] ## BAM from alignment
    ch_bai                 // channel: [ val(meta), [ bai ] ] ## BAI from alignment
    ch_fasta               // channel: [ [:], /path/to/genome.fa]
    ch_fasta_index         // channel: [ val(meta), /path/to/genome.fa.fai]

    main:

    ch_rastair_mbias = Channel.empty()
    ch_rastair_call  = Channel.empty()
    ch_versions      = Channel.empty()

    log.info "Running TAPS conversion module with Rastair to assess C->T conversion as a readout for methylation."

    RASTAIR_MBIAS (
        ch_bam,
        ch_bai,
        ch_fasta.map{ it[1] },
        ch_fasta_index.map{ it[1] },
    )
    ch_rastair_mbias = RASTAIR_MBIAS.out.txt // channel: [ val(meta), [ txt ] ]
    ch_versions      = ch_versions.mix(RASTAIR_MBIAS.out.versions.first())

    RASTAIR_CALL (
        ch_bam,
        ch_bai,
        ch_fasta.map{ it[1] },
        ch_fasta_index.map{ it[1] },
    )
    ch_rastair_call = RASTAIR_CALL.out.txt // channel: [ val(meta), [ txt ] ]
    ch_versions     = ch_versions.mix(RASTAIR_CALL.out.versions.first())

    emit:
    mbias       = ch_rastair_mbias   // channel: [ val(meta), path("*.txt") ]
    call        = ch_rastair_call    // channel: [ val(meta), path("*.txt") ]
    versions    = ch_versions // channel: path("*.version.txt")
}
