/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/methylseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/methylseq
    Website: https://nf-co.re/methylseq
    Slack  : https://nfcore.slack.com/channels/methylseq
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTA_INDEX_BISMARK_BWAMETH   } from './subworkflows/nf-core/fasta_index_bismark_bwameth/main'
include { BWA_INDEX                     } from './modules/nf-core/bwa/index/main'
include { PIPELINE_INITIALISATION       } from './subworkflows/local/utils_nfcore_methylseq_pipeline'
include { PIPELINE_COMPLETION           } from './subworkflows/local/utils_nfcore_methylseq_pipeline'
include { getGenomeAttribute            } from './subworkflows/local/utils_nfcore_methylseq_pipeline'
include { METHYLSEQ                     } from './workflows/methylseq/'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
params.fasta         = getGenomeAttribute('fasta')
params.fasta_index   = getGenomeAttribute('fasta_index')
params.bwameth_index = getGenomeAttribute('bwameth')
params.bwamem_index  = getGenomeAttribute('bwa')
params.bismark_index = params.aligner == 'bismark_hisat' ? getGenomeAttribute('bismark_hisat2') : getGenomeAttribute('bismark')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_METHYLSEQ {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()

    //
    // Initialize file channels or values based on params
    //
    ch_fasta                = params.fasta         ? Channel.fromPath(params.fasta).map{ it -> [ [id:it.baseName], it ] } : Channel.empty()
    ch_or_val_fasta_index   = params.fasta_index   ? Channel.fromPath(params.fasta_index).map{ it -> [ [id:it.baseName], it ] } : []
    ch_or_val_bismark_index = params.bismark_index ? Channel.fromPath(params.bismark_index).map{ it -> [ [id:it.baseName], it ] } : []
    ch_or_val_bwameth_index = params.bwameth_index ? Channel.fromPath(params.bwameth_index).map{ it -> [ [id:it.baseName], it ] } : []
    // ch_or_val_bwamem_index  = params.bwamem_index  ? Channel.fromPath(params.bwamem_index).map{ it -> [ [id:it.baseName], it ] } : []

    //
    // SUBWORKFLOW: Prepare any required reference genome indices
    //
    FASTA_INDEX_BISMARK_BWAMETH(
        ch_fasta,
        ch_or_val_fasta_index,
        ch_or_val_bismark_index,
        ch_or_val_bwameth_index,
        params.aligner,
        params.collecthsmetrics
    )
    ch_versions = ch_versions.mix(FASTA_INDEX_BISMARK_BWAMETH.out.versions)

    // Here we could check if ch_or_val_bwamem_index is empty or not
    // if it is empty, we can run the BWA_INDEX subworkflow
    // if it is not empty, we need to validate the index
    BWA_INDEX (
        ch_fasta
    )
    ch_versions = ch_versions.mix(BWA_INDEX.out.versions)

    //
    // WORKFLOW: Run pipeline
    //

    METHYLSEQ (
        samplesheet,
        ch_versions,
        FASTA_INDEX_BISMARK_BWAMETH.out.fasta,
        FASTA_INDEX_BISMARK_BWAMETH.out.fasta_index,
        FASTA_INDEX_BISMARK_BWAMETH.out.bismark_index,
        FASTA_INDEX_BISMARK_BWAMETH.out.bwameth_index,
        BWA_INDEX.out.index,
    )
    ch_versions = ch_versions.mix(METHYLSEQ.out.versions)

    emit:
    multiqc_report = METHYLSEQ.out.multiqc_report // channel: [ path(multiqc_report.html )  ]
    versions       = ch_versions                  // channel: [ path(versions.yml) ]

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_METHYLSEQ (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_METHYLSEQ.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
