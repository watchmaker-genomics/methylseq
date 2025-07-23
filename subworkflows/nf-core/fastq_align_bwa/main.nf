//
// Alignment with BWA
//
include { PARABRICKS_FQ2BAM        } from '../../../modules/nf-core/parabricks/fq2bam/main'
include { BWA_MEM                  } from '../../../modules/nf-core/bwa/mem/main'
include { BAM_INDEX_STATS_SAMTOOLS } from '../bam_sort_stats_samtools/main'

workflow FASTQ_ALIGN_BWA {
    take:
    ch_reads        // channel (mandatory): [ val(meta), [ path(reads) ] ]
    ch_bwamem_index // channel (mandatory): [ val(meta2), path(index) ]
    ch_fasta        // channel (optional) : [ val(meta3), path(fasta) ]
    val_sort_bam    // boolean (mandatory): true or false
    use_gpu         // boolean (mandatory): whether to use GPU or CPU for bwamem alignment

    main:
    ch_versions = Channel.empty()

    //
    // Map reads with BWA
    //
    if (use_gpu) {
        /*
        * Align with parabricks GPU enabled fq2bammeth implementation of bwameth
        */
        PARABRICKS_FQ2BAM (
            ch_reads,
            ch_fasta,
            ch_bwamem_index,
            [], // interval file
            [], // known sites
            'bam' // output format
        )
        ch_alignment = PARABRICKS_FQ2BAM.out.bam
        ch_versions  = ch_versions.mix(PARABRICKS_FQ2BAM.out.versions)
    } else {
        BWA_MEM ( ch_reads, ch_bwamem_index, ch_fasta, val_sort_bam )
        ch_alignment = BWA_MEM.out.bam
        ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())
    }

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_INDEX_STATS_SAMTOOLS ( ch_alignment, ch_fasta )
    ch_versions = ch_versions.mix(BAM_INDEX_STATS_SAMTOOLS.out.versions)

    emit:
    // bam_orig = BWA_MEM.out.bam                       // channel: [ val(meta), path(bam) ]
    // bam      = BAM_INDEX_STATS_SAMTOOLS.out.bam      // channel: [ val(meta), path(bam) ]
    bam      = ch_alignment                          // channel: [ val(meta), path(bam) ]
    bai      = BAM_INDEX_STATS_SAMTOOLS.out.bai      // channel: [ val(meta), path(bai) ]
    csi      = BAM_INDEX_STATS_SAMTOOLS.out.csi      // channel: [ val(meta), path(csi) ]
    stats    = BAM_INDEX_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), path(stats) ]
    flagstat = BAM_INDEX_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), path(flagstat) ]
    idxstats = BAM_INDEX_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), path(idxstats) ]

    versions = ch_versions                          // channel: [ path(versions.yml) ]
}
