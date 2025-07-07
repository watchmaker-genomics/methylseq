/**
This module describes the Rastair process 'call' to
run Rastair and assessing C->T conversion
as a readout for methylation
in a genome-wide basis
*/

process RASTAIR_CALL {
    label 'process_medium'
    container "docker.io/sbludwig/rastair:version-0.8.1"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(bai)
    tuple val(meta3), path(fasta)
    tuple val(meta4), path(fai)

    output:
    tuple val(meta), path("*.rastair_call.txt"),    emit: txt
    path "versions.yml",                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def nOT_clip = meta.rastair_nOT_clip ?: '0,0,10,0' // [r1_start, r1_end, r2_start, r2_end]
    def nOB_clip = meta.rastair_nOB_clip ?: '0,0,10,0' // [r1_start, r1_end, r2_start, r2_end]

    """
    rastair call \\
        --threads ${task.cpus} \\
        --nOT ${nOT_clip} \\
        --nOB ${nOB_clip} \\
        --fasta-file ${fasta} \\
        ${bam} > ${prefix}.rastair_call.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rastair: \$(rastair --version)
    END_VERSIONS
    """
}
