/**
This module describes the Rastair process 'call' to
run Rastair and assessing C->T conversion
as a readout for methylation
in a genome-wide basis
*/

process RASTAIR_CALL {
    label 'process_medium'
    container "docker.io/sbludwig/rastair:version-0.8.2"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(bai)
    path(fasta)
    path(fai)

    output:
    tuple val(meta), path("*.rastair_call.txt"),    emit: txt
    tuple val(meta), path("*_methylkit.txt.gz"),    emit: gz
    path "versions.yml",                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    rastair call \\
        --threads ${task.cpus} \\
        --nOT ${nOT_clip} \\
        --nOB ${nOB_clip} \\
        --fasta-file ${fasta} \\
        ${bam} | tee ${prefix}.rastair_call.txt | /app/scripts/rastair_call_to_methylkit.sh | gzip -c > ${prefix}.rastair_methylkit.txt.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rastair: \$(rastair --version)
    END_VERSIONS
    """
}
