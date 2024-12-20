process mosdepth {
    publishDir params.aligned_output_dir, mode: 'copy'
    tag "$bam"
    
    input:
    tuple path(bam), path(bam_index)

    output:
    path "${prefix}.mosdepth.summary.txt", emit: summary
    path "${prefix}.regions.bed.gz", emit: region_bed

    script:
    prefix = bam.baseName
    """
    mosdepth --version

    mosdepth \
        --threads ${task.cpus - 1} \
        --by 500 \
        --no-per-base \
        --use-median \
        ${prefix} \
        ${bam}
    """
}