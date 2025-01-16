process fibertools_extract  {

    publishDir params.fibertools_output_dir, mode: 'copy'
    tag "$sample_id"


    input:
    tuple val(sample_id), path(bam)
    output:
    path "${sample_id}.m6a.bed.gz", emit: m6a_bed 

    script:

    """

    ft --version
    ft extract ${bam} --m6a ${sample_id}.m6a.bed.gz

    """


    

}