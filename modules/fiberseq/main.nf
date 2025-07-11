process fibertools_extract  {

    publishDir params.fibertools_output_dir, mode: 'copy'
    tag "$sample_id"
    container 'quay.io/biocontainers/fibertools-rs:0.6.2--h3b373d1_0'


    input:
    tuple val(sample_id), path(bam), path(bai)
    output:
    path "${sample_id}.m6a.bed.gz", emit: m6a_bed 
    

    script:

    """

    ft --version
    ft extract ${bam} --m6a ${sample_id}.m6a.bed.gz
    

    """

}

process fibertools_extract_all  {

    tag "$sample_id"
    container 'quay.io/biocontainers/fibertools-rs:0.6.2--h3b373d1_0'


    input:
    tuple val(sample_id), path(bam), path(bai)
    output:
    tuple val(sample_id), path("${sample_id}.all.bed.gz"), emit: ft_all_bed 


    script:

    """

    ft --version
    ft extract ${bam} --all  ${sample_id}.all.bed.gz --simplify
    
    """

    stub:
    """
    touch ${sample_id}.all.bed.gz
    """
}


process m6a_extract_parquet {
    tag "${input_parquet}"
    publishDir "${params.fibertools_output_dir}", mode: 'copy', overwrite: true
    label 'medium_low_memory'
    
    container 'indapa/indapa-methylation-data-processing:latest'
    
    input:
    path(input_parquet)

    output:
    path "*.csv"

    script:

    """

    /opt/bin/m6a_rust ${input_parquet} .

    """

    stub:
    """
    echo "/opt/bin/m6a_rust ${input_parquet} ."
    touch ${input_parquet}.csv
    """

}

process fibertools_extract_all_to_parquet {

     /* convert fibertools extract bed file to parquet format */
    tag "${bed_file}"
    label 'high_memory_spot'
    publishDir "${params.fibertools_output_dir}", mode: 'copy', overwrite: true

    container 'indapa/indapa-methylation-data-processing:latest'

    input:
    tuple val(sample_id), path(bed_file)

    output:
    path "${sample_id}.all.parquet", emit: ft_all_parquet

    

    script:
    
    """
    python /opt/bin/fiberseq_all_bedToParquet.py -bed ${bed_file} -sampleID ${sample_id}
    """

    stub:
    """
    touch ${sample_id}.all.parquet
    """
}