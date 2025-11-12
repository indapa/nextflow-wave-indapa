process multiqc {
    publishDir "${params.multiqc_output_dir}", mode: 'copy', overwrite: true
    label 'low_memory'
    container 'quay.io/biocontainers/multiqc:1.28--pyhdfd78af_0'

    input:
    path input_files
    
    output:
    path 'multiqc_report.html'
    path 'multiqc_data'

    script:
    """
    multiqc .
    """
}

process multiqc_downsample {
    publishDir "${params.multiqc_downsampled_output_dir}", mode: 'copy', overwrite: true
    label 'low_memory'
    container 'quay.io/biocontainers/multiqc:1.28--pyhdfd78af_0'

    input:
    path input_files
    
    output:
    path 'multiqc_report.html'
    path 'multiqc_data'

    script:
    """
    multiqc .
    """
}