process pbmm2_align {
    label 'high_memory'
    publishDir "${params.aligned_output_dir}", mode: 'copy'
    tag "$sample_id"

    input:
    path reference
    tuple val(sample_id), path(bam)
    val threads
    val sort_threads

    output:
    tuple path("${sample_id}.aligned.bam"), path("${sample_id}.aligned.bam.bai"), emit: aligned_bam
    //path "${$sample_id}.read_length_and_quality.tsv", emit: bam_rl_qual
    
    script:
    """
    pbmm2 --version
    pbmm2 align \\
        --sort \\
        -j $threads \\
        -J $sort_threads \\
        --preset HIFI \\
        --sample ${sample_id} \\
        --log-level INFO \\
        --unmapped \\
        --bam-index BAI \\
        $reference \\
        $bam \\
        ${sample_id}.aligned.bam

    
    """
}


process cpg_pileup {
    
    publishDir "${params.cpg_output_dir}", mode: 'copy'
    tag "$bam"

    input:
    tuple path(bam), path(bam_index)
    path reference
    path reference_index
    path cpgmodel

    output:
    path "${bam.baseName}.*.bed", emit: pileup_beds
    path "${bam.baseName}.*.bw", emit: pileup_bigwigs

    script:
    """
    aligned_bam_to_cpg_scores \\
        --threads 4 \\
        --bam ${bam} \\
        --ref ${reference} \\
        --output-prefix ${bam.baseName} \\
        --min-mapq 1 \\
        --min-coverage 10 \\
        --model /usr/local/bin/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite
    """
}
