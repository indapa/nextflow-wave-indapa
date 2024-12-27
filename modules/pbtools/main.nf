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


process hificnv {
    publishDir params.cnv_output_dir, mode: 'copy'
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam), path(bam_index)
    path reference
    path reference_index
    path exclude_bed
    path expected_bed
    val(cpus)

    output:
    tuple val(sample_id), path("hificnv.${sample_id}.vcf.gz"),  emit: hifi_vcf
    tuple val(sample_id), path("hificnv.${sample_id}.depth.bw"), emit: hifi_bigwig
    tuple val(sample_id), path("hificnv.${sample_id}.copynum.bedgraph"), emit: hifi_bedgraph
    tuple val(sample_id), path("hificnv.log"), emit: hifi_log

    script:
    """
    set -euo pipefail

    hificnv --version
    hificnv --bam ${bam} \\
    --ref ${reference} \\
    --exclude  ${exclude_bed} \\
    --expected-cn ${expected_bed} \\
    --threads ${cpus} 
    """
}


process trgt {
    publishDir params.trgt_output_dir, mode: 'copy'
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam), path(bam_index)
    path reference
    path reference_index
    path tandem_repeat_bed
    val(karyotype)
    val(cpus)

    output:
    tuple val(sample_id), path("${sample_id}.trgt.spanning.sorted.bam"), path("${sample_id}.trgt.spanning.sorted.bam.bai"), emit: spanning_reads
    tuple val(sample_id), path("${sample_id}.trgt.sorted.vcf.gz"), path("${sample_id}.trgt.sorted.vcf.gz.tbi"), emit: repeat_vcf

    script:
    """
    set -euo pipefail

    trgt --version
    trgt genotype \\
        --threads ${cpus} \\
        --karyotype ${karyotype} \\
        --genome ${reference} \\
        --repeats ${tandem_repeat_bed} \\
        --reads ${bam} \\
        --output-prefix ${sample_id}.trgt

    bcftools --version
    bcftools sort \\
        --output-type z \\
        --output-file ${sample_id}.trgt.sorted.vcf.gz \\
        ${sample_id}.trgt.vcf.gz

    bcftools index \\
        --threads ${cpus} \\
        --tbi \\
        ${sample_id}.trgt.sorted.vcf.gz

    samtools --version
    samtools sort \\
        -@ ${cpus} \\
        -o ${sample_id}.trgt.spanning.sorted.bam \\
        ${sample_id}.trgt.spanning.bam

    samtools index \\
        -@ ${cpus} \\
        ${sample_id}.trgt.spanning.sorted.bam
    """
}

/* next set of processes deal with SVs */

process extractRegions {
    input:
    tuple val(sample_id), path(bam), path(bam_index)
    
       
    output:
    tuple val(sample_id), stdout

    script:
    """
    samtools view -H ${bam} | grep '^@SQ' | grep -v chrUn | grep -v random | cut -f2 | cut -d':' -f2
    """
}

process pb_discover {
    publishDir params.sv_output_dir, mode: 'copy'
    tag "$sample_id"
    
    input:
    tuple val(sample_id), val(region), path(bam), path(bam_index)
    path trf_bed
       
    output:
    tuple val(sample_id), path("${sample_id}.${region}.svsig.gz"), emit: pb_discover
        
    script:
    """
    set -euo pipefail

    pbsv --version
    pbsv discover --hifi --tandem-repeats ${trf_bed} --region ${region} ${bam} ${sample_id}.${region}.svsig.gz
    """
}

process pb_call {
    publishDir params.sv_output_dir, mode: 'copy'
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(svsig_files)
    path reference
       
    output:
    tuple val(sample_id), path("${sample_id}.pbsv.vcf"), emit: pb_call
        
    script:
    """
    set -euo pipefail

    pbsv --version
    pbsv call -j 8 ${reference} ${svsig_files.join(' ')} ${sample_id}.pbsv.vcf

    """

    
}