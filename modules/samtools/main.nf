process index_bam {

    /* index bam file */

    tag "${sample_id}"

    publishDir params.bam_index_output_dir, mode: 'copy'

    container 'community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c'

    input:
    tuple val(sample_id), path(bam)

    output:
    path "${bam}.bai", emit: bai

    script:
    """

    samtools index ${bam}

    """
}

/*
process extract_unmapped_reads_fasta {
    
    // extract unmapped reads from bam file and convert to fasta 

    tag "${sample_id}"
    label 'low_memory'

    publishDir params.unmapped_reads_output_dir, mode: 'copy', overwrite: true

    container 'community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    path "${sample_id}.unmapped_reads.fa.gz", emit: unmapped_fasta

    script:
    """

    samtools view -f 4 "${bam}" | awk '{print ">" $1 "\n" $10}' | gzip > "${sample_id}.unmapped_reads.fa.gz"

    """

    stub:
    """
    touch ${sample_id}.unmapped_reads.fa.gz
    """
}
*/
// calculate read lengths for the bam file


process read_lengths {

    /* calculate read lengths for the bam file */

    tag "${bam}"
    label 'low_memory'

    publishDir params.bam_stats_output_dir, mode: 'copy'

    container 'community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c'

    input:
    tuple  path(bam), path(bai)

    output:
    path "${bam}.read_lengths.txt", emit: read_lengths

    script:
    """

    samtools view ${bam} | awk '{print length(\$10)}' | sort -n | uniq -c > ${bam}.read_lengths.txt

    """
}
process filter_mapq_primary_alns {


    //filter mapping quality and primary alignments
    label 'low_memory'
    tag "${bam}"

    publishDir params.aligned_filtered_output_dir, mode: 'copy'

    input:
    tuple path(bam), path(bai)

    output:
    path "${bam.simpleName}.filtered.bam", emit: filtered_bam
    path "${bam.simpleName}.filtered.bam.bai", emit: filtered_bai

    container 'community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c'

    script:
    """

    samtools view -b -q 60 -F 0x900 ${bam} > ${bam.simpleName}.filtered.bam
    samtools index ${bam.simpleName}.filtered.bam

    """
}
// generate stats for the bam file

process bam_stats {

    /* generate stats for the bam file */
    label 'low_memory'

    tag "${bam}"

    publishDir params.bam_stats_output_dir, mode: 'copy'

    container 'community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c'

    input:
    tuple path(bam), path(bai)

    output:
    path "${bam}.stats", emit: stats

    script:
    """

    samtools stats --threads 4 ${bam} > ${bam}.stats

    """

    stub:
    """

    touch ${bam}.stats

    """
}


// process to extract out chr20 from the bam file


process extract_chr20 {

    /* extract chr20 from the bam file; helpful for test data */

    tag "${sample_id}"

    publishDir params.aligned_output_dir, mode: 'copy'

    container 'community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    path "${sample_id}.chr20.aligned.bam", emit: chr20_bam
    path "${sample_id}.chr20.aligned.bam.bai", emit: chr20_bai

    script:
    """

    samtools view -b ${bam} chr20 > ${sample_id}.chr20.aligned.bam
    samtools index ${sample_id}.chr20.aligned.bam

    """
}

process subsample_bam {
    
    /* subsample unaligned  bam file for the first 25k reads
       Useful for making test dataset  
     */

    tag "${sample_id}"

    publishDir params.test_data_dir, mode: 'copy'

    container 'community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c'

    input:
    tuple val(sample_id), path(bam)

    output:
    
    path "${bam.simpleName}.subsampled.bam", emit: subsampled_bam
    

    script:
    def basename = bam.simpleName  // Gets name without extension
    """

    samtools view -h ${bam} | head -n 250000  | samtools view -bS - > ${basename}.subsampled.bam
    
    """
}

// process to conver aligned BAM to unaligned BAM

process unaligned_bam {

    /* convert aligned BAM to unaligned BAM */

    tag "${sample_id}"

    publishDir params.unaligned_output_dir, mode: 'copy'

    container 'community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    path "${sample_id}.unaligned.bam", emit: unaligned_bam

    script:
    """

    samtools view -h input.bam | \
   awk '{if (\$0 ~ /^@/) {print} else {\$2="4"; \$3="*"; \$4="0"; \$5="0"; \$6="*"; \$7="*"; \$8="0"; \$9="0"; print}}' | \
    samtools view -b > ${sample_id}.unaligned.bam

    """
}