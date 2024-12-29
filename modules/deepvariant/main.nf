process deepvariant {
    label 'high_memory'
    tag "${bam.simpleName}"
    publishDir params.deepvariant_output_dir, mode: 'copy'
    
    container "google/deepvariant:1.8.0"
    
    input:
        path ref            // Reference genome
        path ref_index      // reference index
        path bam           // Input BAM file
        val threads        // Number of shards/threads
        
    
    output:
        path "${bam.baseName}.vcf.gz", emit: vcf
        path "${bam.baseName}.vcf.gz.tbi", emit: vcf_tbi
        path "${bam.baseName}.g.vcf.gz", emit: gvcf
        path "${bam.baseName}.g.vcf.gz.tbi", emit: gvcf_tbi
        
    script:
    """
    /opt/deepvariant/bin/run_deepvariant \
        --model_type PACBIO \
        --ref ${ref} \
        --reads ${bam} \
        --output_vcf ${bam.baseName}.vcf.gz \
        --output_gvcf ${bam.baseName}.g.vcf.gz \
        --regions chr20 \
        --num_shards ${threads} 
    """
}