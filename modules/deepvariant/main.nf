process deepvariant {
    tag "${bam.simpleName}"
    publishDir params.deepvariant_output_dir, mode: 'copy'
    
    container "google/deepvariant:1.8.0"
    
    input:
        path ref            // Reference genome
        path ref_index      // reference index
        path bam           // Input BAM file
        val threads        // Number of shards/threads
        
    
    output:
        path "${bam.simpleName}.vcf", emit: vcf
        path "${bam.simpleName}.g.vcf", emit: gvcf
        
    script:
    """
    /opt/deepvariant/bin/run_deepvariant \
        --model_type PACBIO \
        --ref ${ref} \
        --reads ${bam} \
        --output_vcf ${bam.simpleName}.vcf \
        --output_gvcf ${bam.simpleName}.g.vcf \
        --num_shards ${threads} 
    """
}