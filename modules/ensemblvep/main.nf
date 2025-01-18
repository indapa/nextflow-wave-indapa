process annotate_vep {

    tag "${sample_id}"
    container 'quay.io/biocontainers/ensembl-vep:112.0--pl5321h2a3209d_0'
    publishDir params.deepvariant_output_dir, mode: 'copy', pattern: '*.vcf.gz*'

    input:

    path phased_vcf
    tuple val(sample_id), path(haplotagged_bam)
    path pigeon_gtf
    path pigeon_tbi
    path reference
    

    output:
    path "${sample_id}.deepvariant.phased.vep.vcf.gz", emit: vep_phased_deepvariant
    path "${sample_id}.deepvariant.phased.vep.vcf.gz.tbi", emit: vep_phased_deepvariant_tbi

    script:
    """

    

    vep -i ${phased_vcf} -o ${sample_id}.deepvariant.phased.vep.vcf.gz --format vcf --gtf ${pigeon_gtf} --fasta ${reference} --vcf --everything --fork 8 --compress_output bgzip
    tabix -p vcf ${sample_id}.deepvariant.phased.vep.vcf.gz
    """


}