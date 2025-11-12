


process cnv_vcfToBedB {

    tag "$sample_id"
    container "community.wave.seqera.io/library/bcftools:1.21--4335bec1d7b44d11"
    publishDir params.bedtools_output_dir, mode: 'copy', overwrite: true

    input:
    tuple(val(sample_id), path(vcf))

    output:
    tuple (val(sample_id), path("${sample_id}.B.cnv.bed"), emit: cnv_bedB)


    
    script:
   

    """
    bcftools view -f PASS ${vcf} | bcftools query -f '%CHROM\t%POS\t%INFO/END' - | grep -v _ | grep -v GL  | grep -v KI > ${sample_id}.B.cnv.bed
    """
    
}




process cnv_vcfToBedA {

    tag "$sample_id"
    container "community.wave.seqera.io/library/bcftools:1.21--4335bec1d7b44d11"
    publishDir params.bedtools_output_dir, mode: 'copy', overwrite: true

    input:
    tuple(val(sample_id), path(vcf))

    output:
    tuple (val(sample_id), path("${sample_id}.A.cnv.bed"), emit: cnv_bedA)

    
    script:
   

    """
    bcftools view -f PASS ${vcf} | bcftools query -f '%CHROM\t%POS\t%INFO/END' - | grep -v _ | grep -v GL  | grep -v KI > ${sample_id}.A.cnv.bed
    """




}


process cnv_bed {
    
    tag "$sample_id"
    container "community.wave.seqera.io/library/bcftools:1.21--4335bec1d7b44d11"
    publishDir params.bedtools_output_dir, mode: 'copy', overwrite: true
    
    input:
    tuple val(sample_id), path(eremid_vcf), path(internal_vcf)

    output:
    tuple val(sample_id), 
          path("${sample_id}.eremid.hifi.cnv.bed"), 
          path("${sample_id}.internal.hifi.cnv.bed"), 
          emit: cnv_beds

    script: // use bcftools to make a CNV bed file from the VCF
    """
    bcftools view -f PASS ${eremid_vcf} | bcftools query -f '%CHROM\t%POS\t%INFO/END' - | grep -v _ | grep -v GL > ${sample_id}.eremid.hifi.cnv.bed
    bcftools view -f PASS ${internal_vcf} | bcftools query -f '%CHROM\t%POS\t%INFO/END' - | grep -v _ | grep -v GL > ${sample_id}.internal.hifi.cnv.bed
    """

}