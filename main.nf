#!/usr/local/bin/nextflow

nextflow.enable.dsl=2

include { pbmm2_align; cpg_pileup } from './modules/pbtools'
paramsFromSchema(params, 'nextflow_schema.json')


// Replace the input_bams channel definition with this:
def checkSamplesheet(samplesheet_file) {
    if (!file(samplesheet_file).exists()) {
        exit 1, "Samplesheet file not found: ${samplesheet_file}"
    }
    return file(samplesheet_file)
}

ss_status = checkSamplesheet(params.samplesheet)

// Create channels for input BAM files
Channel.fromPath(params.samplesheet)
    .splitCsv(header: true)
    .map { row -> 
        def sample_id = row.sample_id
        def bam_file = file(row.bam_file)
        if (!bam_file.exists()) {
            error "BAM file not found: ${bam_file}"
        }
        return tuple(sample_id, bam_file)
    }
    .set { input_bams_ch }

workflow {

input_bams_ch.view { sample_id, bam -> "Sample ID: $sample_id, BAM: $bam" }


pbmm2_align(
        file(params.reference),
        //input_bams,
        input_bams_ch,
        params.cpu,
        params.sort_threads
    )
    // print the pbmm2_align channel

    pbmm2_align.out.aligned_bam.view { "Aligned BAM: $it" }

 cpg_pileup(
        pbmm2_align.out.aligned_bam,
        file(params.reference),
        file(params.reference_index),
        file(params.cpgmodel)
    )

}

