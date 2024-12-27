#!/usr/local/bin/nextflow

nextflow.enable.dsl=2

include { pbmm2_align; cpg_pileup; hificnv; trgt; extractRegions; pb_discover; pb_call } from './modules/pbtools'
include { mosdepth } from './modules/mosdepth'


def required_params = ['reference', 'samplesheet', 'cpgmodel', 'karyotype', 'sv_output_dir']
for (param in required_params) {
    if (!params[param]) {
        error "Parameter '$param' is required!"
    }
}


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

/* read alignment */
pbmm2_align(
        file(params.reference),
        //input_bams,
        input_bams_ch,
        params.cpu,
        params.sort_threads
    )
    // print the pbmm2_align channel

    pbmm2_align.out.aligned_bam.view { "Aligned BAM: $it" }

 /* cpg calling */
 cpg_pileup(
        pbmm2_align.out.aligned_bam,
        file(params.reference),
        file(params.reference_index),
        file(params.cpgmodel)
    )

/* read depth analysis */
mosdepth(pbmm2_align.out.aligned_bam)


/* aligned bam channel used for cnv, tandem repeat and sv analysis */
bam_bai_ch = pbmm2_align.out.aligned_bam.map { bam, bai -> 
    def sample_id = bam.baseName.replaceFirst(/\..*$/, '')
    tuple(sample_id, bam, bai)
}

bam_bai_ch.view { sample_id, bam, bai -> "Input: $sample_id, BAM: $bam, BAI: $bai" }

/* cnv analysis */
hificnv(
        bam_bai_ch,
        file(params.reference),
        file(params.reference_index),
        file(params.exclude_bed),
        file(params.expected_bed),
        params.cpu
    )


/*tandem repeat analysis */
trgt(   bam_bai_ch,
        file(params.reference),
        file(params.reference_index),
        file(params.trgt_repeats),
        params.karyotype,
        params.cpu
    )




/* SV analysis
        regions
        discover
        call 
*/
  
    regions_ch = extractRegions(bam_bai_ch)
        .flatMap { sample_id, regions -> 
            regions.split('\n').collect { region -> tuple(sample_id, region.trim()) }
        }

    // Combine regions with corresponding BAM and BAI files
    discover_input_ch = regions_ch.combine(bam_bai_ch, by: 0)
    
    discover_input_ch.view { sample_id, region, bam, bai -> "Input: $sample_id, Region: $region, BAM: $bam, BAI: $bai" }


    // Run pb_discover process
    pb_discover_results = pb_discover(discover_input_ch, params.trf_bed)

    // Group svsig files by sample_id
    svsig_files_by_sample = pb_discover_results.groupTuple()

    // Create a channel for the reference genome
    reference_ch = channel.fromPath(params.reference)

    // Run pb_call process
    pb_call(svsig_files_by_sample, reference_ch)


}


