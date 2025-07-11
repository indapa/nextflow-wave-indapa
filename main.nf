#!/usr/local/bin/nextflow

nextflow.enable.dsl=2

include { pbmm2_align; cpg_pileup; hificnv; trgt; pb_discover; pb_call; hiphase } from './modules/pbtools'
include { mosdepth } from './modules/mosdepth'
include { deepvariant; deepvariant_chr20 } from './modules/deepvariant'
include {fibertools_extract; fibertools_extract_all; fibertools_extract_all_to_parquet; m6a_extract_parquet} from './modules/fiberseq'
include {annotate_vep_no_phased; annotate_vep} from './modules/ensemblvep'
include{bam_stats} from './modules/samtools'




def required_params = ['reference', 'samplesheet',  'karyotype']
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

    Channel.fromPath(params.samplesheet)
    .splitCsv(header: true)
    .map { row -> row.sample_id }
    .set { sample_ids_ch }



def REGIONS = [
    'chr1',
    'chr2',
    'chr3',
    'chr4',
    'chr5',
    'chr6',
    'chr7',
    'chr8',
    'chr9',
    'chr10',
    'chr11',
    'chr12',
    'chr13',
    'chr14',
    'chr15',
    'chr16',
    'chr17',
    'chr18',
    'chr19',
    'chr20',
    'chr21',
    'chr22',
    'chrX',
    'chrY'
   
]

// Create a channel from the fixed regions
Channel
    .fromList(REGIONS)
    .set { regions_ch }

workflow {
    
    //input_bams_ch.view { sample_id, bam -> "Sample ID: $sample_id, BAM: $bam" }

    //regions_ch.view { region -> "BED region: $region" } 

/* read alignment */
    pbmm2_align(
        file(params.reference),
        //input_bams,
        input_bams_ch,
        params.cpu,
        params.sort_threads
    )
    // print the pbmm2_align channel

    //pbmm2_align.out.aligned_bam.view { "Aligned BAM: $it" }

 /* cpg calling  */
    cpg_pileup(
        pbmm2_align.out.aligned_bam
        
    )

/* read depth analysis  */
    mosdepth(pbmm2_align.out.aligned_bam)





/* aligned bam channel used for cnv, tandem repeat and sv analysis */
    bam_bai_ch = pbmm2_align.out.aligned_bam.map { bam, bai -> 
    def sample_id = bam.baseName.replaceFirst(/\..*$/, '')
    tuple(sample_id, bam, bai)
}

    //bam_bai_ch.view { sample_id, bam, bai -> "Input: $sample_id, BAM: $bam, BAI: $bai" }

/* cnv analysis */
    hificnv(
        bam_bai_ch,
        file(params.reference),
        file(params.reference_index),
        file(params.exclude_bed),
        file(params.expected_bed),
        params.cpu
    )


/*tandem repeat analysis  */
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
  
       bam_bai_ch
        .combine(regions_ch)
        .map { sample_id, bam, bai, region -> 
            // For each region, create the tuple: [sample_id, region, bam, bai]
            [sample_id, region, bam, bai]
        }
        .set { bam_regions_ch }

        
        //bam_regions_ch.view { sample_id, region, bam, bai -> "sample: $sample_id, bed: $region, BAM: $bam, BAI: $bai"}

        pb_discover_results=pb_discover(bam_regions_ch, params.trf_bed )

         // Group svsig files by sample_id
        svsig_files_by_sample = pb_discover_results.groupTuple()

        

        // Run pb_call process
        pb_call(svsig_files_by_sample, file(params.reference))

        

        //deepvariant

        deepvariant(params.reference, params.reference_index, bam_bai_ch, params.deepvariant_threads )
        //deepvariant_chr20(params.reference, params.reference_index, bam_bai_ch, params.deepvariant_threads )

    


        //hiphase

        /*hiphase( skip for now because the SM headers and VCF headers are not matching and it throws an error
            deepvariant.out.vcf,
            deepvariant.out.vcf_tbi,
            pb_call.out.pb_call,
            trgt.out.repeat_vcf,
            pbmm2_align.out.aligned_bam,
            file(params.reference)
        )*/


        //fibertools extract generates bed6 

        fibertools_extract_all (
            bam_bai_ch
        )

        // convert bed6 to parquet format
        fibertools_extract_all_to_parquet(
            fibertools_extract_all.out.ft_all_bed
        )

        //m6a extract parquet data to CSV
        m6a_extract_parquet(
            fibertools_extract_all_to_parquet.out.ft_all_parquet
        )


        //vep annotate

        annotate_vep_no_phased(
            deepvariant.out.vcf,
            deepvariant.out.vcf_tbi,
            sample_ids_ch,
            file(params.pigeon_gtf_bgzip),
            file(params.pigeon_gtf_tbi),
            file(params.reference)
        )


        bam_stats(
            pbmm2_align.out.aligned_bam
        )

        

}


