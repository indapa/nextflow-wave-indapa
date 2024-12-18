include { hello } from './modules/foo'
include { pb_cpg_tools; trgt; hificnv; pbsv; pbmm2 } from './modules/pbtools'

workflow {
  hello()
  pb_cpg_tools()
  trgt()
  hificnv()
  pbsv()
  pbmm2()
}
