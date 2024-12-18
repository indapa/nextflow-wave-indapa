process pb_cpg_tools {
  debug true
  """
  aligned_bam_to_cpg_scores --version
  """
}

process trgt {
  debug true

  """
  trgt --version 
  """
}

process pbsv {
  debug true

  """
  pbsv --version
  """

}

process hificnv {
  debug true
  """
  hificnv --version
  """
}

process pbmm2 {
  debug true
  """
  pbmm2 --version
  """
}