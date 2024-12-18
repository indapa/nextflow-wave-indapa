process samtools {
  container 'community.wave.seqera.io/library/samtools:08a02f71fa4f25f4'
  debug true
  """
  samtools --version-only
  """
}



process tabix {
    container 'community.wave.seqera.io/library/tabix:fb9b456586dfe79d'
    debug true
    """
    tabix 
    """
}

process mosdepth {
  container 'community.wave.seqera.io/library/mosdepth:c3f6e8b37658dce8'
  debug true
  """
  mosdepth -h

  """
}

process bedtools {
  container 'community.wave.seqera.io/library/bedtools:5bd0542220837f43'
  debug true

  """
  bedtools --version
  """
}

workflow {
  samtools()
  //bwToBedGraph()
  mosdepth()
  bedtools()
}