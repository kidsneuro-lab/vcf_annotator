Feature: Splice distance annotations

  Scenario Outline: Annotate splice distances for a single variant
    Given a transcript index
    And a variant at "<chrom>" <pos> with ref "<ref>" and alt "<alt>"
    When I annotate splice distances with tag "sj" and include_mane
    Then the SJ field should match the expected rows built from "<allele>" "<transcript_id>" "<gene>" "<variant_type>" "<donor_distance>" "<donor_region_type>" "<donor_region_number>" "<acceptor_distance>" "<acceptor_region_type>" "<acceptor_region_number>" <mane_flag>
    And the TSV SJ should equal the INFO SJ

    Examples:
      | chrom | pos      | ref | alt | allele | transcript_id  | gene   | variant_type | donor_distance | donor_region_type | donor_region_number | acceptor_distance | acceptor_region_type | acceptor_region_number | mane_flag |
      | chr12 | 48006054 | C   | A   | A      | XM_017018828.1 | COL2A1 | snp          | 1              | intron            | 1                   | -325              | intron               | 1                      | 0         |
      | chr12 | 48006055 | C   | A   | A      | XM_017018828.1 | COL2A1 | snp          | -1             | exon              | 1                   | NULL              | NULL                 | NULL                   | 0         |