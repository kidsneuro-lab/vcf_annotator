Feature: Splice distance annotations

  Scenario Outline: Annotate splice distances for a single variant
    Given a transcript index
    And a variant at "<chrom>" <pos> with ref "<ref>" and alt "<alt>"
    When I annotate splice distances with tag "sj" and include_mane
    Then the SJ field should match the expected rows built from "<allele>" "<transcript_id>" "<gene>" "<variant_type>" "<donor_exon_rank>" "<donor_region>" "<donor_phase>" "<donor_distance>" "<acceptor_region>" "<acceptor_phase>" <mane_flag>
    And the TSV SJ should equal the INFO SJ

    Examples:
      | chrom | pos      | ref | alt | allele | transcript_id   | gene   | variant_type | donor_exon_rank | donor_region | donor_phase | donor_distance | acceptor_region | acceptor_phase | mane_flag |
      | chr12 | 48006054 | C   | A   | A      | XM_017018828.1  | COL2A1 | snp          | 1               | intron       | 1           | -325           | intron          | 1              | 0        |