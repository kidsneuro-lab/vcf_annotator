



# Test data

## Sample GenePred

```bash
awk -F"\t" '$12=="COL2A1" && ($1=="NM_033150.3" || $1=="NM_001844.5" || $1=="XM_017018828.1")' GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.genepred > GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.sample.genepred
awk -F"\t" '$12=="EMD"' GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.genepred >> GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.sample.genepred
awk -F"\t" '$12=="HBB"' GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.genepred >> GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.sample.genepred
awk -F"\t" '$12=="DYSF" && ($1=="NM_003494.4" || $1=="NM_001130987.2" || $1=="NM_001130986.2")' GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.genepred >> GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.sample.genepred
```

## Sample MANE list
```bash
awk -F"\t" 'NR==1 || $4=="COL2A1" || $4=="DYSF" || $4=="EMD" || $3=="HBB"' MANE.GRCh38.v1.4.summary.txt > MANE.GRCh38.v1.4.summary.sample.txt
```