# Test data
## Sample GenePred

### Source
- https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/GCF_000001405.40-RS_2025_08/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz

### 

**Convert annotation GTF to genepred format**
```bash
gtfToGenePred -ignoreGroupsWithoutExons -genePredExt -includeVersion GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.genepred
```

**Obtain subset of genepred for EMD, HBB, COL2A1 and DYSF**
```bash
awk -F"\t" '$12~/^EMD$|^HBB$|^COL2A1$|^DYSF$/' GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.genepred > GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.sample.genepred

rm GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.genepred
```
## Sample MANE list

### Source
- https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/MANE.GRCh38.v1.5.summary.txt.gz


```bash
gunzip -c MANE.GRCh38.v1.5.summary.txt.gz | awk -F"\t" 'NR==1 || $4=="COL2A1" || $4=="DYSF" || $4=="EMD" || $3=="HBB"' > MANE.GRCh38.v1.5.summary.sample.txt

rm MANE.GRCh38.v1.5.summary.txt.gz
```