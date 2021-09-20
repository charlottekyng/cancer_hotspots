# Cancer Hotspots

### A set of cancer hotspots for the analysis of somatic mutations

The hotspots in this resource derive from [www.cancerhotspots.org](https://www.cancerhotspots.org/#/download) and [www.3dhotspots.org](https://www.3dhotspots.org/#/download), plus four additional non‑coding hotspots from [Weinhold *et al*. 2014](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4217527/).

Original hotspots were aggregated into `resources/hg19.hotspots_changv2_gao_nc.vcf.gz` (and the liftover `resources/hg38.hotspots_changv2_gao_nc.snpeff.tab.gz`). These VCFs were used in our [usb-modules-v2 pipeline](https://github.com/charlottekyng/usb-modules-v2) through `SnpSift annotate`. We realised this strategy was not ideal because `SnpSift annotate` would tag any variant that falls within a region in the hotspot VCF file, resulting in some false positives (especially due to large hotspot indels that cover wider genomic ranges), or false negatives (variants that do not overlap with the genomic regions in the hotspot VCF but affect a hotspot amino acid).

We therefore sought to make the hotspot annotation more stringent by accounting for different mutation types (SNVs and indels), and genomic effects (missense, frameshifts, splice variants etc.), and by referring to the hotspot amino acid position for SNV and indel variants, and codon positions for splice variants, instead of relying on genomic coordinates.

### Overview of the creation of the new hotspot resources
#### 1. Hotspot annotation based on the current `SnpEff` used in our pipelines
In order to base our lookup on amino acid or codon positions, it is necessary to annotate the hotspot VCFs in the same way we annotate our somatic variants, i.e. using `SnpEff v.4.3t` with the `-canon` option.
```
resources/hg19.hotspots_changv2_gao_nc.snpeff.tab.gz
resources/hg38.hotspots_changv2_gao_nc.snpeff.tab.gz
```
We also annotated the VCFs with `dbNSFP v.4.1a` and extracted the fields:
```
dbNSFP_gnomAD_exomes_AF
dbNSFP_ExAC_nonTCGA_AF
dbNSFP_ExAC_nonTCGA_Adj_AF
dbNSFP_1000Gp3_AF
```

The above files were manually curated in a spreadsheet program, resulting in
```
resources/hg19.hotspots_changv2_gao_nc.snpeff.xlsx
resources/hg38.hotspots_changv2_gao_nc.snpeff.xlsx
```
Edge cases included updating some gene names and gene effects.

#### 2. Single amino acid hotspots
These hotspots derive from SNV hotspots and comprise the vast majority of our hotspots. For each gene we extracted the amino acid reference and position and keep only the unique entries.
We omitted hotspots with `AF_all_max > 0.01` (which is the `max` of the population frequencies from the dbNSFP fileds listed above). Note that we did not omit the amino acid if there were other hotspots affecting it. In most cases, especially those with very high `AF_all_max`, the filtering resulted in discarding the amino acid.

Final hotspots:
```
hotspots_single_aa_change.hg19.txt
hotspots_single_aa_change.hg38.txt
```

Lookup rules:
1. Variant has to be annotated as `missense_variant`, `stop_gained`, `initiator_codon_variant`, or `start_lost`.
2. Variant's gene name, reference amino acid code, and its position have to match the hotspot's entry.


#### 3. Splice hotspots
A modest fraction of hotspots affect splicing sites. We extracted the codon position from the `SnpEff` annotations `splice_acceptor_variant` and `splice_donor_variant`.

Final hotspots:
```
hotspots_splice.hg19.txt
hotspots_splice.hg38.txt
```

Lookup rules:
1. Variant has to be annotated as `splice_acceptor_variant` or `splice_donor_variant`.
2. Variant's gene name and codon position have to match the hotspot's entry.

#### 4. Indel hotspots
Indel hotspots were extracted from `resources/cancerhotspots.org_hotspots_v2.xls`, sheet `INDEL-hotspots`. For every gene we checked whether the amino acid position corresponds to the new `SnpEff -canon` annotation, and corrected it accordingly.

The amino acid positions in the indel entries are the same between `hg19` and `hg38`.

Final hotspots:
```
hotspots_indel.txt
```

Lookup rules:
1. Variant's annotation has to contain the string `inframe`.
2. Variant's gene name and amino acid start and end positions have to match the hotspot's entry.


#### 5. Non-coding hotspots
Non-coding hotspots remain in their original VCF format because they affect specific genomic regions.

Final hotspots:
```
hotspots_non_coding.hg19.vcf
hotspots_non_coding.hg38.vcf
```

Lookup rules:
1. Variant has to occur in the same region as the non-coding hotspot.

### Caveats
Some canonical transcripts used by `SnpEff v.4.3t` differ between `hg19` and `hg38`. A noteworthy case is `CDKN2A` whose canonical transcript in `hg19` corresponds to an alternative reading frame. A possible workaround would be to use `VEP` instead of `SnpEff` (whose canonical transcripts might be biologically more relevant), or to force `SnpEff` to use `VEP` canonical transcripts with the option `-canonList`.
