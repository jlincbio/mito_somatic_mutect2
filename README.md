# mito_somatic_mutect2
## Mitochondrial somatic mutation detection pipeline by Mutect2 

Scripts used for Ikeda et. al, "Novel immune evasion involving mitochondrial transfer in the tumour microenvironment" submitted for review.
The main R script (__mutect_mito_process.R__) contains the workflow, with comment blocks separating each section, and the accompanying shell script is the qsub script for Mutect2.

* Updated Feb. 16, 2024 to reflect the additional batch of sequencing results for updated manuscript
* Updated Jun. 22, 2024 
  
[![DOI](https://zenodo.org/badge/590244285.svg)](https://zenodo.org/doi/10.5281/zenodo.10685680)


### System Requirements
* Hardware required: a computing cluster or a modern multi-threaded computer/workstation; for a single workstation the recommended specs include an 8-core/16-thread CPU with 64GB or 128GB of memory. We performed the analysis on a computing cluster with SGE for job dispatch, but you can simply execute the shell script by itself without a cluster.
* Operating system: Linux (tested here on a multi-machine cluster running Ubuntu Server 20.04 LTS); the script itself should also run in macOS or Windows Subsystem for Linux (WSL) provided that bash is the shell of choice
* A modern version of [R](https://www.r-project.org/) with [data.table](https://cran.r-project.org/web/packages/data.table/index.html)
* [GATK toolkit](https://gatk.broadinstitute.org/hc/en-us): version 4.1.8.1 as tested
* [EAGLE](https://github.com/tony-kuo/eagle)
* [Haplogrep](https://haplogrep.i-med.ac.at/haplogrep2/index.html): CLI version 2.4.0 as tested

### Analytical Procedure

* Steps highlighted here largely follow __Total mtDNA sequencing for clinical samples__ in the manuscript and are provided here for reference.
* TILs sequenced: matched LCLs/peripheral blood and/or cancer cells from cohort A, FFPE tumour samples from 95 patients from cohort B, 197 patients with NSCLC (cohort C1/2) and matched normal tissues from 45 patients for total mtDNA sequencing.
* Library/Instrumentation: Precision ID mtDNA Whole Genome Panel on IonTorrent Proton

#### Preparation Step

* BaseCaller from Ion Suite for adapter trimming
  ```--trim-qual-cutoff 15 --barcode-filter-minreads 10 --phasing-residual-filter=2.0 --max-phasing-levels 2 --num-unfiltered 1000 --barcode-filter-postpone 1```
* Mapping
  Reference: Revised Cambridge Reference Sequence (rCRS) mitochondrial sequence
  Tool: TMAP 5.12.3-i bam -v -Y -u --prefix-exclude 5–o 2 stage1 map4
  This generates binary alignment maps (BAMs) that are used for the script provided here
  
#### Somatic Mutation Calling
For each sample, variants were called by Mutect2 in the Genome Analysis Toolkit (GATK, version 4.1.8) under the mitochondrial mode with the read filter marked as duplicate disabled. The detected variants were annotated with GRCh38.p14 as reference by SnpEff (version 5.1d)84 and subjected to additional candidate validation by EAGLE (version 1.1.1)85. Allele frequencies computed by EAGLE (--hetbias=0 and –omega=1e-6) were used as the basis for subsequent variant analyses. After pooling the EAGLE results from all samples, we used Z-scores determined from variant allele frequency (VAF) per position as a metric to evaluate potential variants. Alterations with Z-scores > 3, VAF > 0.2, and read depths exceeding 100 were selected as measures to minimise false positives from sequencing errors. In cohorts B and C1/2, variants with VAF > 0.85 were initially considered polymorphisms and excluded from final reporting under the assumption that somatic variants would not exhibit high VAFs because of the presence of a substantial fraction of non-cancerous cells within the tumour tissue in the specimens collected. Haplogrep (version 2.4.0, Kulcyznski classification mode with 17_FU1 phylotree, 10 top hits, extended final output, and heteroplasmy levels set to 0.85)86 was used to characterise potential haplogroups of individual cases and possible unlabelled single nucleotide polymorphisms from our frequency-based selection criteria. Variants labelled as “hotspot” and “local private variant” were regarded as polymorphic variants; filtered variants were visually inspected to exclude probable sequencing errors at the termini of PCR amplicons or within homopolymer sequence stretches. For performance assessment, we compared the presence of true variants, confirmed by a paired analysis of matched tumour and normal tissue samples (n = 45), with variants that were called only from tumour samples. Under these criteria, true variants were called with a false positive rate of 0% and a false negative rate of 12.2%. In cohort A, where a subject presented LCL, PBL, TIL, and/or cancer cell sequencing data, the results were compared with LCLs or PBLs as controls for somatic calling; mutations in TILs or cancer cells were not considered somatic if they existed in the matched LCLs. Under this criterion, 10 mtDNA mutations were confirmed in TILs and cancer cells (Extended Data Table 1 and Supplementary Data 1). In cohorts B and C1/2, we identified 158 variants, including 48 frameshift and 110 substitution variants, from 237 FFPE tumour samples (Supplementary Data 4). For each subject, the overall mtDNA variant status was classified as truncating, missense, silent, tRNA, rRNA, D-loop, or intergenic sites. For survival analysis, subjects were divided into mtDNA mutations (+) if they presented truncating, missense, tRNA, or rRNA variants (melanoma, 32.6%; NSCLC [C1], 51.2%; NSCLC [C2], 35.7%) and mtDNA mutations (−) if they had a D-loop or an intergenic site or had silent or no variants.

### References

1. R Core Team (2024). _R: A Language and Environment for Statistical Computing_. R Foundation for Statistical Computing, Vienna, Austria.
  <https://www.R-project.org/>.
2. Dowle M, Srinivasan A (2023). _data.table: Extension of `data.frame`_. R package version 1.14.8, <https://CRAN.R-project.org/package=data.table>.
3. Van der Auwera GA & O'Connor BD. (2020). Genomics in the Cloud: Using Docker, GATK, and WDL in Terra (1st Edition). O'Reilly Media.
4. Tony Kuo and Martin C Frith and Jun Sese and Paul Horton. EAGLE: Explicit Alternative Genome Likelihood Evaluator. BMC Medical Genomics. 11(Suppl 2):28. https://doi.org/10.1186/s12920-018-0342-1
5. Weissensteiner H, Pacher D, Kloss-Brandstätter A, Forer L, Specht G, Bandelt HJ, Kronenberg F, Salas A, Schönherr S. HaploGrep 2: mitochondrial haplogroup classification in the era of high-throughput sequencing. Nucleic Acids Res. 2016 Jul 8;44(W1):W58-63. doi: 10.1093/nar/gkw233. Epub 2016 Apr 15. PMID: 27084951; PMCID: PMC4987869.