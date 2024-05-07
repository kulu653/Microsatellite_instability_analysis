
# Microsatellite_instability_analysis

Results and methodology of this analysis is published in: 
		 	 	 		
1. Shrestha KS, Aska EM, Tuominen MM, Kauppi L. Tissue-specific reduction in MLH1 expression induces microsatellite instability in intestine of Mlh1+/- mice. DNA Repair (Amst). 2021 Oct;106:103178. doi: 10.1016/j.dnarep.2021.103178. Epub 2021 Jul 9. PMID: 34311271
 							
2. Shrestha KS, Tuominen MM, Kauppi L. Mlh1 heterozygosity and promoter methylation associates with microsatellite instability in mouse sperm. Mutagenesis. 2021 Mar 19:geab010. doi: 10.1093/mutage/geab010. Epub ahead of print. PMID: 33740045 
 							
## Introduction to microsatellite and microsatellite instability (MSI)

Microsatellites are defined as tandem repeat DNA sequence up to a total of 100 nucleotides in length, consisting of 1-6 bp repeat units. These are among the most abundant repetitive genomic entity. Microsatellite instability (MSI) is a hallmark of Lynch syndrome (LS) associated cancer, these cancer types show missmatch repair defects. Selected panel of microsatellite are tested for MSI to diagnose LS- cancers. MSI, as reported by the microsatellite markers, confer the global genome instability and LS-cancers. Microsatellite length analysis serves a multifaceted role beyond cancer diagnosis, extending to mutation detection, linkage analysis, quantitative trait loci (QTL) mapping, genetic diversity assessment, pedigree analysis, and heterozygosity detection.

## Novelty of my analysis: 

There are few commercial softwares for MSI analysis (also known as fragment analysis), however, none of them have high-throughput applicability. Also there is no analysis pipeline to perform end-to-end MSI analysis. This analysis utilizes an R package for fragment analysis (Fragman) and performs end-to-end MSI analysis using FSA extension files (FASTA-type file) containing DNA fragment intensities read by capillary electrophoresis machine as input, and produces report (as .CSV) files reporting MSI in single DNA molecule level, along with visualizations for the same. 


## MSI data analysis pipeline:

This MSI analysis tests MSI at 3 microsatellite loci (A27 (green panel), A33 (blue panel) and D14Mit15 (yellow panel)) at single DNA molecule level. The data analysis pipeline is as follows, each section of analysis is also annotated in the R codes provided. 

1. Data import: import .fsa files to the R environment
2. Data analysis: single DNA molecule MSI analysis
3. Data export: final report - summary MSI score table of the three microsatellites and ladder scores.

