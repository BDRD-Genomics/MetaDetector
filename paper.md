---
title: 'MetaDetector: Tactical Genomic Biosurveillance in Air-Gapped Environments'
tags:
  - bioinformatics
  - metagenomics
  - taxonomy
authors:
- family-names: "Rice"
  given-names: "Gregory K"
  orcid: "https://orcid.org/0000-0001-8509-6379"
  affiliation: "1, 2"
- family-names: "Paskey"
  given-names: "Adrian C"
  orcid: "https://orcid.org/0000-0003-4575-3092"
- family-names: "Thomas"
  given-names: "Quinn K"
  orcid: "https://orcid.org/0009-0008-2920-4484"
  affiliation: "1, 2"
- family-names: "Long"
  given-names: "Kyle A."
  orcid: ""
  affiliation: "1, 2"
- family-names: "Cicalo"
  given-names: "Anthony"
  orcid: ""
  affiliation: "1, 2"
- family-names: "Cer"
  given-names: "Regina Z."
  orcid: "https://orcid.org/0000-0002-2395-980X"
  affiliation: "1"
- family-names: "Bishop-Lilly"
  given-names: "Kimberly A."
  orcid: "https://orcid.org/0000-0002-5744-8493"
  affiliation: "1"
affiliations:
 - name: Genomics and Bioinformatics Department, Biological Defense Research Directorate, Naval Medical Research Command-Frederick, United States
   index: 1
 - name:  Leidos, USA
   index: 2
 - name:  Defense Threat Reduction Agency, USA
   index: 3
date: 3 September 2026
bibliography: paper.bib
---

# Summary

MetaDetector is a versatile bioinformatic pipeline developed for complementary read- and contig-based classification of microbial sequences within complex metagenomic samples. Orchestrated by Slurm workload manager (1), the pipeline provides an end-to-end workflow that accepts both short-read (Illumina) and long-read (Oxford Nanopore, PacBio) sequencing data. It automates quality control (FastQC (2), BBDuk (3), fastp (4)), optional host sequence removal (BBMap (3), minimap2 (5)), de novo assembly (metaSPAdes (6), SPAdes (7), Unicycler (8), Dragonflye (9)), and taxonomic classification of both reads and contigs using DIAMOND (10) and MEGAN (11). The primary outputs, including DIAMOND alignment archive (.daa) files and simplified tabular summary files, provide a flexible and reproducible framework for processing large metagenomic datasets generated during biosurveillance activities.

# Statement of need

The development of MetaDetector was driven by the need for an efficient and modular pipeline capable of handling diverse metagenomic datasets in secure, high-performance computing (HPC) environments with limited or no internet connectivity. Despite their capabilities, many contemporary pipelines are unsuitable for deployment in air-gapped settings. MetaDetector addresses this by packaging open-source tools into portable containers (Docker/Singularity) optimized for offline execution. Its modular architecture allows users to select appropriate assembly strategies for their data, while its use of both read- and contig-based classification provides a more comprehensive view of a sample's microbial content. MetaDetector reduces the manual effort required to process complex sequencing datasets and generate analysis-ready output. 


# State of the field  

MetaDetector performs as well as existing pipelines for metagenomic analysis, building the capacity for high performance in air gapped settings. To evaluate MetaDetector’s capacity to characterize challenging viral genomes using both short-read and hybrid sequencing data, we benchmarked it against four established tools: geNomad v1.8.0 (30), Chan Zuckerberg ID (CZID v6.0 (31, 32), Kraken2 v2.1.3 (33), and Mash v2.3 (34).
Each tool selected for this benchmarking exercise represents a distinct methodological framework, and was executed using default parameters: 
  - geNomad uses a predefined marker dataset to identify viruses and plasmids.
  - CZID performs de novo short-read assembly via SPAdes (16), then classifies contigs via BLAST (35). 
  - Kraken2 uses k-mer matching approach for classification of the lowest common ancestor (LCA).
  - Mash screen estimates hash containment against RefSeq genome sketches.

Two datasets from NCBI’s Sequence Read Archive (SRA) were selected for benchmarking: 
  1. SRR10179613, Illumina short reads generated from a fruit bat body swab, selected to test sensitivity in detecting a low-abundance, recently discovered virus (Dawn bat paramyxovirus, DbPV) against a background of host and bacterial reads.  
  2. PRJNA587334, matched Oxford Nanopore MinION and Illumina MiSeq read sets generated from an mpox virus (MPXV) clinical isolate (29), selected to assess hybrid long/short read assembly and classification. Because CZID does not currently support hybrid assembly, only short reads were evaluated for that platform.

## Benchmarking Results
`Table 1` summarizes the degree to which all evaluated tools identified the target species or a closely related near-neighbor. 
- In the environmental bat-swab sample (Table 1A), MetaDetector classified the highest number of short reads (63 reads) and assembled a corresponding contig. geNomad resolved DbPV only to the family level (Paramyxoviridae), whereas MetaDetector, CZID, and Kraken2 achieved species-level resolution by identifying closely-related near-neighbor lineages.
- All evaluated tools successfully classified Mpox virus in the hybrid assembly dataset (Table 1B). CZID and MetaDetector identified the highest number of MPXV reads, with 54,456 and 49,047 reads respectively. MetaDetector, CZID, and Kraken2 each produced comparable numbers of MPXV contigs (12, 13, and 9 contigs respectively).
- These results demonstrate that MetaDetector’s sensitivity and taxonomic specificity are on par with existing state-of-the-art tools, while providing an offline-capable and customizable pipeline tailored for secure and hybrid metagenomic workflows.

<i> **Table 1.** Summary of benchmarking results via MetaDetector, geNomad, CZID, kraken2, and mash.</i> \
1A. Detection of DbPV
| Tool | Specificity of DbPV detection | Specific assignment | # reads classified to lowest assignment | # contigs classified to lowest assignment |
| --- | --- | --- | --- | --- |
| MetaDetector | moderate | Mumps/bat mumps orthorubulavirus (TaxIDs: 2560602; 2560195) | 63 | 1 |
| geNomad | low | Paramyxoviridae (TaxID: 11158) | n/a | n/a |
| CZID | moderate | Bat mumps orthorubulavirus (TaxID: 2560340) | 16 | 0 |
| kraken2 | moderate | Bat Paramyxovirus Epo_spe/AR1/DRC/2009 (NC_038271) | 2 | 1 |
| mash | moderate | Bat Paramyxovirus Epo_spe/AR1/DRC/2009n(NC_038271) | n/a | n/a |

1B. Handling of hybrid mpox virus data
| Tool | Specificity of MPXV detection | Specific assignment | # reads classified as MPXV | # contigs classified as MPXV |
| --- | --- | --- | --- | --- |
| MetaDetector | high | MPXV (TaxID: 10244) | 49047 | 12 |
| geNomad | low | Poxviridae (TaxID: 10240) | n/a | n/a |
| CZID | high | MPXV (TaxID: 10244) | 54456 | 13 |
| kraken2 | high | MPXV (TaxID: 10244) | 35595 | 9 |
| mash | high | MPXV (NC_063383) | n/a | n/a |

# Software design
The architecture of MetaDetector was designed to address three computational challenges inherent to unbiased biosurveillance and metagenomic characterization: handling large, high-throughput datasets (>10GB) with low viral target abundance (as low as ~0.0003%), integrating multi-platform sequencing input, and maintaining high performance within an air-gapped high-performance computing (HPC) environment. 
  - The Slurm workload manager was chosen over other frameworks to enable automated, multi-thread batch execution on on-premises HPC server clusters without external dependencies.  
  - Dual-platform input and preprocessing tools BBDuk/BBMap for Illumina and fastp/minimap2 for ONT/PacBio reads balance computational throughput with sequencing technology constraints to allow for automated trimming, host depletion, and rRNA/contaminant removal in single execution steps.
  - Support for multiple assemblers (metaSPAdes, SPAdes, Unicycler, Dragonflye), coupled with post-assembly read mapping prioritizes contig accuracy over raw processing speed. Trimmed reads are mapped back to contigs using BBMap/minimap2 to verify assembly validity prior to classification.
  - Taxonomic classification via Diamond BLASTX (NCBI NR) and MegaBLAST (NCBI core_nt) binned via MEGAN LCA prevents misclassification of divergent sequences by incorporating protein analysis while also suppressing false-positive calls via MEGAN’s weighted LCA algorithms.

# Research impact statement

In biosurveillance and military medical operations, sequencing often occurs in settings where cloud-only platforms are inaccessible. MetaDetector’s design prioritizes parameter modularity, enabling researchers to tune quality thresholds and assembler parameters depending on sample complexity, host genomic contribution, and compute limits, ultimately generating standardized MEGAN (.daa) and tabular outputs compatible with downstream visualization pipelines such as Pavian. MetaDetector has been adopted for routine use in five peer reviewed publications describing biosurveillance efforts (12-16). MetaDetector was implemented using Docker (45) to provide reach back support and train 19 personnel at a U.S. Department of Defense laboratory overseas. Trainees were able to individually analyze a sample from long and short read data through identifying the organism(s) of interest and performed advanced characterization, including phylogenetic analysis, of the target organism(s).

# Acknowledgements and Disclaimers
This work was supported by Navy WUN A1417 and Global Emerging Infections Surveillance (GEIS) Branch ProMIS ID P0054_23_NM to KAB-L.
The views expressed in this article are those of the authors and do not necessarily reflect the official policy or position of the Department of the Navy, Department of Defense, nor the U.S. Government. Some authors are employees of the U.S. Government. This work was prepared as part of their official duties. Title 17 U.S.C. §105 provides that “Copyright protection under this title is not available for any work of the United States Government”. Title 17 U.S.C. §101 defines a U.S. Government work as a work prepared by a military service member or employee of the U.S. Government as part of that person’s official duties.

# AI usage disclosure

Generative AI was not used in any part of the software creation or documentation, which was designed and prepared by human authors. GenAI.mil was used to proofread sections of the written documentation and any edits resulting from those recommendations were vetted and implemented by human authors.

# Availability

MetaDetector is available at https://github.com/BDRD-Genomics/MetaDetector. The repository includes comprehensive documentation, example data, and a minimal test profile to verify installation. To demonstrate the pipeline’s performance in a real-world scenario, we have included a case study in the documentation showing the successful detection of known viruses from publicly available datasets (SRR10179613, PRJNA587334) using both long and short reads. 

# References
1.	Yoo AB, Jette MA, Grondona M. Slurm: Simple linux utility for resource management, p 44-60. In (ed),  Springer, 
2.	Andrews S. 2017. FastQC: a quality control tool for high throughput sequence data. 2010.
3.	Bushnell B. 2014. BBMap: a fast, accurate, splice-aware aligner. Lawrence Berkeley National Lab (LBNL).
4.	Chen S. 2023. Ultrafast one‐pass FASTQ data preprocessing, quality control, and deduplication using fastp. Imeta 2:e107.
5.	Li H. 2018. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics 34:3094-3100.
6.	Nurk S, Meleshko D, Korobeynikov A, Pevzner PA. 2017. metaSPAdes: a new versatile metagenomic assembler. Genome Res 27:824-834.
7.	Bankevich A, Nurk S, Antipov D, Gurevich AA, Dvorkin M, Kulikov AS, Lesin VM, Nikolenko SI, Pham S, Prjibelski AD. 2012. SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. Journal of computational biology 19:455-477.
8.	Wick RR, Judd LM, Gorrie CL, Holt KE. 2017. Unicycler: resolving bacterial genome assemblies from short and long sequencing reads. PLoS computational biology 13:e1005595.
9.	Petit III RA. 2024. Dragonflye: assemble bacterial isolate genomes from nanopore reads. Github2021.
10.	Buchfink B, Reuter K, Drost HG. 2021. Sensitive protein alignments at tree-of-life scale using DIAMOND. Nat Methods 18:366-368.
11.	Huson DH, Auch AF, Qi J, Schuster SC. 2007. MEGAN analysis of metagenomic data. Genome research 17:377-386.
12.	Adhikari BN, Paskey AC, Frey KG, Bennett AJ, Long KA, Kuhn JH, Hamilton T, Glang L, Cer RZ, Goldberg TL, Bishop-Lilly KA. 2024. Virome profiling of fig wasps (Ceratosolen spp.) reveals virus diversity spanning four realms. Virology 591:109992.
13.	Bennett AJ, Paskey AC, Kuhn JH, Bishop-Lilly KA, Goldberg TL. 2020. Diversity, Transmission, and Cophylogeny of Ledanteviruses (Rhabdoviridae: Ledantevirus) and Nycteribiid Bat Flies Parasitizing Angolan Soft-Furred Fruit Bats in Bundibugyo District, Uganda. Microorganisms 8.
14.	Paskey AC, Lim XF, Ng JHJ, Rice GK, Chia WN, Philipson CW, Foo R, Cer RZ, Long KA, Lueder MR, Glang L, Frey KG, Hamilton T, Mendenhall IH, Smith GJ, Anderson DE, Wang LF, Bishop-Lilly KA. 2023. Genomic Characterization of a Relative of Mumps Virus in Lesser Dawn Bats of Southeast Asia. Viruses 15.
15.	Paskey AC, Ng JHJ, Rice GK, Chia WN, Philipson CW, Foo RJH, Cer RZ, Long KA, Lueder MR, Frey KG, Hamilton T, Mendenhall IH, Smith GJ, Wang LF, Bishop-Lilly KA. 2020. The temporal RNA virome patterns of a lesser dawn bat (Eonycteris spelaea) colony revealed by deep sequencing. Virus Evol 6:veaa017.
16.	Bennett AJ, Paskey AC, Ebinger A, Pfaff F, Priemer G, Hoper D, Breithaupt A, Heuser E, Ulrich RG, Kuhn JH, Bishop-Lilly KA, Beer M, Goldberg TL. 2020. Relatives of rubella virus in diverse mammals. Nature 586:424-428.
