###MToolBox (v.1.0)
###MTOOLBOX
MToolBox is a highly automated bioinformatics pipeline to reconstruct and analyze human mitochondrial DNA from high throughput sequencing data. MToolBox includes an updated computational strategy to assemble mitochondrial genomes from Whole Exome and/or Genome Sequencing (PMID: 22669646) and an improved fragment-classify tool (PMID:22139932) for haplogroup assignment, functional and prioritization analysis of mitochondrial variants. MToolBox provides pathogenicity scores, profiles of genome variability and disease-associations for mitochondrial variants. MToolBox provides also a Variant Call Format file (version 4.0) featuring, for the first time, allele-specific heteroplasmy.  
  
The MToolBox pipeline includes:

- an extended version of a previously published computational strategy for mtDNA genome assembly (PMID: 22669646). The pipeline has been integrated with the detection of insertions and deletions (indels) and the assessment of the heteroplasmic fraction (HF) and related confidence interval (CI) for each mt variant. HF and CI are integrated as genotype specific meta-information in a Variant Call Format (version 4.0) file;
- the mt-classifier tool, for haplogroup prediction, mt variant functional annotation and prioritization.

####CHANGELOG - 11 December 2016

`GenomeAnalysisTK.jar` has been removed from the `MToolBox/ext_tools` directory.  
Users that would like to run GATK IndelRealigner are now asked to download a newer version of GATK and place it in the
`MToolBox/ext_tools` folder:

```
mv GenomeAnalysisTK.jar /path/to/MToolBox/MToolBox/ext_tools/
```

####For all changelogs, please visit [changelogs](https://github.com/mitoNGS/MToolBox/blob/master/changelog.md). 

####For the previous version of MToolBox (v.0.3.3), please go to https://github.com/mitoNGS/MToolBox/releases. (DEPRECATED)

####For further information, please visit the official [MToolBox Wiki](https://github.com/mitoNGS/MToolBox/wiki).

