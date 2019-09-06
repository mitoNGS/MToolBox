### MToolBox (v.1)

MToolBox is a highly automated bioinformatics pipeline to reconstruct and analyze human mitochondrial DNA from high throughput sequencing data. MToolBox includes an updated computational strategy to assemble mitochondrial genomes from Whole Exome and/or Genome Sequencing (PMID: 22669646) and an improved fragment-classify tool (PMID:22139932) for haplogroup assignment, functional and prioritization analysis of mitochondrial variants. MToolBox provides pathogenicity scores, profiles of genome variability and disease-associations for mitochondrial variants. MToolBox provides also a Variant Call Format file (version 4.0) featuring, for the first time, allele-specific heteroplasmy.  
  
The MToolBox pipeline includes:

- an extended version of a previously published computational strategy for mtDNA genome assembly (PMID: 22669646). The pipeline has been integrated with the detection of insertions and deletions (indels) and the assessment of the heteroplasmic fraction (HF) and related confidence interval (CI) for each mt variant. HF and CI are integrated as genotype specific meta-information in a Variant Call Format (version 4.0) file;
- the mt-classifier tool, for haplogroup prediction, mt variant functional annotation and prioritization.


The MToolBox development and manteinance relies on the work of undegraduates, PhD and Postdocs, which we would like to thank for their invaluable support:
- [Claudia Calabrese](https://github.com/clody23)
- [Domenico Simone](https://github.com/domenico-simone)
- [Mariangela Diroma](https://github.com/ma-diroma)
- [Roberto Preste](https://github.com/robertopreste)
- [Rosanna Clima](https://github.com/Ros85)
- Mariangela Santorsola
- Ornella Vitale
- Cristiano Gutta' (past)
- Francesco Rubino (past)

We would also like to thank Professor Marcella Attimonelli (University of Bari), who inspired and supervises the MToolBox project.


**As of September 2019**

Update to MToolBox v.1.2

Please check the changelog for more details.

**As of April 2019**


We have accidentally pushed fixes to master that have not been adequately tested. Fixes are the following, pushed from the 4th of March 2019 onwards:

```
4186b5add55e7aa3eb2560565bdc65b62c677f85
075200422ae616c8f9d532756ae5a266e67bf929
46ac15cf32f06a55b67002608e1f2856d94cf697
0d83ee821db1d41422db26c116eb028af8d02bec
03d1c37ad375d3f47fefac5f37b478641d9ef057
```

For those users who update the MToolBox pipeline by doing `git pull` we ask to repeat the `git pull` and then undo local changes by doing `git checkout 8dae1134ae20f0ee3d41296cadb47433ed8118d8`.

This will reset your local MToolBox repo to the MToolBox tag release v.1.1.

We apologize for this inconvenient.

**As of February 2019**

A bug in the phylotree build 17 tree parsing has been fixed, causing replacement of `haplogroups.txt` and `phylotree_r17.pickle` files. The bug was causing wrong predictions mostly for X haplogroups. Please update your MToolBox installation with `git pull`.


**As of April 2018**    
Update to MToolBox v.1.1

Update of the [Sitevar nucleotide variability](http://www.hmtdb.uniba.it/siteVariability). The nucleotide variability of each variant site was estimated on the multi-alignment of 30,860 healthy genomes reported in [HmtDB](http://www.hmtdb.uniba.it/hmdb/). These genomes were derived from primary INSDC databases and individual submissions. Credits to [Rosanna Clima](https://github.com/Ros85) and Ornella Vitale. 

for more details please check the changelog file.

**As of 24 October 2017**    

MToolBox can be now also installed on Mac OS X by specifing `-o` in the `install.sh` command line. For further details please visit:
[MToolBox installation](https://github.com/mitoNGS/MToolBox/wiki/Installation). Credits to [Roberto Preste](https://github.com/robertopreste).


**As of 11 December 2016**

`GenomeAnalysisTK.jar` has been removed from the `MToolBox/ext_tools` directory.  
Users that would like to run GATK IndelRealigner are now asked to download a newer version of GATK and place it in the
`MToolBox/ext_tools` folder:

```
mv GenomeAnalysisTK.jar /path/to/MToolBox/MToolBox/ext_tools/
```

To keep track of the **MToolBox updates** please visit [changelogs](https://github.com/mitoNGS/MToolBox/blob/master/changelog.md) and take a look at the MToolBox [latest releases](https://github.com/mitoNGS/MToolBox/releases).


### For further information, please visit the official [MToolBox Wiki](https://github.com/mitoNGS/MToolBox/wiki).

