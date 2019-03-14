# tcga_mc3_pipeline
Pipeline developed to analyze intratumor heterogeneity.


### Part 1: Running [TitanCNA](https://github.com/gavinha/TitanCNA.git)

- [MC3 project](https://gdc.cancer.gov/about-data/publications/mc3-2017)

- Download available bam files, get target region bed file, liftover to hg38

- Run TitanCNA
  - Use https://github.com/ATLi2001/titancna.git

- Use Google Cloud to increase computation power
- Need Jinpeng to fill in rest of information


### Part 2: Generating [PhyloWGS](https://github.com/morrislab/phylowgs.git) Inputs

##### Filter VCF file mutations

- Download vcf files and public maf file onto virtual instance
  - Use GDC download tool https://gdc.cancer.gov/access-data/gdc-data-transfer-tool
- Install R onto virtual instance; install packages maftools, vcfR, dplyr
- In reading the maf, the virtual instance will run out of memory. To remedy, trim the maf file; use 37 to keep the HGVSp_Short column
  ```
  cut -f 1-37 mc3.v0.2.8.PUBLIC.maf > mc3.v0.2.8.PUBLIC.filtered.maf
  ```

- Create tsv file with list of all vcf files and the corresponding patient id
- Run vcf_pan12_filter.R
  - Segfault errors may occur; manually delete those patients and skip them

##### Generate SSM and CNV data.txt files

- Use the PhyloWGS program's parser to take in the TitanCNA output and the VCF files
  - Use https://bitbucket.org/merckey/phylowgs.git


### Part 3: Run PhyloWGS

- Use the multievolve.py with default parameters
