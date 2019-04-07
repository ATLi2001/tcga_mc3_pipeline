# tcga_mc3_pipeline
Pipeline developed to analyze intratumor heterogeneity.


### Part 1: Running [TitanCNA](https://github.com/gavinha/TitanCNA.git)

- [MC3 project](https://gdc.cancer.gov/about-data/publications/mc3-2017)

- Download available bam files, get target region bed file, liftover to hg38

- Run TitanCNA
  - Use https://github.com/ATLi2001/titancna.git

- Use Google Cloud 
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
  - Filter syn mutations out
  - Segfault errors may occur; manually delete those patients and skip them

##### Generate SSM and CNV data.txt files

- Use the PhyloWGS program's parser to take in the TitanCNA output and the VCF files
  - Use https://bitbucket.org/merckey/phylowgs.git


### Part 3: Run PhyloWGS

- Set up Google Cloud environment: [Instructions](https://bitbucket.org/merckey/google_cloud/src/6500515126779350a301f327fcc0e5f92a455d57/austin_project.md?fileviewer=file-view-default)
  - for building docker images, use a .dockerignore file: [Instructions](https://docs.docker.com/engine/reference/builder/#dockerignore-file)

- Use the multievolve.py with default parameters
  - n = 4, I = inf, B = 1000, s = 2500, i = 5000, random seed

- Submit using dsub; use --provider google-v2 as this is the most recent version
  ```
  dsub \
    --provider google-v2 \
    --disk-size 0 \
    --project nih-commons-credit-project \
    --zones "us-east1-*" \
    --logging gs://austin_pipeline_mc3/logging/ \
    --image us.gcr.io/nih-commons-credit-project/phylowgs_multi:latest \
    --command '/tmp/phylowgs/run_phylowgs_multi.sh' \
    --tasks google_cloud/phylowgs/pan12_multi.tsv \
    --wait
    ```
### Part 4: Generate Summary Statistics

- Set up Google Cloud environment

- Note that the following scripts are needed and were written by the author; printo.py was edited from the original
  - tree_reader_multi_latest.py - overall call
  - printo2_multi_latest.py - calculate the summary statistics
  - r_medicc_ce.py - calculate the [MEDICC](https://bitbucket.org/rfs/medicc) Clonal Expansion Index
  - printo.py - create text files representing the tree structure
  
- Submit using dsub; note that the file name prefix for each patient includes both the cancer type and patient id
  ```
  dsub \
  --provider google-v2 \
  --disk-size 0 \
  --project nih-commons-credit-project \
  --zones "us-east1-*" \
  --logging gs://austin_pipeline_mc3/logging/ \
  --image us.gcr.io/nih-commons-credit-project/summ_stats_multi:mc3 \
  --command '/tmp/phylowgs/run_tree_reader_multi_latest.sh' \
  --tasks google_cloud/phylowgs/pan12_tree_reader_multi.tsv \
  --wait
  ```
- List of summary statistics
  - Clonal expansion index (ce)
  - Depth of tree
  - Maximum mutations in a single node
  - Number of branches
  - Number of nodes
  - Total mutation number
  - Proportion of mutations in the trunk (trunk proportion)
  - Proportion of cnv's in the trunk (trunk proportion cnv)

### Part 5: Analyze Summary Statistics

- Run scripts submit_summ_stats_multi_analysis.R and submit_ordering_matrix.R
  - the summ_stats_multi_analysis will do clinical analyses including survival, recurrence, stage, gender, age, and smoking
    - some cancer types may be missing available data (i.e. ovarian cancer doesn't have two genders to analyze)
  - the ordering_matrix will create a heatmap with the temporal ordering of the genes
    - genes used in each cancer were determined by PANCAN paper
