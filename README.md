# Scripts
-**car001_s001.R** -> build genetic map with OneMap combining the two F1 carnation populations, **input/car001_d001.vcf.gz** and **input/car001_d002.vcf.gz**. Make plot of built genetic map **input/car001_d004.Rds** with *LinkageMapView*.<br>

-**car001_s002.R** -> calculate adjusted entry means for phenotypic data of Boxriker *et al*. 2018 and get the ones for individuals of population 1 and 2 only.<br>

-**car001_s003.R** -> carry out GWAS with *rrBLUP* and make manhattan and QQ plots, using built genetic map **input/car001_d004.Rds**, genotypes of both populations put together **input/car001_d005.pop** filtered with *SelectionTools*, phenotypes of both populations combined.<br>

-**car001_s004.R** -> carry put genomic prediction including top associated GWAS markers as fixed effects, for both only additive and additive + dominance effects.<br>

-**car001_s005.R** -> make plot of principal coordinate analysis and heatmap genetic distances between individuals. Make correlation plot between genetic and physical distances.<br>

-**car001_s006.R** -> carry out QTL mapping for the two populations individually and make plots.<br>   

-**car001_s010.R** -> modified functions of *LinkageMapView*.<br>
# input
-**car001_d001.vcf.gz** -> population 1 in vcf format to read into OneMap.

-**car001_d002.vcf.gz** -> population 2 in vcf format to read into OneMap.

-**car001_d003.Rds** -> built genetic map in OneMap format with populations 1 and 2. 5412 SNP markers.

-**car001_d004.Rds** -> built genetic map in R dataframe format for GWAS calculations.

-**car001_d005.pop** -> SNP marker data to read into *SelectionTools* for filtering, 13917 SNP markers and 163 individuals.

-**car001_d006.Rds** -> phenotypic data of complete trial (Boxriker *et al*. 2018) to estimate entry means with *ASReml-R*. 3368 entries, 12 variables: vase life in days (y), greenhouse temperature in C° (GT), replication (REP), type of carnation (CT), storage treatment or not (STO), genotype (GENO), box in the greenhouse (BOX), position inside box (P), day of starting second phase in laboratory (DV), tray in the laboratory (TRAY), position in tray (POS), stem harvested, first or second (W). 556 Genotypes (GENO) in total. 

-**car001_d007.Rds** -> list individuals in populations 1 and 2 only.

-**car001_d008.Rds** -> adjusted entry means of vase life for individuals of populations 1 and 2 only to use for GWAS.  

-**car001_d009.Rds** -> top scoring markers log10(0.1) from GWAS to include in genomic prediction model. 

-**car001_d010.Rds** -> genetic and physical map positions of mapped markers. 

-**car001_d011.vcf.gz** -> genotypic data of population 1 after being filtered for quality SNPs.

-**car001_d012.vcf.gz** -> genotypic data of population 2 after being filtered for quality SNPs.

-**car001_d021.Rdata** -> permutation test results for CIM scan of population 1. 1000 permutations. 

-**car001_d022.Rdata** -> permutation test results for CIM scan of population 2. 1000 permutations. 

# output
