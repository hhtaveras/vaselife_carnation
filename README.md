# Scripts
-**car001_s001.R** -> build genetic map with OneMap combining the two F1 carnation populations, **car001_d001.vcf** and **car001_d002.vcf**. Make plot of built genetic map car001_d004.Rds with *LinkageMapView*.<br>

-**car001_s002.R** -> calculate adjusted entry means for phenotypic data of Boxriker *et al*. 2018 and get the ones for individuals of population 1 and 2 only.<br>

-**car001_s003.R** -> carry out GWAS with *rrBLUP* and make manhattan and QQ plots, using created genetic map **car001_d004.Rds**, genotypes of both populations put together **car001_d005.pop** filtered with *SelectionTools*, phenotypes of both populations combined.<br>

-**car001_s004.R** -> carry put genomic prediction including top associated GWAS markers as fixed effects, for both only additive or additive + dominance effects.<br>

-**car001_s005.R** -> make plot of principal coordinate analysis and heatmap genetic distances between individuals, modified roger’s distance. Make correlation plot between genetic and physical distances.<br>

-**car001_s006.R** -> carry out QTL mapping for the two populations individually and make plots.<br>   

-**car001_s010.R** -> modified functions of *LinkageMapView*.<br>
