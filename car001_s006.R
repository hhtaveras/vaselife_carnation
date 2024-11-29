library(onemap)
library(fullsibQTL)
library(ggplot2)
library(gridExtra)

###################################################################################################################################################################################################################

pop1 <- onemap_read_vcfR(vcf = "input/car001_d011.vcf.gz",  # filtered markers with ind. and parents pop1
                          parent1 = "i089",                              
                          parent2 = "i088",
                          cross = "outcross")

# correct name of last ind. 
rownames(pop1$geno)[nrow(pop1$geno)]<-gsub(" ", "", rownames(pop1$geno)[nrow(pop1$geno)])

# read in phenotypic data
pheno<-readRDS("input/car001_d008.Rds")
colnames(pheno)<-c("ind", "VL")

# match phenotypes with genotypic data
pheno1<-pheno[pheno$ind %in% rownames(pop1$geno),]

# add phenotypic data to OneMap object 
pheno1<-data.frame(VL=pheno1$VL, row.names=pheno1$ind)
pop1$pheno<-pheno1

# read in genotype and phenotype object together with genetic map
fsib1<-create_fullsib( pop1,
                       map.list = readRDS("input/car001_d003.Rds"),
                       step = 1,
                       error.prob = 0.1,
                       map.function = "kosambi")

# select a maximum of 10 cofactors for CIM 
set.seed(123)
cofs_fs1 <- cof_selection( fsib1, 
                           pheno.col = 1,
                           thres.effect = 0.05,
                           n.cofactor = 10,
                           k=2)

# carry out composite interval mapping
cim1 <- cim_scan( fullsib = cofs_fs1, lg = "all", ws = 15, pheno.col = 1, LOD = TRUE )

# permutation test for significant threshold, takes a long time! Use saved file
cim_perm1 <- cim_scan( fullsib = cofs_fs1, lg = "all", pheno.col = 1, 
                       LOD = TRUE, n.perm = 1000, ws = 15,
                       write.perm = "input/car001_d021.txt" ) # save CIM scan

# save permutation test results
save( cim_perm1, file = "input/car001_d021.Rdata" )

# load permutation test results
load( "input/car001_d021.Rdata" )

# plot CIM results
(maxm<-round(max(cim1[,3]), 1))

p1<-plot_fullsibQTL(fullsib = fsib1, fullsib.scan = cim1, grayscale = T, lgs = 1:5, thr = summary(cim_perm1, alpha = 0.05)[1,1]) + 
  theme(legend.position='none') +
  labs(y = " ") +
  labs(x = " ") +
  scale_y_continuous(limits = c(-1, maxm), labels = c(as.character(seq(0, maxm))), breaks = c(seq(0, maxm))) +
  theme(plot.margin=unit(c(1,1,0,1),"cm"))
p2<-plot_fullsibQTL(fullsib = fsib1, fullsib.scan = cim1, grayscale = T, lgs = 6:10, thr = summary(cim_perm1, alpha = 0.05)[1,1]) +
  theme(legend.position='none') +
  labs(title = " ") + 
  labs(x = " ") +
  scale_y_continuous(limits = c(-1, maxm), labels = c(as.character(seq(0, maxm))), breaks = c(seq(0, maxm))) +
  theme(plot.margin=unit(c(0.5,1,0.5,1),"cm"))
p3<-plot_fullsibQTL(fullsib = fsib1, fullsib.scan = cim1, grayscale = T, lgs = 11:15, thr = summary(cim_perm1, alpha = 0.05)[1,1]) +
  theme(legend.position='none') +
  labs(y = " ") +
  labs(title = " ") +
  scale_y_continuous(limits = c(-1, maxm), labels = c(as.character(seq(0, maxm))), breaks = c(seq(0, maxm))) +
  theme(plot.margin=unit(c(0,1,1,1),"cm"))

p4<-grid.arrange(p1, p2, p3)

ggsave(p4,
       file = "figures/figS3.pdf", width=10, height=8, units="in")


# QTL characterization 
qtl1<-as.data.frame(summary(cim1))
qtls.cim1 <- r2_ls( fsib1, pheno.col = 1, lg = c(qtl1$lg),
                   pos = c(rownames(qtl1)))
qtls.cim1
# for LG 5 only

which(cim1[,3] > summary(cim_perm1, alpha = 0.05)[1,1])
qtl1_5_all<-as.data.frame(cim1[which(cim1[,3] > summary(cim_perm1, alpha = 0.05)[1,1]),])

# five peaks on chromosome 5
qtl1_5_peaks<-c("loc4", "M2999911757", "M3001144616", "M3001744512", "M3001782642")
qtl1_5_all[rownames(qtl1_5_all) %in% qtl1_5_peaks,]
qtls.cim1_5 <- r2_ls( fsib1, pheno.col = 1, lg = rep(5, length(qtl1_5_peaks)),
                    pos = qtl1_5_peaks)

# phenotypic variation explaine by QTL
qtls.cim1_5 


##########################################################################################################################
##########################################################################################################################

pop2 <- onemap_read_vcfR(vcf = "input/car001_d012.vcf.gz",   # filtered markers with ind. and parents pop2  
                          parent1 = "i020",                 
                          parent2 = "i073",
                          cross = "outcross")
# correct name of last ind. 
rownames(pop2$geno)[nrow(pop2$geno)]<-gsub(" ", "", rownames(pop2$geno)[nrow(pop2$geno)])

# read in phenotypic data
pheno<-readRDS("input/car001_d008.Rds")
colnames(pheno)<-c("ind", "VL")

# match phenotypes with genotypic data
pheno2<-pheno[pheno$ind %in% rownames(pop2$geno),]

# add phenotypic data to OneMap object 
pheno2<-data.frame(VL=pheno2$VL, row.names=pheno2$ind)
pop2$pheno<-pheno2

# read in genotype and phenotype object together with genetic map
fsib2<-create_fullsib( pop2,
                       map.list = readRDS("input/car001_d003.Rds"),
                       step = 1,
                       error.prob = 0.1,
                       map.function = "kosambi")

# select a maximum of 10 cofactors for CIM 
set.seed(123)
cofs_fs2 <- cof_selection( fsib2, 
                           pheno.col = 1,
                           thres.effect = 0.05,
                           n.cofactor = 10,
                           k=2)

# carry out composite interval mapping
cim2 <- cim_scan( fullsib = cofs_fs2, lg = "all", ws = 15, pheno.col = 1, LOD = TRUE )

# permutation test for significant threshold, takes a long time! Use saved file
cim_perm2 <- cim_scan( fullsib = cofs_fs2, lg = "all", pheno.col = 1, 
                       LOD = TRUE, n.perm = 1000, ws = 15,
                       write.perm = "input/car001_d022.txt" ) # save CIM scan

# save permutation test results
save( cim_perm2, file = "input/car001_d022.Rdata" )

# load permutation test results
load( "input/car001_d022.Rdata" )

# plot CIM results
maxm<-round(max(cim2[,3]), 0)

p5<-plot_fullsibQTL(fullsib = fsib2, fullsib.scan = cim2, grayscale = T, lgs = 1:5, thr = summary(cim_perm2, alpha = 0.05)[1,1]) + 
  theme(legend.position='none') +
  labs(y = " ") +
  labs(x = " ") +
  scale_y_continuous(limits = c(-1, maxm), labels = c(as.character(seq(0, maxm))), breaks = c(seq(0, maxm))) +
  theme(plot.margin=unit(c(1,1,0,1),"cm"))
p6<-plot_fullsibQTL(fullsib = fsib2, fullsib.scan = cim2, grayscale = T, lgs = 6:10, thr = summary(cim_perm2, alpha = 0.05)[1,1]) +
  theme(legend.position='none') +
  labs(title = " ") + 
  labs(x = " ") +
  scale_y_continuous(limits = c(-1, maxm), labels = c(as.character(seq(0, maxm))), breaks = c(seq(0, maxm))) +
  theme(plot.margin=unit(c(0.5,1,0.5,1),"cm"))
p7<-plot_fullsibQTL(fullsib = fsib2, fullsib.scan = cim2, grayscale = T, lgs = 11:15, thr = summary(cim_perm2, alpha = 0.05)[1,1]) +
  theme(legend.position='none') +
  labs(y = " ") +
  labs(title = " ") +
  scale_y_continuous(limits = c(-1, maxm), labels = c(as.character(seq(0, maxm))), breaks = c(seq(0, maxm))) +
  theme(plot.margin=unit(c(0,1,1,1),"cm"))

p8<-grid.arrange(p5, p6, p7)

ggsave(p8,
       file = "figures/figS4.pdf", width=10, height=8, units="in")


# for LG 11 only
qtl2_11_all<-as.data.frame(cim2[which(cim2[,3] > summary(cim_perm2, alpha = 0.05)[1,1]),])
qtl2_11_all

# five peaks on chromosome 5
qtl2_11_peaks<-c("M2999852250", "M3000121965")
qtls.cim2_11 <- r2_ls( fsib2, pheno.col = 1, lg = rep(11, length(qtl2_11_peaks)),
                      pos = qtl2_11_peaks)
qtls.cim2_11
