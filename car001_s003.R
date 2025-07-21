library(onemap)
library(fullsibQTL)
library(ggplot2)
library(gridExtra)
library(knitr)

###################################################################################################################################################################################################################

pop1 <- onemap_read_vcfR(vcf = "input/car001_d011.vcf.gz",  # filtered markers with ind. and parents pop1
                         parent1 = "i089",                              
                         parent2 = "i088",
                         cross = "outcross")

table(pop1$geno, useNA = "always")
pop1$geno[pop1$geno == 0] <- 3

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

# select and define a maximum of 10 cofactors for CIM 
set.seed(123)
cofs_fs1 <- cof_selection( fsib1, 
                           pheno.col = 1,
                           thres.effect = 0.05,
                           n.cofactor = 10,
                           k=2)

cofs_fsdef1 <- cof_definition( fsib1, pheno.col = 1, 
                               thres.effect = 0.05, 
                               cof.pos = as.matrix(cofs_fs1$cofactors$names.cof))

# carry out composite interval mapping
cim1 <- cim_scan( fullsib = cofs_fsdef1, lg = "all", ws = 15, pheno.col = 1, LOD = TRUE )

# permutation test for significant threshold, takes a long time! Use saved file
# cim_perm1 <- cim_scan( fullsib = cofs_fs1, lg = "all", pheno.col = 1, 
#                        LOD = TRUE, n.perm = 1000, ws = 15,
#                        write.perm = "input/car001_d021.txt" ) # save CIM scan

# save permutation test results
#save( cim_perm1, file = "input/car001_d021.Rdata" )

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

 # ggsave(p4,
 #        file = "figures/fig2.pdf", width=10, height=8, units="in")


# QTL characterization 
qtl1<-as.data.frame(summary(cim1))
qtls.cim1 <- r2_ls( fsib1, pheno.col = 1, lg = c(qtl1$lg),
                    pos = c(rownames(qtl1)))
qtls.cim1
# for LG 5 only

which(cim1[,3] > summary(cim_perm1, alpha = 0.05)[1,1])
(qtl1_5_all<-as.data.frame(cim1[which(cim1[,3] > summary(cim_perm1, alpha = 0.05)[1,1]),]))

# five peaks on chromosome 5
(qtl1_5_all<-cim1[cim1[,1] == 5 & cim1[,3] > 5,])
qtl1_5_peaks<-c("M2999911757", "M3001144616", "M3001744512")
qtl1_5_all[rownames(qtl1_5_all) %in% qtl1_5_peaks,]
(qtls.cim1_5 <- r2_ls( fsib1, pheno.col = 1, lg = rep(5, length(qtl1_5_peaks)),
                      pos = qtl1_5_peaks))

# for LG 10 only # one peak on chromosome 10
(qtl1_10_all<-cim1[cim1[,1] == 10 & cim1[,3] > 5,])

qtl1_10_peaks<-c("M3001191510")
qtl1_10_all[rownames(qtl1_10_all) %in% qtl1_10_peaks,]
qtls.cim1_10 <- r2_ls( fsib1, pheno.col = 1, lg = rep(10, length(qtl1_10_peaks)),
                      pos = qtl1_10_peaks)

# phenotypic variation peaks lgs 5 and 10 

(qtl1_510<-c(qtl1_5_peaks, qtl1_10_peaks))
qtls.cim1_510 <- r2_ls( fsib1, pheno.col = 1, lg = c(5, 5, 5, 10),
                       pos = qtl1_510)


# get effects of peak QTLs

effs1<-cim_char( cofs_fsdef1, pheno.col = 1, lg = qtls.cim1_510$lg[2], pos = qtls.cim1_510$loc[2])

effs1<-cim_char( cofs_fs1, pheno.col = 1, lg = 5, pos = qtl1_5_peaks[1])
for(i in 3:length(qtls.cim1_510$lg)){
  effs1<- cbind(effs1, cim_char( cofs_fs1, pheno.col = 1, lg = qtls.cim1_510$lg[i], pos = qtls.cim1_510$loc[i]))
}
kable(effs1, digits = 2)

# get MAF of peak QTLs
mafs<-c()
for (i in qtls.cim1_510[2:length(qtls.cim1_510$lg),2]){
print(i)
if(length(table(pop1$geno[,i])) == 2){

p<-(table(pop1$geno[,i])[1]*2)+(table(pop1$geno[,i])[2]*1)
q<-(table(pop1$geno[,i])[1]*0)+(table(pop1$geno[,i])[2]*1)
maf<-min(p, q)/(p+q)
mafs<-append(mafs, maf)
} else if (length(table(pop1$geno[,i])) == 3){
  
  p<-(table(pop1$geno[,i])[1]*2) + (table(pop1$geno[,i])[2]*1) + (table(pop1$geno[,i])[3]*0) 
  q<-(table(pop1$geno[,i])[1]*0) + (table(pop1$geno[,i])[2]*1) + (table(pop1$geno[,i])[3]*2)
  maf<-min(p, q)/(p+q)
  mafs<-append(mafs, maf)
}
}  
  
mafs

# summary table of QTL for pop1

tab1<-data.frame(SNP_ID = qtls.cim1_510[2:length(qtls.cim1_510$lg),2], lg=qtls.cim1_510[2:length(qtls.cim1_510$lg),1],
           pos=effs1[2,1:dim(effs1)[2]], maf = mafs, 
           pval = cim1[qtls.cim1_510[2:length(qtls.cim1_510$lg),2],3],
           eff = effs1[8,1:dim(effs1)[2]],
           pve = qtls.cim1_510[2:length(qtls.cim1_510$lg),3]
           )
rownames(tab1)<-NULL
kable(tab1)

# sort and save top scoring SNP markers for QTL-assisted GP 
cim1
cim1o<-cim1[order(cim1[,3], decreasing = T),]
cim1oo<-cim1o[nchar(rownames(cim1o))>=9,]
cim1oo<-cbind(SNP=rownames(cim1oo), cim1oo)
dim(cim1oo)
head(cim1oo)

#saveRDS(cim1oo, "input/res_cim1_gmap.Rds")

##########################################################################################################################
##########################################################################################################################

pop2 <- onemap_read_vcfR(vcf = "input/car001_d012.vcf.gz",   # filtered markers with ind. and parents pop2  
                         parent1 = "i020",                 
                         parent2 = "i073",
                         cross = "outcross")
table(pop2$geno, useNA = "always")
pop2$geno[pop2$geno == 0] <- 3

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

cofs_fsdef2 <- cof_definition( fsib2, pheno.col = 1, 
                               thres.effect = 0.05, 
                               cof.pos = as.matrix(cofs_fs2$cofactors$names.cof))

# carry out composite interval mapping
cim2 <- cim_scan( fullsib = cofs_fsdef2, lg = "all", ws = 15, pheno.col = 1, LOD = TRUE )

# permutation test for significant threshold, takes a long time! Use saved file
# cim_perm2 <- cim_scan( fullsib = cofs_fs2, lg = "all", pheno.col = 1, 
#                        LOD = TRUE, n.perm = 1000, ws = 15,
#                        write.perm = "input/car001_d022.txt" ) # save CIM scan

# save permutation test results
# save( cim_perm2, file = "input/car001_d022.Rdata" )

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

# ggsave(p8,
#        file = "figures/fig3.pdf", width=10, height=8, units="in")


# for LG 11 only
qtl2_11_all<-as.data.frame(cim2[which(cim2[,3] > summary(cim_perm2, alpha = 0.05)[1,1]),])
qtl2_11_all

# five peaks on chromosome 5
qtl2_11_peaks<-c("M2999852250", "M3000121965")
qtls.cim2_11 <- r2_ls( fsib2, pheno.col = 1, lg = rep(11, length(qtl2_11_peaks)),
                       pos = qtl2_11_peaks)
qtls.cim2_11


# get effects of peak QTLs

effs2<-cim_char( cofs_fs2, pheno.col = 1, lg = qtls.cim2_11$lg[2], pos = qtls.cim2_11$loc[2])

for(i in 3:length(qtls.cim2_11$lg)){
  effs2<- cbind(effs2, cim_char( cofs_fs2, pheno.col = 1, lg = qtls.cim2_11$lg[i], pos = qtls.cim2_11$loc[i]))
}
kable(effs2, digits = 2)

# get MAF of peak QTLs
mafs<-c()
for (i in qtls.cim2_11[2,2]){
  print(i)
  if(length(table(pop2$geno[,i])) == 2){
    
    p<-(table(pop2$geno[,i])[1]*2)+(table(pop2$geno[,i])[2]*1)
    q<-(table(pop2$geno[,i])[1]*0)+(table(pop2$geno[,i])[2]*1)
    maf<-min(p, q)/(p+q)
    mafs<-append(mafs, maf)
  } else if (length(table(pop2$geno[,i])) == 3){
    
    p<-(table(pop2$geno[,i])[1]*2) + (table(pop2$geno[,i])[2]*1) + (table(pop2$geno[,i])[3]*0) 
    q<-(table(pop2$geno[,i])[1]*0) + (table(pop2$geno[,i])[2]*1) + (table(pop2$geno[,i])[3]*2)
    maf<-min(p, q)/(p+q)
    mafs<-append(mafs, maf)
  }
}  

mafs

# summary table of QTL for pop2

tab2<-data.frame(SNP_ID = qtls.cim2_11[2,2], lg=qtls.cim2_11[2,1],
                 pos=effs2[2,1], maf = mafs, 
                 pval = cim2[qtls.cim2_11[2,2],3],
                 eff = effs2[8,1],
                 pve = qtls.cim2_11[2,3]
)
rownames(tab2)<-NULL
kable(tab2)

# both populations together 

tab<-rbind(tab1, tab2)
tab<-cbind(pop=c(rep("pop1", dim(tab1)[1]), "pop2"), tab)
kable(tab, digits = 2, format = "latex")

#write.table(tab, "input/significant_both_pops.txt")

# sort and save top scoring SNP markers for QTL-assisted GP 
cim2
cim2o<-cim2[order(cim2[,3], decreasing = T),]
cim2oo<-cim2o[nchar(rownames(cim2o))>=9,]
cim2oo<-cbind(SNP=rownames(cim2oo), cim2oo)
dim(cim2oo)
head(cim2oo)

#saveRDS(cim2oo, "input/res_cim2_gmap.Rds")

# both populations together
cim12oo<-rbind(cbind(cim1oo, pop=rep("pop1", nrow(cim1oo))), cbind(cim2oo, pop=rep("pop2", nrow(cim2oo))))
cim12oo<-cim12oo[order(as.numeric(cim12oo[,4]), decreasing = T),]

cim12ooo<-cim12oo[!duplicated(cim12oo[,1]),]

dim(cim12ooo)

#saveRDS(cim12ooo, "input/car001_d023.Rds")








