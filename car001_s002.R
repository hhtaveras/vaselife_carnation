########################
# calculate entry means
########################
library(asreml)
set.seed(1)


data<-readRDS("input/car001_d006.Rds")

# mixed model 
m1<-asreml(fixed = Y ~ GT + REP + STO + CT + CT:GENO,
           random = ~ STO:CT:GENO + REP:BOX + REP:BOX:P + DV + DV:TRAY + DV:TRAY:POS,
           na.action = na.method("omit"),
           start.values = F,
           workspace = 80e6,
           data = data)

# get entry means
p1<-predict(m1, classify = "CT:GENO", sed = T, 
            pworkspace = 80e6)

# get entry means of only the invividuals on populations 1 and 2

ct1<-subset(p1$pvals, p1$pvals$CT == 1) # get genotypes corresponding to CT 1 

names<-readRDS("input/car001_d007.Rds")

pheno<-ct1[!ct1$GENO %in% setdiff(ct1$GENO, names),2:3]
pheno$GENO<-as.character(pheno$GENO)
rownames(pheno)<-NULL
pheno

saveRDS(pheno, "input/car001_d008.Rds")
