library(SelectionTools)
library(sommer)
library(WeightIt)
library(Matrix)
library(MASS)
library(crayon)
library(dplyr)
library(ggplot2)


###################################
# read in and filter genotypic data

st.input.dir <- "input"
st.output.dir <- "output"

# Read marker data
st.read.marker.data("car001_d005.pop", format = "m")

# Preprocess data
st.restrict.marker.data (NoAll.MAX = 2) # Maximum number of alleles
st.restrict.marker.data (MaMis.MAX = 0.1) # Max missing at a marker
st.restrict.marker.data (ExHet.MIN = 0.1) # Minimum gene diversity
st.restrict.marker.data (InMis.MAX = 0.1) # Max missing per individual
x<-st.marker.data.statistics()

# make random additive effects matrix
ZZ <- gs.build.Z(data.set="default", out.filename="Z.matrix", auxfiles=T)

# recode data to rrBLUP format 
ZZ[ZZ == 0] <- -1
ZZ[ZZ == 1] <- 0
ZZ[ZZ == 2] <- 1
ZZ[1:10,1:10]

##########################
# read in phenotypic data  

pheno<-readRDS("input/car001_d008.Rds")
colnames(pheno)<-c("ind", "VL")

# remove genotypes in phenotypic data that are not present in genotypic data
pheno<-pheno[!pheno$ind %in% setdiff(pheno$ind, rownames(ZZ)),]

# make yy
yy<-data.frame(y = pheno$VL)
rownames(yy)<-pheno$ind

# remove genotypes in genotypic data that are not present in phenotypic data
ZZ<-ZZ[!rownames(ZZ) %in% setdiff(rownames(ZZ), rownames(yy)),]

# get same order of genotypes on both 
yy$name<-rownames(yy)
yy<-yy[match(rownames(ZZ), rownames(yy)),]
yy<-data.frame(y=yy$y, row.names = rownames(yy))

# read in top associated markers from GWAS 
top<-readRDS("input/car001_d009.Rds") 

# define number of top markers to be included as fixed effects 
n.tops<-c(0, 1, 2, 5, 10, 25, 50, 75, 100, 120)

# number of cross-validation runs 
NRUNS <- 1000

# dataframe to save correlations per run
corrsAD<-data.frame(matrix(NA, nrow = NRUNS, ncol = 0))

# begin loop
###################################################################################

for (j in n.tops) {
  
  # for no markers added as fixed effects
  if (j == 0) {
    crrs <- c()
    for (i in 1:NRUNS) {
      
      set.seed(i)
      # set training and validation sets,
      #33 individuals for validation set, make validation set y values and additive effect matrix  
      vv <- sample(rownames(yy), round(nrow(yy) / 5))
      yy.val <- data.frame(VL = yy[vv,])
      rownames(yy.val) <- vv
      
      yy.trn<-data.frame(VL = yy[setdiff(rownames(yy), rownames(yy.val)),])
      rownames(yy.trn) <- setdiff(rownames(yy), rownames(yy.val))
      
      ZZ.val <- ZZ[vv,]
      ZZ.trn <- ZZ[setdiff(rownames(ZZ), rownames(ZZ.val)),]
      
      dim(yy.val)
      dim(yy.trn)
      dim(ZZ.val)
      dim(ZZ.trn)
      
      # create dominance effects matrices for training an validation sets
      ZZd <- ZZ
      ZZd.trn <- ZZ.trn
      ZZd.val <- ZZ.val
      
      # enconde for dominance effects
      ZZd[ZZd == 0] <- 2
      ZZd[ZZd == -1] <- 0
      ZZd[ZZd == 1] <- 0
      ZZd[ZZd == 2] <- 1
      
      ZZd.val <- ZZd[vv,]
      ZZd.trn <- ZZd[setdiff(rownames(ZZd), rownames(ZZd.val)),]
      dim(ZZd.trn)
      dim(ZZd.val)
      dim(ZZ.trn)
      dim(ZZd.trn)
      yy.trn$X1 <- yy.trn$VL
      
      # rrBLUP
      tryCatch({
        system.time(ans <- mmer(
          X1 ~ 1,
          random = ~ vsr(list(ZZ.trn), buildGu = F) + vsr(list(ZZd.trn), buildGu = F),  # additive and dominace random effects
          rcov = ~ units,
          nIters = 3,
          data = yy.trn,
          verbose = FALSE
        ))
        
        # predict VL for validation set 
        y.hat <- ZZ.val %*% as.matrix(ans$U$`u:ZZ.trn`$X1) + ZZd.val %*% as.matrix(ans$U$`u:ZZd.trn`$X1)
        
        # compare observed and predicted VL for validation set
        (rr <- cor(yy.val, y.hat))       
        crrs <- append(crrs, rr)
        
      }, error = function(e) {
        cat("ERROR :", conditionMessage(e), "\n")
      })
      
    }
    corrsAD <- cbind(corrsAD, crrs)
    
    # for 1 to 120 markers added as fixed effects  
  } else {
    crrs <- c()
    
    # make fixed effects design matrix with top associated markers
    tops <- head(top, j)
    st.copy.marker.data("top", "default")
    st.restrict.marker.data(mar.list = tops$name, data.set = "top")
    XX <- gs.build.Z(data.set = "top",
                     out.filename = "X.matrix",
                     auxfiles = T)
    
    namecols <- c(colnames(XX))
    
    # recode data to rrBLUP format
    XX[XX == 0] <- -1
    XX[XX == 1] <- 0
    XX[XX == 2] <- 1
    
    # find missing genotypes in phenotypic or genotypic data
    dif <- setdiff(rownames(XX), rownames(yy))
    
    # remove genotypes in genotypic date that are not present in phenotypic data
    XX <- as.matrix(XX[!rownames(XX) %in% dif,], colnames = T)
    colnames(XX) <- namecols
    
    for (i in 1:NRUNS) {
      
      # set training and validation sets
      set.seed(i)
      vv <- sample(rownames(yy), round(nrow(yy) / 5))
      yy.val <- data.frame(VL = yy[vv,])
      rownames(yy.val) <- vv
      yy.trn <- data.frame(VL = yy[setdiff(rownames(yy), rownames(yy.val)),])
      rownames(yy.trn) <- setdiff(rownames(yy), rownames(yy.val))
      
      
      XX.val <- as.matrix(XX[vv,])
      colnames(XX.val) <- namecols
      XX.trn <- as.matrix(XX[setdiff(rownames(XX), rownames(XX.val)),])
      colnames(XX.trn) <- namecols
      ZZ.val <- ZZ[vv,]
      ZZ.trn <- ZZ[setdiff(rownames(ZZ), rownames(ZZ.val)),]
      dim(yy.val)
      dim(yy.trn)
      dim(XX.val)
      dim(XX.trn)
      dim(ZZ.val)
      dim(ZZ.trn)
      
      # check if fixed effects matrix has full rank, if not make into a full rank matrix
      dim(XX.trn)
      if (qr(XX.trn)$rank < ncol(XX.trn)) {
        XX.trn <- make_full_rank(XX.trn)
      }
      dim(XX.trn)
      nc <- c(colnames(XX.trn))
      dim(XX.val)
      XX.val <- as.matrix(XX[vv, colnames(XX.trn)])
      colnames(XX.val) <- nc
      dim(XX.val)
      
      # create dominance effects matrices for training an validation sets
      ZZd <- ZZ
      ZZd.trn <- ZZ.trn
      ZZd.val <- ZZ.val
      ZZd[ZZd == 0] <- 2
      ZZd[ZZd == -1] <- 0
      ZZd[ZZd == 1] <- 0
      ZZd[ZZd == 2] <- 1
      ZZd.val <- ZZd[vv,]
      ZZd.trn <- ZZd[setdiff(rownames(ZZd), rownames(ZZd.val)),]
      dim(yy.trn)
      dim(XX.trn)
      dim(ZZ.trn)
      dim(ZZd.trn)
      
      yy.trn$X1 <- yy.trn$VL
      
      # rrBLUP
      tryCatch({
        system.time(ans <- mmer(
          X1 ~ 1 + XX.trn,                                                              # top associated markers as fixed effects 
          random = ~ vsr(list(ZZ.trn), buildGu = F) + vsr(list(ZZd.trn), buildGu = F),  # additive and dominace random effects
          rcov = ~ units,
          nIters = 3,
          data = yy.trn,
          verbose = FALSE
        ))
        
        # predict VL for validation set
        y.hat <- cbind(rep(1, dim(XX.val)[1]), XX.val) %*% as.matrix(ans$Beta$Estimate) +         
          ZZ.val %*% as.matrix(ans$U$`u:ZZ.trn`$X1) + 
          ZZd.val %*% as.matrix(ans$U$`u:ZZd.trn`$X1)
        
        # compare observed and predicted VL for validation set
        (rr <- cor(yy.val, y.hat))
        crrs <- append(crrs, rr)
        
      }, error = function(e) {
        cat("ERROR :", conditionMessage(e), "\n")
      })
    }
    
    # account for NAs
    if (length(crrs) < NRUNS){
      NNA <- NRUNS-length(crrs)
      crrs<-c(crrs, rep(NA, NNA))
    }
    corrsAD <- cbind(corrsAD, crrs)
  }
  
}        

colnames(corrsAD)<- c(paste0("mks_", n.tops))

################################################################################
# only additive effects
################################################################################

corrsA<-data.frame(matrix(NA, nrow = NRUNS, ncol = 0))


# begin loop
###################################################################################

for (j in n.tops) {
  if (j == 0) {
    crrs <- c()
    
    for (i in 1:NRUNS) {
      
      set.seed(i)
      # set training and validation sets,
      #33 individuals for validation set, make validation set y values and additive effect matrix
      vv <- sample(rownames(yy), round(nrow(yy) / 5))
      yy.val <- data.frame(VL = yy[vv,])
      rownames(yy.val) <- vv
      yy.trn <- data.frame(VL = yy[setdiff(rownames(yy), rownames(yy.val)),])
      rownames(yy.trn) <- setdiff(rownames(yy), rownames(yy.val))
      ZZ.val <- ZZ[vv,]
      ZZ.trn <- ZZ[setdiff(rownames(ZZ), rownames(ZZ.val)),]
      dim(yy.val)
      dim(yy.trn)
      dim(ZZ.val)
      dim(ZZ.trn)
      yy.trn$X1 <- yy.trn$VL
      
      # rrBLUP
      tryCatch({
        system.time(ans <- mmer(
          X1 ~ 1,
          random = ~ vsr(list(ZZ.trn), buildGu = F),  # only additive effects
          rcov = ~ units,
          nIters = 3,
          data = yy.trn,
          verbose = FALSE
        ))
        
        # estimate VL for validation set
        y.hat <- ZZ.val %*% as.matrix(ans$U$`u:ZZ.trn`$X1) 
        
        # compare observed and predicted VL for validation set
        (rr <- cor(yy.val, y.hat))      
        crrs <- append(crrs, rr)
        
      }, error = function(e) {
        cat("ERROR :", conditionMessage(e), "\n")
      })
      
    }
    corrsA <- cbind(corrsA, crrs)
    
    # for 1 to 120 markers added as fixed effects
  } else {
    
    crrs <- c()
    # make fixed effects design matrix with top associated markers
    tops <- head(top, j)
    st.copy.marker.data("top", "default")
    st.restrict.marker.data(mar.list = tops$name, data.set = "top")
    XX <- gs.build.Z(data.set = "top",
                     out.filename = "X.matrix",
                     auxfiles = T)
    namecols <- c(colnames(XX))
    
    # recode data to rrBLUP format
    XX[XX == 0] <- -1
    XX[XX == 1] <- 0
    XX[XX == 2] <- 1
    
    # find missing genotypes in phenotypic or genotypic data
    setdiff(rownames(yy), rownames(XX))
    setdiff(rownames(XX), rownames(yy))
    dif <- setdiff(rownames(XX), rownames(yy))
    
    # remove genotypes in genotypic date that are not present in phenotypic data
    XX <- as.matrix(XX[!rownames(XX) %in% dif,], colnames = T)
    colnames(XX) <- namecols
    
    # check for correct dimesions
    
    dim(yy)
    dim(XX)
    dim(ZZ)
    
    for (i in 1:NRUNS) {
      
      set.seed(i)
      # set training and validation sets,
      #33 individuals for validation set, make validation set y values and additive effect matrix
      vv <- sample(rownames(yy), round(nrow(yy) / 5))
      yy.val <- data.frame(VL = yy[vv,])
      rownames(yy.val) <- vv
      yy.trn <-data.frame(VL = yy[setdiff(rownames(yy), rownames(yy.val)),])
      rownames(yy.trn) <- setdiff(rownames(yy), rownames(yy.val))
      
      XX.val <- as.matrix(XX[vv,])
      colnames(XX.val) <- namecols
      XX.trn <-
        as.matrix(XX[setdiff(rownames(XX), rownames(XX.val)),])
      colnames(XX.trn) <- namecols
      ZZ.val <- ZZ[vv,]
      ZZ.trn <- ZZ[setdiff(rownames(ZZ), rownames(ZZ.val)),]
      dim(yy.val)
      dim(yy.trn)
      dim(XX.val)
      dim(XX.trn)
      dim(ZZ.val)
      dim(ZZ.trn)
      
      # check if fixed effects matrix has full rank, if not make into a full rank matrix
      if (qr(XX.trn)$rank < ncol(XX.trn)) {
        XX.trn <- make_full_rank(XX.trn)
      }
      dim(XX.trn)
      nc <- c(colnames(XX.trn))
      dim(XX.val)
      XX.val <- as.matrix(XX[vv, colnames(XX.trn)])
      colnames(XX.val) <- nc
      dim(XX.val)
      yy.trn$X1 <- yy.trn$VL
      
      
      # rrBLUP
      tryCatch({
        system.time(ans <- mmer(
          X1 ~ 1 + XX.trn,                            # top associated markers as fixed effects 
          random = ~ vsr(list(ZZ.trn), buildGu = F),  # additive random effects only
          rcov = ~ units,
          nIters = 3,
          data = yy.trn,
          verbose = FALSE
        ))
        
        # estumate VL for validation set
        y.hat <- cbind(rep(1, dim(XX.val)[1]), XX.val) %*% as.matrix(ans$Beta$Estimate) +         
          ZZ.val %*% as.matrix(ans$U$`u:ZZ.trn`$X1) 
        
        # compare observed and predicted VL for validation set
        (rr <- cor(yy.val, y.hat))
        crrs <- append(crrs, rr)
        
      }, error = function(e) {
        cat("ERROR :", conditionMessage(e), "\n")
      })
      
    }
    # account for NAs
    if (length(crrs) < NRUNS){
      NNA <- NRUNS-length(crrs)
      crrs<-c(crrs, rep(NA, NNA))
    }
    corrsA <- cbind(corrsA, crrs)
  }
  
}        

colnames(corrsA)<- c(paste0("mks_", n.tops))

ggcorrs<-data.frame(matrix(NA, nrow = 3, ncol = 0))

for(i in 1:10){
  
  ys<-corrsA[,i]
  xs<-rep(i, length(corrsA[,i]))
  gp<-rep(1, length(corrsA[,i]))
  
  df<-data.frame(xs, ys, gp)
  
  ggcorrs<-rbind(ggcorrs, df)
  
}

dim(ggcorrs)

for(i in 1:10){
  
  ys<-corrsAD[,i]
  xs<-rep(i, length(corrsAD[,i]))
  gp<-rep(2, length(corrsAD[,i]))
  
  df<-data.frame(xs, ys, gp)
  
  ggcorrs<-rbind(ggcorrs, df)
  
}

dim(ggcorrs)

# plot correlations

# get means and medians
corrs<-data.frame(mks0A=corrsA$mks_0, mks0AD=corrsAD$mks_0,
                  mks1A=corrsA$mks_1, mks1AD=corrsAD$mks_1,
                  mks2A=corrsA$mks_2, mks2AD=corrsAD$mks_2,
                  mks5A=corrsA$mks_5, mks5AD=corrsAD$mks_5,
                  mks10A=corrsA$mks_10, mks10AD=corrsAD$mks_10,
                  mks25A=corrsA$mks_25, mks25AD=corrsAD$mks_25,
                  mks50A=corrsA$mks_50, mks50AD=corrsAD$mks_50,
                  mks75A=corrsA$mks_75, mks75AD=corrsAD$mks_75,
                  mks100A=corrsA$mks_100, mks100AD=corrsAD$mks_100,
                  mks120A=corrsA$mks_120, mks120AD=corrsAD$mks_120)

means<-c()
meds<-c()
for(i in 1:20){
  
  means<-append(means, mean(corrs[,i], na.rm = T))
  meds<-append(meds, median(corrs[,i], na.rm = T))
  
}

# position and numbers after decimal of means and medians
ggcorrs2 = distinct(ggcorrs, gp, xs) %>%
  arrange(gp, xs)

ggcorrs2$yloc = max(ggcorrs$ys, na.rm = T) + .065
ggcorrs2$yloc2 = max(ggcorrs$ys, na.rm = T) + .03
ggcorrs2$xloc = sort(c(seq(0.8, 9.8, 1), seq(1.2, 10.2, 1)), decreasing = F)

(ggcorrs2$means = c(substr(means, 1, 4)))
(ggcorrs2$meds = c(substr(meds, 1, 4)))


# plot
p <- ggplot(ggcorrs, aes(x=factor(xs), y=ys))
p2<-p + geom_boxplot(aes(fill = factor(gp)), width = 0.69, show.legend = T) +
  labs(fill = " ") +
  #scale_fill_manual(values=c("#FB9A99", "#A6CEE3"), labels = c("A", "A + D")) +
  scale_fill_manual(values=c("gray85", "gray69"), labels = c("A", "A + D")) +
  theme_bw() +
  scale_y_continuous(name = "r(y,ŷ)",
                     #expand = expansion(mult = c(0.5, 0.1)),
                     limits = c(min(ggcorrs$ys, na.rm = T), 1),
                     breaks = seq(0, 1, 0.5), 
                     labels = c("0.0", "0.5", "1.0")
  ) +
  scale_x_discrete(name ="No. of markers used as fixed effects",
                   expand = expansion(mult = c(0.12, 0)),
                   labels=c("None","1","2", "5", "10", "25", 
                            "50", "75", "100", "120")) +
  geom_hline(yintercept = median(corrsA[,1]), colour = "#FF6347",
             size= 0.25, linetype = "longdash") + 
  geom_text(data = ggcorrs2, aes(y = yloc, x = xloc, label = means), 
            position = position_dodge(width = .75)) + 
  geom_text(data = ggcorrs2, aes(y = yloc2, x = xloc, label = meds)) + 
  
  annotate("text", x = 0.3, y=ggcorrs2$yloc[1]+0.0022, label = "μ ",
           fontface = 2) +
  annotate("text", x = 0.3, y=ggcorrs2$yloc2[1]+0.0013, label = "Z ",
           fontface = 2) +  
  theme(axis.line = element_line(color='black'), 0,
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color = NA),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        axis.ticks.x=element_blank(),
        #axis.title.y = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(vjust=2), 
        axis.title.x = element_text(margin = margin(t = 12)),
        legend.position = c(0.043,0.14),
        panel.grid.major = element_blank())

p2

# ggsave( p2 ,
#         file="figures/fig5.png", 
#         width=12, height=6, units="in")








