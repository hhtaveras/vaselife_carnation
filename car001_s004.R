library(qqman)
library(GWASTools)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(SelectionTools)
library(rrBLUP)
#####################
# read in genetic map 

gmap<-readRDS("input/car001_d004.Rds")

head(gmap)

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
x <- st.marker.data.statistics()
geno <- x$genotypes

geno[1:5, 1:5]

# combine genetic map with genotype data 
mapgeno<-merge(x=gmap, y=geno, by.x = "mk.names", by.y = "Mar/Ind")
mapgeno<-cbind(mk.names=mapgeno$mk.names, lg=mapgeno$lg, dist=mapgeno$dist, mapgeno[,4:dim(mapgeno)[2]])
mapgeno<-mapgeno[order(as.numeric(mapgeno$lg), as.numeric(mapgeno$dist)),]
rownames(mapgeno)<-NULL

mapgeno[1:15, 1:15]

# recode genotypes for rrBLUP package

mapgeno[mapgeno == "1/1"] <- -1
mapgeno[mapgeno == "2/2"] <- 1
mapgeno[mapgeno == "1/2"] <- 0
mapgeno[mapgeno == "2/1"] <- 0
mapgeno[mapgeno == "-1/-1"] <- NA
mapgeno[,2:dim(mapgeno)[2]]<-lapply(mapgeno[,2:dim(mapgeno)[2]], as.numeric)

mapgeno[1:15, 1:15]

# calculate relationship matrix and inpute markers 

X <- data.matrix (mapgeno[,4:ncol(mapgeno)])
rownames(X) <- mapgeno$mk.names

imputed<-rrBLUP::A.mat(t(X), impute.method = "Mean", return.imputed = T)
mapgenoi<-data.frame(mk.names=mapgeno$mk.names, lg=mapgeno$lg, dist=mapgeno$dist, t(imputed$imputed), row.names = NULL)
A<-imputed$A
mapgenoi[1:15, 1:15]

##########################
# read in phenotypic data  

pheno<-readRDS("input/car001_d008.Rds")

# remove individuals missing in mapgeno 
mapgenoi<-mapgenoi[,!colnames(mapgenoi) %in% setdiff(colnames(mapgenoi)[4:dim(mapgenoi)[2]], pheno$GENO)]
pheno<-pheno[!pheno$GENO %in% setdiff(pheno$GENO, colnames(mapgenoi)[4:dim(mapgenoi)[2]]),]
A<-A[!rownames(A) %in% setdiff(rownames(A), pheno$GENO),!colnames(A) %in% setdiff(colnames(A), pheno$GENO)]

colnames(pheno)<-c("ind", "VL")

head(pheno,15)
dim(pheno)
A[1:10, 1:10]
dim(A)

##############
# actual GWAS 

m11 <- rrBLUP::GWAS ( pheno = pheno[,c("ind", "VL")],
                      geno = mapgenoi,
                      P3D = T,
                      plot = F,
                      K = A,
                      n.PC = 1)
# put GWAS results as data frame
p.v <- data.frame( SNP = m11$mk.names,
                   CHR = m11$lg,
                   BP = m11$dist,
                   P = 10^ (-m11[,4] ))
# manhattan plot 
p1<-~{qqman::manhattan ( p.v, main="A",
                         adj = 0,             
                         ylim = c(0, 4),
                         #suggestiveline = 10,
                         #genomewideline = (-log10(0.01)),
                         #col = c("gray10", "gray60")
                         col = c("#000080", "#6495ED"),
                         xlab = "Linkage group")
#abline(a = -log10(0.01), b = 0, col = "#FF6347", lty = 5, lwd = 1 )
}

# QQ Plot
p2<-~{qqPlot( 10^(-m11[,4]),
              ylim = c(0,3.61),
              xlim = c(0,3.61),
              xaxt="n", 
              yaxt="n",
              main="B", truncate = F, adj =0)
axis(1, at=c(0, 1, 2, 3))
axis(2, at=c(0, 1, 2, 3))
}

# put both plots together
(p3<-plot_grid(p1, p2, rel_widths = c(2,1)))

# save plot
# ggsave(p3,
#        file = "figures/figS3.png", 
#        width=12, height=4.5, units="in")

########################################################################
# get top associated markers from GWAS to include in genomic prediction

log<- -log10(0.1) #  -log 10 (p-value) 

idx <- (m11[,4] > log)
sig <- m11[idx ,"mk.names"]
pval <- m11[idx ,"VL"]
lg <- m11[idx ,"lg"]
dist <- m11[idx ,"dist"]

top<-data.frame(name = sig, lg=lg, position = dist,  pvalue = pval)
top<-top[order(-top$pvalue),]
dim(top)

# get marker maf
genoo<-mapgeno[,-c(2:3)]
mafs<-c()
for(i in 1:nrow(genoo)){
  if( length(table(as.matrix(genoo[i,2:ncol(genoo)]))) == 3 ){  
    (p<-table(as.matrix(genoo[i,2:ncol(genoo)]))[1]*2+table(as.matrix(genoo[i,2:ncol(genoo)]))[2])
    (q<-table(as.matrix(genoo[i,2:ncol(genoo)]))[3]*2+table(as.matrix(genoo[i,2:ncol(genoo)]))[2])
    (maf<-q/(p+q))
    mafs<-append(mafs,maf)
  } else if ( length(table(as.matrix(genoo[i,2:ncol(genoo)]))) == 2 ){
    (p<-table(as.matrix(genoo[i,2:ncol(genoo)]))[1]*2+table(as.matrix(genoo[i,2:ncol(genoo)]))[2])
    (q<-table(as.matrix(genoo[i,2:ncol(genoo)]))[2])
    (maf<-q/(p+q))
    mafs<-append(mafs,maf)
  }
}
genomaf<-as.data.frame(cbind(name=genoo$mk.names, maf=mafs))
rownames(genomaf)<-NULL
tgenomaf<-genomaf[genomaf$name %in% top$name,]
tgenomafo<-tgenomaf[match(top$name, tgenomaf$name),]
top$maf<-tgenomafo$maf

# save
#saveRDS(top, "input/car001_d009.Rds")



top2<-top[top$pvalue>2,]
top2[order(top2$lg, top2$position),]

















