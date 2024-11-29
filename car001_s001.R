################################################################
# construct genetic linkage map wiht OneMap 
################################################################

library(onemap)
library(dplyr)
set.seed(123)


pop1 <- onemap_read_vcfR(vcf = "input/car001_d001.vcf",  # only pop1 parents, 
                          parent1 = "i089",                                # all markers 13917
                          parent2 = "i088",
                          cross = "outcross")


pop2 <- onemap_read_vcfR(vcf = "input/car001_d002.vcf.gz",  # only pop2 parents, 
                          parent1 = "i020",                                # all markers 13917  
                          parent2 = "i073",
                          cross = "outcross")

# filter out markers with missing data 
pop1_filtered <- filter_missing(pop1, threshold = 0.25)
pop2_filtered <- filter_missing(pop2, threshold = 0.25)

# combine the two populations into one
comb_pop <- combine_onemap(pop1_filtered, pop2_filtered)

#find and create bin markers for combinen population
find_comb_pop_bins <- find_bins(comb_pop, exact = F)
comb_pop_bins <- create_data_bins(comb_pop, find_comb_pop_bins)

# get suggested LOD score for declaring statistical significance
LOD_sug <- suggest_lod(comb_pop_bins)

# set global error rate
comb_pop_bins_01<-create_probs(comb_pop_bins, global_error = 0.1)

# estimate two point recombination fraction
twopts <- rf_2pts(comb_pop_bins_01, rm_mks = T, LOD = LOD_sug)

# change linkage group names from chrn to just n or NA
for (i in 1:length(twopts$CHROM)){
  if (nchar(twopts$CHROM[i]) == 4){
    cn<-(substr(twopts$CHROM[i], 4,4))
    twopts$CHROM[i]<-cn
  } else if (nchar(twopts$CHROM[i]) == 5){
    cn<-substr(twopts$CHROM[i], 4,5)
    if (cn == "NA"){
      twopts$CHROM[i]<-NA
    } else {
      twopts$CHROM[i]<-cn
    }
  }
}

# assigning unlinked markers to linkage groups if possible
CHR_mks <- group_seq(input.2pts = twopts, 
                     seqs = "CHROM", 
                     unlink.mks = "all", 
                     repeated = FALSE, 
                     LOD = LOD_sug)


# get number total number of markers that will be used for the map

sum(sapply(CHR_mks$sequences, function(z) length(z$seq.num)))

# order markers and calculate genetic distances, takes a long time
set_map_fun(type = "kosambi")
for (i in 1:15){
  CH<-paste0("CH", i)
  arg<- as.character(i)
  CHR<-CHR_mks$sequences[[i]]
  CHR_ord <- order_seq(CHR, touchdown = T, n.init = 5)
  CHR_ord$ord.all[3]<-CHR_ord$ord.all[3]/10
  CHR_frame <- make_seq(CHR_ord, "force")
  assign(paste0("CHR", i, "_final"), CHR_frame)
}

# put linkage groups together
map<-list(CHR1_final, CHR2_final, CHR3_final,
          CHR4_final, CHR5_final, CHR6_final,
          CHR7_final, CHR8_final, CHR9_final,
          CHR10_final, CHR11_final, CHR12_final,
          CHR13_final, CHR14_final, CHR15_final)

saveRDS(map, "input/car001_d003.Rds")

###############################################################
# extract and save only map information 
###############################################################

map<-readRDS("input/car001_d003.Rds")

# export all LGs
for (i in 1:length(map)){
  assign(paste0("lg",i), parents_haplotypes(map[[i]])[,3:4])
}
# put extracted LGs together in a list
lmap<-list(lg1, lg2, lg3, lg4, lg5,
           lg6, lg7, lg8, lg9, lg10,
           lg11, lg12, lg13, lg14, lg15)

# make map dataframe including LG
lgmap<-lmap
for (i in 1:15){
  lgmap[[i]]$lg<-c(rep(i, length(lmap[[i]]$dist)))
}

# save LG 1 
gmap<-lgmap[[1]]

# put all LGs together in one data frame
for (i in 2:15){
  gmap<-union_all(gmap, lgmap[[i]])
}

# summary final map
options(pillar.sigfig = 5)
gmap %>% group_by(lg) %>% summarise(n=n(), 
                                    length=max(dist), 
                                    density=n()/max(dist),
                                    average_distance=max(dist)/n())

saveRDS(gmap, "input/car001_d004.Rds")


################################################################
# draw map with LinkageMapView
################################################################

library(LinkageMapView)
source("car001_s010.R") # read in modified functions for LinkageMapView

gmap<-readRDS("input/car001_d004.Rds")

# rearrage map dataframe
cmap<-data.frame(lg=sprintf("%02d", gmap$lg), pos=gmap$dist, locus=gmap$mk.names)

# name the output path and file 
outfile <- "figures/fig1.pdf"

# define color palette 
sectcoldf <- lmvdencolor(cmap,colorin = c("#9E0142", "#D53E4F", "#F46D43", "#FEE08B", "#ABDDA4", "#3288BD"))

lmv.d.plot(cmap,outfile,denmap=TRUE,sectcoldf=sectcoldf, 
           pdf.width = 7,
           pdf.height = 9,
           cex.axis = 1,
           pdf.pointsize = 13,
           cex.lgtitle = 1, 
           font.axis = 1,
           ylab = "cM")

