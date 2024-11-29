library(SelectionTools)
library(gridExtra)
library(ggplot2)
library(reshape)
library(cowplot)
library(dplyr)


###################################
# read in and filter genotypic data

st.input.dir <- "input"
st.output.dir <- "output"

# Read marker data
st.read.marker.data("car001_d005.pop", format = "m")

#################
# Preprocess data 
st.restrict.marker.data ( NoAll.MAX=2 ) # Maximum number of alleles
st.restrict.marker.data ( MaMis.MAX=0.1 ) # Max missing at a marker
st.restrict.marker.data ( ExHet.MIN=0.1 ) # Minimum gene diversity
st.restrict.marker.data ( InMis.MAX=0.1 ) # Max missing per individual

# Write data in another format
st.write.marker.data(format="n", nfilename = "car001_d015")

# Load marker data
rgs <- read.table("output/car001_d015.npo", header=T)

# order based on names, this puts all ind. from population 4 before population 5
rgso<-rgs[,order(colnames(rgs))]

# Calculate Modified Roger's distance
d1 <- gd.genetic.distance(rgso, measure="mrd")
head(d1)

# Bring data frame into matrix format
p <- gd.data.parameters (rgso)
aa <- d1[,3]
bb <- p$no.pop
cc <- p$pop.list
d2 <- structure(aa, Size = bb, Labels = cc, Diag = FALSE,
                Upper = FALSE, method = "euclidean", class = "dist")
d2

##################################
# Principal coordinates analysis 

d3 <- as.matrix(d2)
n <- ncol(d3)
I <- diag(n)
K <- I-1/n
A <- d3 * d3 * -0.5
B <- K %*% A %*% K
E <- eigen(B)
norm <- t(matrix(sqrt(abs(E$values)),ncol=n,nrow=n))
Y <- norm * E$vectors
Y3 <- Y[,1:3]
Y3 <- data.frame(Y3)
Y3$RGS <- rownames(d3)
colnames(Y3) <- c("PCo1","PCo2","PCo3", "car")
# Explained variance
var.exp <- round (100* (E$values / sum (E$values) ),2)

d4 <- melt(d3)
colnames(d4) <- c("X1", "X2", "MRD")

##########################################
# PCoA and heatmap combined with ggplot2 

# x coordinates for heatmap with ggplot2

p1<-0
p2<-0
for (i in 1:ncol(d3)){
  
  if(as.numeric(substr(colnames(d3)[i], 2,4)) < 224){
    p1<-p1+1
  }else{
    p2<-p2+1
  }  
}
x.coord <- c(0.1, p1+0.3, p1+p2+0.5)

# set theme with white background

theme_set(theme_bw())

# Heatmap
cols <- c("#56B4E9", "#009E73")

xx <- -2 # minimum x and y values for heatmap
sz <- 3  # size of the color-coded axes in heatmap

plot1 <- ggplot(d4, aes(x = X1, 
                        y = X2, fill = MRD)) +
  geom_tile() +
  scale_fill_gradient(low="darkred", high="yellow", name="MRD") + # color scale for heatmap
  
  xlab(NULL) + # no axis labels
  ylab(NULL) +
  
  coord_fixed(ratio = 1, xlim = c(xx, max(x.coord)), ylim = c(xx, max(x.coord)), 
              expand = TRUE, clip = "on") +
  
  # color-coded x axis
  geom_segment(aes(y = xx, yend = xx, x = x.coord[1], xend = x.coord[2]), color = cols[1], size=sz) + 
  geom_segment(aes(y = xx, yend = xx, x = x.coord[2], xend = x.coord[3]), color = cols[2], size=sz) +
  
  
  # color-coded y axis
  geom_segment(aes(y = x.coord[1], yend = x.coord[2], x = xx, xend = xx), color = cols[1], size=sz) + 
  geom_segment(aes(y = x.coord[2], yend = x.coord[3], x = xx, xend = xx), color = cols[2], size=sz) +
  
  ggtitle("b") +
  
  theme(panel.border = element_blank(), # remove black border around the figure
        panel.grid.major = element_blank(), # remove grid
        panel.grid.minor = element_blank(), # remove grid
        
        text=element_text(size=8,  family="sans"), # Helvetica
        legend.text = element_text(size = 8),
        plot.title = element_text(size = 12),
        axis.title = element_text(size = 8),
        axis.ticks = element_blank(), # remove axis ticks
        axis.text = element_blank())  # remove tick labels

# PCoA plot

# assign individuasl to respective population 
pop<-c()
for (i in 1:ncol(d3)){
  
  if(as.numeric(substr(colnames(d3)[i], 2,4)) < 224){
    pop<-append(pop, "p1")
  }else{
    pop<-append(pop, "p2")
  }  
}
Y3$pop<-(pop=as.factor(pop))

# plot
cols <- c("#56B4E9", "#009E73")

plot2 <- ggplot(Y3, aes(x = PCo1, y = PCo2)) + 
  geom_point(color=cols[Y3$pop], shape=Y3$pop, size=2) +
  xlab(label = paste0("PCo 1 (", var.exp[1], " %)")) +
  ylab(label = paste0("PCo 2 (", var.exp[2], " %)")) +
  coord_fixed() + # square plot
  
  ggtitle("a") +
  
  theme(panel.border = element_blank(), # remove black border around the figure
        aspect.ratio=0.9,
        text=element_text(size=8,  family="sans"),
        legend.text = element_text(size = 8),
        plot.title = element_text(size = 12),
        axis.title = element_text(size = 8),
        axis.text = element_text(size=8))

# put both plots together 
plot3 <- plot_grid(plot2, plot1, rel_widths = c(1,1.1))
plot3

# save plot
ggsave( plot3 ,
        file="figures/figS1.pdf", 
        width=8, height=5, units="in")


###############################################################################
# Plot correlation between genetic distance and physical position of SNPs
###############################################################################

maps<-readRDS("input/car001_d010.Rds")
head(maps)

maps %>% group_by(chromosome) %>% summarise(n=n(), 
                                            min_g=min(genetic_map),
                                            max_g=max(genetic_map), 
                                            min_p=min(physical_map),
                                            max_p=max(physical_map),
                                            corr=cor(genetic_map, physical_map))

plot4<-ggplot(data=maps, aes(x=genetic_map,y=physical_map)) + 
    facet_wrap(~chromosome, ncol = 3) +
    geom_point(size = 2, color = "cornflowerblue") +
    theme_bw() +
    xlab(paste0("Genetic distance (cM)")) + 
    ylab(paste0("Physical distance (Mb)")) +
    theme(axis.text.y = element_text(size = 11)) +
    theme(axis.text.x = element_text(size = 11)) +
    theme(axis.title.y = element_text(vjust = 4, size = 14)) +
    theme(axis.title.x = element_text(vjust = -1.5, size = 14)) +    
    theme(plot.margin = unit(c(0.2,0.2,0.5,0.5), "cm")) +
    theme(legend.position="none", strip.background=element_rect(colour="black", fill="white")) +    
    scale_y_continuous( breaks = seq(70000,51000000, 8000000), labels = as.character(seq(70000,51000000, 8000000)/1000000))

plot4

# save plot
ggsave(plot4,
       file = "figures/figS2.png", 
       width=7, height=8, units="in")






