# preparation
library(dplyr)
library(ggplot2)
library(dendextend)
library(tidyverse)  
library(cluster)    
library(factoextra)
library(readxl)
library(sysfonts)
library(showtextdb)
library(showtext)
library(plotrix)
library(RColorBrewer)
library(readxl)
library(GOplot)
library(ggdendro)
library(gridExtra)
library(ggrepel)

# plot setting
myfree<-theme_set(theme_bw())
windowsFonts(HEL=windowsFont("Helvetica CE 55 Roman"),
             RMN=windowsFont("Times New Roman"),
             ARL=windowsFont("Arial"))
old_theme <- theme_update(
  plot.title=element_text(family="RMN", size=18, colour="black"),
  legend.title=element_text(family = "RMN", size = 18, colour = "black"),
  legend.text=element_text(family = "RMN", size = 18, colour = "black"),
  legend.key.size=unit(.2,"inches"),
  axis.title.x=element_text(family="RMN", size=18, colour="black"),
  axis.title.y=element_text(family="RMN", size=18,angle=90, colour="black"),
  axis.text.x=element_text(family="RMN", size=14),
  axis.text.y=element_text(family="RMN", size=14),
  axis.ticks=element_line(colour="black"),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  panel.background=element_blank(),
  axis.line=element_line(size=1)
)

############## hierarchical clustering analysis (HCA) ########################

##################### fecundity ###################

# DATA PREPARATION 

# data import
Fecundity <- read_excel("C:/Users/Yahan/Desktop/LR/SturgeonDataSet_YH_0604.xlsx", 
                                      sheet = "Fecundity", col_names = TRUE)
# convert data from tibble into a data frame
Fecundity <- as.data.frame(Fecundity)
rownames(Fecundity) <- Fecundity$Species
# tiding data
# remove useless columns
Fecundity <- Fecundity[,-c(4:6)]
Fecundity <- Fecundity[,-1]
# clear NA
Fecundity = na.omit(Fecundity) 
# scaling/standardizing the data
Fecundity <- scale(Fecundity)
head(Fecundity)

# Agglomerative Hierarchical Clustering
# methods to assess
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")

# function to compute coefficient
ac <- function(m) {
  agnes(Fecundity, method = m)$ac
}
map_dbl(m, ac)
# based on the results, choose ward as method
# visualize the dendrogram
# Hierarchical clustering using Ward
hc_Fe <- agnes(Fecundity, method = "ward")
pltree(hc_Fe, cex = 0.6, hang = -1, main = "Dendrogram of agnes for Fecundity")

# Divisive Hierarchical Clustering
# Cut tree into 3 groups
sub_grp_Fe <- cutree(hc_Fe, k = 3)
# Number of members in each cluster
table(sub_grp_Fe)
# specify the border colors 
rect.hclust(hc_Fe, k = 3, border = 2:5)
# visualize the result in a scatter plot
fviz_cluster(list(data = Fecundity, cluster = sub_grp_Fe))

##################### growth ###################

# DATA PREPARATION 

# data import
Growth_o <- read_excel("C:/Users/Yahan/Desktop/LR/SturgeonDataSet_YH_0604.xlsx", 
                        sheet = "Growth(VGBF)", col_names = TRUE)
# convert data from tibble into a data frame
Growth_o <- as.data.frame(Growth_o)
rownames(Growth_o) <- Growth_o$Species
# tiding data
# remove useless columns
Growth_o <- Growth_o[,-c(4:6)]
Growth_o <- Growth_o[,-1]
# clear NA
Growth_o = na.omit(Growth_o) 
# scaling/standardizing the data
Growth <- scale(Growth_o)
head(Growth)

# Agglomerative Hierarchical Clustering
# function to compute coefficient
ac <- function(m) {
  agnes(Growth, method = m)$ac
}
map_dbl(m, ac)
# based on the results, choose ward as method
# visualize the dendrogram
# Hierarchical clustering using Ward
hc_Gr <- agnes(Growth, method = "ward")
pltree(hc_Gr, cex = 0.6, hang = -1, main = "Dendrogram of agnes for Growth")

# Divisive Hierarchical Clustering
# Cut tree into 3 groups
sub_grp_Gr <- cutree(hc_Gr, k = 4)
# Number of members in each cluster
table(sub_grp_Gr)
# specify the border colors 
rect.hclust(hc_Gr, k = 4, border = 2:5)
# visualize the result in a scatter plot
fviz_cluster(list(data = Growth, cluster = sub_grp_Gr))

# there is a trend for k value and infinite body length, so I calculate the lg value and plot
Growth_o$LgK <- log(Growth_o$K)
Growth_o$LgMaxLen <- log(Growth_o$MaxLen)
# sum of squared rediduals 
# select the best line going through the cloud of points
m1 <- lm(LgK ~ MaxLen, data = Growth_o) # to fit the linear model 
summary(m1) # show all the info that we need 
# prediction and prediction error
# create a scatterplot with the least squares line laid on top
plot(Growth_o$LgK ~ Growth_o$MaxLen, 
     xlab = "Maximum body length/cm", ylab = "Lg K", 
     family="RMN", cex.lab =1.5, cex.main = 2,cex = 2,
     main = "A). Relationship for all sturgeons" ) 
abline(m1) # plot a line based on its slope and intercept
anova(m1)
par(mfrow=c(2,2))
plot(m1)

# extract huso 
huso_Gr <- as.data.frame(Growth_o[c(51,52,53,54,46,56,57,41,36),])
huso_Gr <- huso_Gr[order(huso_Gr$MaxLen), ]
m3 <- lm(K~MaxLen +I(MaxLen^(2)), data = huso_Gr)
summary(m3)
fitted(m3)
residuals(m3)
plot(x = huso_Gr$MaxLen, y = huso_Gr$K, 
     xlab = "Maximum body length/cm", ylab = "K",
     family="RMN", cex.lab =1.5, cex.main = 2,cex = 2,
     main = "C). Relationship for Huso genus")
lines(huso_Gr$MaxLen, fitted(m3))
par(mfrow=c(2,2))
plot(m3) 

# extract A. oxyrinchus
Aoxy_Gr <- as.data.frame(Growth_o[c(20,21,27,32,38,40,45,49),])
m4 <- lm(LgK ~ MaxLen, data = Aoxy_Gr) # to fit the linear model 
summary(m4) # show all the info that we need 
# prediction and prediction error
# create a scatterplot with the least squares line laid on top
plot(Aoxy_Gr$LgK ~ Aoxy_Gr$MaxLen, 
           xlab = "Maximum body length/cm", ylab = "Lg K", 
           family="RMN", cex.lab =1.5, cex.main = 2,cex = 2,
           main = "D). Relationship for A.oxyrinchus") 
abline(m4) # plot a line based on its slope and intercept
anova(m4)
par(mfrow=c(2,2))
plot(m4)

# other sturgeons
St_Gr <- as.data.frame(Growth_o[-c(51,52,53,54,46,56,57,41,36,20,21,27,32,38,40,45,49,1,7),])

m2 <- lm(K~MaxLen +I(MaxLen^(2)), data = St_Gr)
summary(m2)
fitted(m2)
residuals(m2)
plot(x = St_Gr$MaxLen, y = St_Gr$K, 
           xlab = "Maximum body length/cm", ylab = "K",
           family="RMN", cex.lab =1.5, cex.main = 2, cex = 2,
           main = "B). Improved for rest majority")
lines(St_Gr$MaxLen, fitted(m2))
par(mfrow=c(2,2))
plot(m2) 

##################### maturity ###################

# DATA PREPARATION 

# data import
Maturity <- read_excel("C:/Users/Yahan/Desktop/LR/SturgeonDataSet_YH_0604.xlsx", 
                         sheet = "Maturity", col_names = TRUE)
# convert data from tibble into a data frame
Maturity <- as.data.frame(Maturity)
rownames(Maturity) <- Maturity$Species
# tiding data
# remove useless columns
Maturity <- Maturity[,-c(6:10)]
Maturity <- Maturity[,-1]
# clear NA
Maturity = na.omit(Maturity) 
# scaling/standardizing the data
Maturity <- scale(Maturity)
head(Maturity)

# Agglomerative Hierarchical Clustering
# function to compute coefficient
ac <- function(m) {
  agnes(Maturity, method = m)$ac
}
map_dbl(m, ac)
# based on the results, choose ward as method
# visualize the dendrogram
# Hierarchical clustering using Ward
hc_Mt <- agnes(Maturity, method = "ward")
pltree(hc_Mt, cex = 0.6, hang = -1, main = "Dendrogram of agnes for Maturity")

# Divisive Hierarchical Clustering
# Cut tree into 3 groups
sub_grp_Mt <- cutree(hc_Mt, k = 4)
# Number of members in each cluster
table(sub_grp_Mt)
# specify the border colors 
rect.hclust(hc_Mt, k = 4, border = 2:5)
# visualize the result in a scatter plot
fviz_cluster(list(data = Maturity, cluster = sub_grp_Mt))


######################## MaxLen vs Fecundity ################################

# data import
MaxLen_Fe <- read_excel("C:/Users/Yahan/Desktop/LR/SturgeonDataSet_YH_0604.xlsx", 
                                      sheet = "MaxLen_Fe")
# data tiding
# convert into data frame
MaxLen_Fe <- as.data.frame(MaxLen_Fe)
# set rownames
row.names(MaxLen_Fe) <- MaxLen_Fe$Species
# remove useluess columns 
MaxLen_Fe <- MaxLen_Fe[, -1]
# remove NAs
MaxLen_Fe = na.omit(MaxLen_Fe)
head(MaxLen_Fe)

# visualization Max_Fe and MaxLen
plot(x = MaxLen_Fe$MaxLen, y = MaxLen_Fe$MaxFe)
cor(MaxLen_Fe$MaxLen, MaxLen_Fe$MaxFe)
# sum of squared rediduals 
# select the best line going through the cloud of points
Max1 <- lm(MaxFe ~ MaxLen, data = MaxLen_Fe) # to fit the linear model 
summary(Max1) # show all the info that we need 
# prediction and prediction error
# create a scatterplot with the least squares line laid on top
plot(x = MaxLen_Fe$MaxLen, y = MaxLen_Fe$MaxFe)
abline(Max1) #plot a line based on its slope and intercept
anova(Max1)
par(mfrow=c(2,2))
plot(Max1)

# visualization LgMax_Fe and MaxLen
MaxLen_Fe$LgMaxFe <- log(MaxLen_Fe$MaxFe)
plot(x = MaxLen_Fe$MaxLen, y = MaxLen_Fe$LgMaxFe)
cor(MaxLen_Fe$MaxLen, MaxLen_Fe$LgMaxFe)
# sum of squared rediduals 
# select the best line going through the cloud of points
Max2 <- lm(LgMaxFe ~ MaxLen, data = MaxLen_Fe) # to fit the linear model 
summary(Max2) # show all the info that we need 
# prediction and prediction error
# create a scatterplot with the least squares line laid on top
plot(x = MaxLen_Fe$MaxLen, y = MaxLen_Fe$LgMaxFe)
abline(Max2) #plot a line based on its slope and intercept
anova(Max2)
par(mfrow=c(2,2))
plot(Max2)
# 
Max3 <- lm(LgMaxFe~MaxLen +I(MaxLen^2), data = MaxLen_Fe)
summary(Max3)
fitted(Max3)
residuals(Max3)
plot(x = MaxLen_Fe$MaxLen, y = MaxLen_Fe$LgMaxFe, 
     xlab = "Maximum body length", ylab = "Lg value of Maximum Fecundity",
     family="RMN", cex.lab =1.5, cex.main = 1.5, cex = 2,
     main = "A). Relationship between fecundity and max body length")
lines(MaxLen_Fe$MaxLen, fitted(Max3))
par(mfrow=c(2,2))
plot(Max3) 

# exclude huso 
MaxLen_Fe <- MaxLen_Fe[-15,]

# visualization LgMax_Fe and MaxLen
plot(x = MaxLen_Fe$MaxLen, y = MaxLen_Fe$LgMaxFe)
cor(MaxLen_Fe$MaxLen, MaxLen_Fe$LgMaxFe)
# sum of squared rediduals 
# select the best line going through the cloud of points
Max2 <- lm(LgMaxFe ~ MaxLen, data = MaxLen_Fe) # to fit the linear model 
summary(Max2) # show all the info that we need 
# prediction and prediction error
# create a scatterplot with the least squares line laid on top
plot(x = MaxLen_Fe$MaxLen, y = MaxLen_Fe$LgMaxFe)
abline(Max2) #plot a line based on its slope and intercept
anova(Max2)
par(mfrow=c(2,2))
plot(Max2)
# 
Max3 <- lm(LgMaxFe~MaxLen +I(MaxLen^2), data = MaxLen_Fe)
summary(Max3)
fitted(Max3)
residuals(Max3)
plot(x = MaxLen_Fe1$MaxLen, y = MaxLen_Fe1$LgMaxFe, 
     xlab = "Maximum body length", ylab = "Lg value of Maximum Fecundity",
     family="RMN", cex.lab =1.5, cex.main = 1.5, cex = 2,
     main = "B). Improved model for fecundity and max body length")
lines(MaxLen_Fe$MaxLen, fitted(Max3))
par(mfrow=c(2,2))
plot(Max3) 

############ MAX SIZE VS IUCN ############################

# data import
MAX_SIZE_IUCN <- read_excel("C:/Users/Yahan/Desktop/LR/SturgeonDataSet_YH_0604.xlsx", 
                            sheet = "MaxSize_IUCN")
# data tiding
# clear useless columns
MAX_SIZE_IUCN <- MAX_SIZE_IUCN[,-c(6:7)]
# convert into data frame
MAX_SIZE_IUCN <- as.data.frame(MAX_SIZE_IUCN)
# remove the NAs
MAX_SIZE_IUCN = na.omit(MAX_SIZE_IUCN) 
# set rownames for MaxLen and IUCN table
row.names(MAX_SIZE_IUCN) <- MAX_SIZE_IUCN$Species
MAX_SIZE_IUCN = MAX_SIZE_IUCN[,-1] 
head(MAX_SIZE_IUCN)

# visualization
# prepare the colour in ggplot
colourCount = length(unique(MAX_SIZE_IUCN$Names))
getPalette = colorRampPalette(brewer.pal(9, "Paired"))
# plot
ggplot(MAX_SIZE_IUCN, aes(x = factor(MAX_SIZE_IUCN$IUCN, levels = c("CR", "EN", "VU","NT","LC")), y = Leng_max)) +
  stat_boxplot(geom = "errorbar",width=0.15,aes(color="black"))+ 
  geom_boxplot(size=0.5,fill="white",outlier.fill="white",outlier.color="white")+ 
  geom_jitter(aes(fill = Names),width =0.2,shape = 21, size=2.5)+ 
  scale_fill_manual(values = getPalette(colourCount))+
  scale_color_manual(values=c("black","black","black","black","black"))+ 
  theme_update()+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  ylab("Maximum Size/cm")+xlab("IUCN Red List Category") 


##################### MaxLen vs Age range at first maturity ###################

# data import
AgeRange <- read_excel("C:/Users/Yahan/Desktop/LR/SturgeonDataSet_YH_0604.xlsx", 
                       sheet = "MaxLen_AgeRange")

# visualization
plot(AgeRange$Leng_max,AgeRange$Age, 
     xlab = "Maximum Body Length/cm",
     ylab = "Age Range at the first maturity of Sturgeons/years",
     cex = 0, 
     col = getPalette(colourCount))
i=seq(length(AgeRange$Leng_max))
segments(AgeRange$Leng_max[i], AgeRange$Age[i],
         AgeRange$Leng_max[i+17], AgeRange$Age[i+17], 
         col = getPalette(colourCount))
text(Age ~Leng_max, labels=abnm, data=AgeRange, cex=0.8, font=1, col = getPalette(colourCount))


################### MaxLen vs EcoThreats ######################

# data import
# import data containing the degree of ecological threats for each sturgeon species
EcoT <- read_excel("C:/Users/Yahan/Desktop/LR/SturgeonDataSet_YH_0604.xlsx", 
                   sheet = "MaxLen_EcoThreats")

# visualization
# set maximum body length as y axise, threats as x axise
ggplot(EcoT, aes(x = ecothr, y = maxlen), group = continent)+
  geom_point(aes(cex = grade, color = sturgeon),alpha=0.6) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10)) +
  labs(title= '', x= 'Ecological Threats', y= 'Maximum Body Length') +
  facet_grid(.~continent, scales="free_x") +
  geom_text_repel(
    aes(label = name), 
    size = 3,
    segment.color = "black") +
  theme(legend.position = "none") 

