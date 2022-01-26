library(vcfR)
library(adegenet)
library(hierfstat)
library(geodist)
library(ecodist)
library(adespatial)
library(vegan)

### genetic differenciation using Euclidian distance:

# load coordinates of Ht
xyHam <- read.delim("~/bioinformatics/GeographicPhylogeny/coords.txt")
names(xyHam) <- c("ID", "latitude", "longitude")
xyHam <- xyHam[-c(1:3),]

# load SNP data
vcf <- read.vcfR("hamiltosporidium_10122020.SNP.tvaer.recode.vcf", verbose = FALSE)

# convert data to genind
genind_H <- vcfR2genind(vcf)
# add pop labels
pop(genind_H)<-substr(indNames(genind_H),1,5)
# convert to genpop object to calculate distance
genpop_H <- genind2genpop(genind_H, pop= indNames(genind_H))
distgenEUCL <- dist(genpop_H, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)

CoAncestH <- as.matrix(distgenEUCL)

# mantel tests
ecodist::mantel(as.dist(CoAncestH)~as.dist(geodist(xyHam)),nperm = 100000)
#mantelr      pval1      pval2      pval3  llim.2.5% ulim.97.5%
#  0.4121365  0.0082600  0.9917500  0.0082600  0.3557073  0.5139311
ecodist::mantel(as.dist(CoAncestH_N)~as.dist(geodist(xyHam_N)),nperm = 100000)
#mantelr      pval1      pval2      pval3  llim.2.5% ulim.97.5%
#  0.7394963  0.0008400  0.9991700  0.0008400  0.5740292  0.9005247
ecodist::mantel(as.dist(CoAncestH_S)~as.dist(geodist(xyHam_S)),nperm = 100000)
#mantelr      pval1      pval2      pval3  llim.2.5% ulim.97.5%
#  0.9773313  0.0014400  1.0000000  0.0014400  0.9557155  0.9991905

# dbMEM by RDA
CoAncestH.pc <- prcomp(as.dist(CoAncestH)) # r-stats
xyHam.dbmem <- dbmem(as.dist(geodist(xyHam)))
Ham.rda<-rda(CoAncestH.pc$x, xyHam.dbmem)
RsquareAdj(Ham.rda)
#$r.squared
#[1] 0.1524055
#
#$adj.r.squared
#[1] 0.08720588
anova(Ham.rda, perm=1000)
#Permutation test for rda under reduced model
#Permutation: free
#Number of permutations: 999
#
#Model: rda(X = CoAncestH.pc$x, Y = xyHam.dbmem)
#Df Variance      F Pr(>F)
#Model     1   329440 2.3375  0.197
#Residual 13  1832163

CoAncestH_N.pc <- prcomp(as.dist(CoAncestH_N)) # r-stats
xyHam_N.dbmem <- dbmem(as.dist(geodist(xyHam_N)))
Ham_N.rda<-rda(CoAncestH_N.pc$x, xyHam_N.dbmem)
RsquareAdj(Ham_N.rda)
#$r.squared
#[1] 0.4332768
#
#$adj.r.squared
#[1] 0.3523163
anova(Ham_N.rda, perm=1000)
#Permutation test for rda under reduced model
#Permutation: free
#Number of permutations: 999
#
#Model: rda(X = CoAncestH_N.pc$x, Y = xyHam_N.dbmem)
#         Df Variance      F Pr(>F)
#Model     1   115256 5.3517  0.001 ***
#Residual  7   150754
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

CoAncestH_S.pc <- prcomp(as.dist(CoAncestH_S)) # r-stats
xyHam_S.dbmem <- dbmem(as.dist(geodist(xyHam_S)))
Ham_S.rda<-rda(CoAncestH_S.pc$x, xyHam_S.dbmem)
RsquareAdj(Ham_S.rda)
#$r.squared
#[1] 0.8231823
#
#$adj.r.squared
#[1] 0.7789778
anova(Ham_S.rda, perm=1000)
#'nperm' >= set of all permutations: complete enumeration.
#Set of permutations < 'minperm'. Generating entire set.
#Permutation test for rda under reduced model
#Permutation: free
#Number of permutations: 719
#
#Model: rda(X = CoAncestH_S.pc$x, Y = xyHam_S.dbmem)
#Df Variance      F   Pr(>F)
#Model     1   661448 18.622 0.002778 **
#  Residual  4   142077
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
