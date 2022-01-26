setwd("~/pascal/Master/Master_publication/data")

library(plotrix)

usfs <- read.delim("Ht_usfs_all.txt", header=FALSE)
usfs_syn <- read.table("Ht_usfs_all_syn.txt", quote="\"", comment.char="")
usfs_nonsyn <- read.table("Ht_usfs_all_nonsyn.txt", quote="\"", comment.char="")
gap.barplot(y=usfs[,3], xaxlab =c("invariant", 1:18), ytics = c(0, 5000, 10000, 15000, 2272409), gap=c(20000,2270000), col=c(rep("grey", 19)), ylab = "", xlab = "", main = "SFS - north", xtics = c(1.5, 2:19), las=0)
axis.break(2, 20000, breakcol="white", style="gap")
axis.break(2, 20000*(1+0.01), breakcol="black", style="slash")
axis.break(4, 20000*(1+0.01), breakcol="black", style="slash")

vector <- NULL
vector <- usfs[1,3]
for (i in 1:nrow(usfs_syn)) {
  vector <- c(vector,usfs_syn[i,3])
  vector <- c(vector,usfs_nonsyn[i,3])
}

gap.barplot(y=vector, ytics = c(0, 5000, 10000, 2169996), gap=c(13000,2167500), col=c("grey", rep(c("#1f78b4", "#a6cee3"), 30)), ylab = "Frequency", xlab = "Number of derived alleles", main = "", las=0, xtics = c(0, 1.75,2.25, 3.75,4.25, 5.75,6.25, 7.75,8.25, 9.75,10.25, 11.75,12.25, 13.75,14.25, 15.75,16.25, 17.75,18.25, 19.75,20.25, 21.75,22.25, 23.75,24.25, 25.75,26.25, 27.75,28.25, 29.75,30.25, 31.75,32.25, 33.75,34.25,35.75,36.25, 37.75,38.25, 39.75,40.25, 41.75,42.25, 43.75,44.25, 45.75,46.25, 47.75,48.25, 49.75,50.25, 51.75,52.25,53.75,54.25, 55.75,56.25,57.75,58.25, 59.75,60.25), xaxt='n')
axis(1,at=c(seq(4,60,4)), tck =0, labels = c(seq(2,30,2)), col = NA)
axis(1,at=0, tck =0, labels = "invariant", col = NA, las = 2)
axis.break(2, 13000, breakcol="white", style="gap")
axis.break(2, 13000*(1+0.01), breakcol="black", style="slash")
axis.break(4, 13000*(1+0.01), breakcol="black", style="slash")
legend("topright",legend=c("synonymous", "nonsynonymous"),col = c("#1f78b4", "#a6cee3"),bty = "n", fill = c("#1f78b4", "#a6cee3"), cex = 1)

usfs <- read.delim("HtN_usfs_all.txt", header=FALSE)
usfs_syn <- read.table("HtN_usfs_syn.txt", quote="\"", comment.char="")
usfs_nonsyn <- read.table("HtN_usfs_nonsyn.txt", quote="\"", comment.char="")

gap.barplot(y=vector, ytics = c(0, 5000, 10000, 2182499), gap=c(15000,2180000), col=c("grey", rep(c("#1f78b4", "#a6cee3"), 18)), ylab = "Frequency", xlab = "Number of derived alleles", main = "", las=0, xtics = c(0, 1.75,2.25, 3.75,4.25, 5.75,6.25, 7.75,8.25, 9.75,10.25, 11.75,12.25, 13.75,14.25, 15.75,16.25, 17.75,18.25, 19.75,20.25, 21.75,22.25, 23.75,24.25, 25.75,26.25, 27.75,28.25, 29.75,30.25, 31.75,32.25, 33.75,34.25,35.75,36.25), xaxt='n')
axis(1,at=c(seq(2,36,2)), tck =0, labels = c(1:18), col = NA)
axis(1,at=0, tck =0, labels = "invariant", col = NA, las = 2)
axis.break(2, 15000, breakcol="white", style="gap")
axis.break(2, 15000*(1+0.01), breakcol="black", style="slash")
axis.break(4, 15000*(1+0.01), breakcol="black", style="slash")
legend("topright",legend=c("synonymous", "nonsynonymous"),col = c("#1f78b4", "#a6cee3"),bty = "n", fill = c("#1f78b4", "#a6cee3"), cex = 1)

usfs <- read.delim("HtS_usfs_all.txt", header=FALSE)
usfs_syn <- read.table("HtS_usfs_syn.txt", quote="\"", comment.char="")
usfs_nonsyn <- read.table("HtS_usfs_nonsyn.txt", quote="\"", comment.char="")

gap.barplot(y=vector, ytics = c(0, 5000, 10000, 2176447), gap=c(15000,2174000), col=c("grey", rep(c("#1f78b4", "#a6cee3"), 18)), ylab = "Frequency", xlab = "Number of derived alleles", main = "", las=0, xtics = c(0, 1.75,2.25, 3.75,4.25, 5.75,6.25, 7.75,8.25, 9.75,10.25, 11.75,12.25, 13.75,14.25, 15.75,16.25, 17.75,18.25, 19.75,20.25, 21.75,22.25, 23.75,24.25), xaxt='n')
axis(1,at=c(seq(2,24,2)), tck =0, labels = c(1:12), col = NA)
axis(1,at=0, tck =0, labels = "invariant", col = NA, las = 2)
axis.break(2, 15000, breakcol="white", style="gap")
axis.break(2, 15000*(1+0.01), breakcol="black", style="slash")
axis.break(4, 15000*(1+0.01), breakcol="black", style="slash")
legend(x = 13, y=18000,legend=c("synonymous", "nonsynonymous"),col = c("#1f78b4", "#a6cee3"),bty = "n", fill = c("#1f78b4", "#a6cee3"), cex = 1)

# magnivora
usfs <- read.table("Hm_usfs_all.txt", quote="\"", comment.char="")
usfs_syn <- read.table("Hm_usfs_all_syn.txt", quote="\"", comment.char="")
usfs_nonsyn <- read.table("Hm_usfs_all_nonsyn.txt", quote="\"", comment.char="")
gap.barplot(y=vector, ytics = c(0, 5000, 10000, 2176675), gap=c(15000,2174500), col=c("grey", rep(c("orangered4", "orangered"), 9)), ylab = "Frequency", xlab = "Number of derived alleles", main = "", las=0, xtics = c(0, 1.75,2.25, 3.75,4.25, 5.75,6.25, 7.75,8.25, 9.75,10.25, 11.75,12.25), xaxt='n')
axis(1,at=c(seq(2,18,2)), tck =0, labels = c(1:9), col = NA)
axis(1,at=0, tck =0, labels = "invariant", col = NA, las = 2)
axis.break(2, 15000, breakcol="white", style="gap")
axis.break(2, 15000*(1+0.01), breakcol="black", style="slash")
axis.break(4, 15000*(1+0.01), breakcol="black", style="slash")
legend(x = 7, y=18000,legend=c("synonymous", "nonsynonymous"),col = c("orangered4", "orangered"),bty = "n", fill = c("orangered4", "orangered"), cex=1)

