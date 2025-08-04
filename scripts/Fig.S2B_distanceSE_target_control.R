#SuperEnhancer dots distances to target: SEL1 or control: WNT4 
#Compare these 2 distributions
#Ways:
#1a. Count how many around within a certain radius
#1b. Calculate mean distance within a certain radius
#2a. Get the minimum distance SEL1 or WNT4 to superEnhancer dots
#2b. Calculate the average distance of N-closest 

# targetName <- 'SELENBP1'
# contrlName <- 'WNT4'
# targetName <- 'CD36'
# contrlName <- 'FAM20C'
targetName <- 'VMP1'
contrlName <- 'NoGene'

library(dplyr)
targetcolor='#03a5fc'
ctrlcolor='grey'

#load data from  csvs
ftsv <- paste0("./data/misc/Imaging_analysis/",targetName, '_', contrlName, '_SEdist.micron.tsv')
dist.df <- read.delim(ftsv, stringsAsFactors = F);
hist(dist.df$dist) #Check distribution of all distances

dth=1.2 # threshold in micron
d.contrl = subset(dist.df, dist < dth & Gene == contrlName)$dist
d.target = subset(dist.df, dist < dth & Gene == targetName)$dist

wilcox.test(d.target, d.contrl, alternative = 'less')

pdf(paste0(targetName,"_vs_",contrlName,"_SE_density_plot.pdf"), width=6, height=5)
plot(density(d.contrl), col=ctrlcolor, main='SE density', xlab='Distance to gene (nm)', lwd=3, lty=2)
lines(density(d.target), col=targetcolor, lwd=3)
legend('topleft', bty='n', lwd=c(3,3), col=c(targetcolor, ctrlcolor), lty=c(1,2), legend=c(targetName, contrlName))
dev.off()


