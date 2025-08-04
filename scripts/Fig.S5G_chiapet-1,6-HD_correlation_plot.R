library(gtools)
library(pheatmap)

exprName='ChIAPET_PC3DM_HDtreat'
exprName='ChIAPET_HDtreat'
mapsize='1000000'
bin='10k'
rangekb=as.numeric(mapsize)/1000
fname=paste0('range',rangekb,'k.bin',bin)
datasuf=paste0('.maxbp',mapsize,'.bin',bin,'.txt')

glen=read.table('hg38.chrom.size', row.names=1)
chroms = paste0('chr', c(1:22, 'X')) 
nchr = length(chroms)

weights = glen[chroms,1]
weights = weights/sum(weights)

libnames=paste0('ACP00',c(10:13, 16:23))
nlib = length(libnames)
pairnames = NULL; 
sccs = NULL;
for ( i in 1:(nlib-1)){
    for (j in (i+1):nlib ){
        pname =  paste0(libnames[i],'_v_',libnames[j])
        f = paste0('scc.',libnames[i],'.',libnames[j], datasuf)
        sccs = c(sccs, sum(read.table(f)[1:nchr,1] * weights))
        pairnames = c(pairnames, pname)
    }
}

names(sccs) = pairnames
write.table(sccs, col.names=F, row.names=T, sep='\t', quote=F, file=paste0('weightedAverage_scc.',fname,'.',exprName,'.txt'))

cormat = matrix(1, ncol=nlib, nrow=nlib)
colnames(cormat) = rownames(cormat) = libnames;
k = 1
for ( i in 1:(nlib-1) ){
    for ( j in (i+1):nlib ){
        cormat[i,j] = cormat[j,i] = sccs[k]
        k=k+1
    }
}

pdf(paste0('weightedAverage_scc.',fname,'.',exprName,'test.pdf'))
pheatmap(cormat, display_numbers=F, border_color = '#c0c0c0c0', fontsize_number=12, fontsize=13, main=exprName)
dev.off()

# with numbers
pdf(paste0('weightedAverage_scc.',fname,'.',exprName,'.numbers.6-6.pdf'), width=6, height=6)
pheatmap(cormat, 
        display_numbers=T, 
        border_color = '#c0c0c0c0', 
        cluster_cols = FALSE,
        cluster_rows = FALSE,
        # number_color = "black",
        fontsize_number=8, 
        fontsize=13, 
        main=exprName)
dev.off()



