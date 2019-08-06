library(optparse)
# make options
opt_list = list(
    make_option(c('-e', '--expr'), help='expr matrix file'),
    make_option(c('-g', '--group'), help='sample group file'),
    make_option(c('-c', '--contrast'), help='comparison info file'),
    make_option(c('-f', '--fold_change'), default=2.0, help='fold change cutoff', type='numeric'),
    make_option(c('-s', '--stat_cutoff'), default=0.05, help='pvalue or adjust pvalue cutoff', type='numeric'),
    make_option(c('-t', '--type'), default='padjust', help='use uncorrected pvalue if set to pvalue'),
    make_option(c('-o', '--out_prefix'), default='', help='out prefix')
)
# make options
opt_parser = OptionParser(option_list=opt_list)
opt = parse_args(opt_parser)
opt$f = log2(opt$f)
if (is.null(opt$e)){
    stop("expr file must be provided!")
}
if (is.null(opt$g)){
    stop("group file must be provided!")
}
if (is.null(opt$c)){
    stop('contrast info file must be provided!')
}


suppressMessages(library(DESeq2))
counts <- read.table(opt$e, header=T, row.names=1, sep="\t")
# counts = floor(counts+0.5)
group <- read.table(opt$g, header=T, row.names=1, sep='\t')
counts <- counts[, rownames(group)]
colnames(group)[1] = 'group'
## Calculation for P50O vs P50M
dds <- DESeqDataSetFromMatrix(countData=counts, colData=group, design= ~group)
rlogCounts = rlog(dds, blind=T)
# pdf(paste(opt$o, "PCA.pdf", sep=''))
# plotPCA(rlogCounts, 'group')
# dev.off()
write.table(assay(rlogCounts), paste(opt$o, "rlogCounts.matrix.xls", sep=''), sep='\t', quote=F, col.names = NA)
dds <- DESeq(dds)

compare <- read.table(opt$c, header=F, sep='\t')
for (i in 1:nrow(compare)){
    ctrl = as.character(compare[i, 1])
    test = as.character(compare[i, 2])
    ctrl_num = sum(group$group==ctrl)
    test_num = sum(group$group==test)
    print(paste(ctrl,'(', ctrl_num, ')', ' vs ', test, '(', test_num, ')', sep=''))
    res <- results(dds, contrast<-c('group', test, ctrl), pAdjustMethod="BH", parallel=T, alpha=opt$s, lfcThreshold=opt$f)
    summary(res)
    if (opt$t == 'padjust'){
        res$significant = (abs(res$log2FoldChange) >= opt$f) & (res[, 'padj'] <= opt$s)
        res[order(res$padj), ]
        print(paste('cutoff: |log2fc| >= ', opt$f, " & pajust", ' <= ', opt$s, sep=''))
    }else{
        res$significant = (abs(res$log2FoldChange) >= opt$f) & (res[, 'P.Value'] <= opt$s)
        res[order(res[, 'pvalue']), ]
        print(paste('cutoff: |log2fc| >= ', opt$f, " & pvalue", ' <= ', opt$s, sep=''))
    }
    res$regulate = 'up'
    res[res$log2FoldChange < 0, 'regulate'] = 'down'
    out_name = paste(opt$o, ctrl, '_vs_', test, ".deseq2.xls", sep='')
    write.table(res, file=out_name, sep="\t", quote=F, row.names=T, col.names=NA)
    ##根据p值筛选
    out_name = paste(opt$o, ctrl, '_vs_', test, '.DE.list', sep='')
    deg = rownames(res[(res$significant==T) %in% T, ])
    print(paste('DEG number: ', length(deg), sep=''))
    deg_reg = cbind(deg, res[(res$significant==T) %in% T, 'regulate'])
    write.table(deg_reg, out_name, sep='\t', col.names=F, quote=FALSE, row.names=F)
}
