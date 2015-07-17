Fisher <- function(M, N, k, r)
{
    sum(sapply(r:k, dhyper, M, N, k))
}

GeneSet <- function()
{
    require(org.Hs.eg.db)
    require(GO.db)

    go2gene <- as.list(org.Hs.egGO2ALLEGS)
    go2gene <- go2gene['GO:0006631']

    genes <- unique(unlist(go2gene))
    gsc <- piano::loadGSC(
        data.frame(gene = genes,
                   set = 'fatty acid metabolic process',
                   stringsAsFactors = FALSE)
        )

    list(genes = genes, gsc = gsc)
}

GSA <- function(deList, gsc)
{
    require(piano)
    lapply(deList, function(de) {
        stats <- de$t
        names(stats) <- de$gene
        GSAsummaryTable(runGSA(stats, geneSetStat = 'maxmean', gsc = gsc))
    })
}

StratifyByReceptorStatus <- function(eSet, er, pr, her2, pos, neg)
{
    require(Biobase)

    ## Define receptor status classes
    pheno <- pData(eSet)
    pData(eSet)$receptor_status <- factor(
        ifelse(
            pheno[, er] == pos |
            pheno[, pr] == pos |
            pheno[, her2] == pos, 'RP',
            ifelse(
                pheno[, er] == neg &
                pheno[, pr] == neg &
                pheno[, her2] == neg, 'TN',
                NA
                )
            ),
        levels = c('RP', 'TN')
        )

    ## Remove samples without receptor status
    eSet[, !is.na(pData(eSet)$receptor_status)]
}

DrawHeatmap <- function(data, groups, centroid.group = 'TN', ...)
{
    require(cluster)
    require(gplots)
    require(RColorBrewer)

    ## Prepare expression matrix
    d <- t(scale(t(data.matrix(data[complete.cases(data), ]))))

    ## Calculate the centroid for TN tumors
    centroid <- rowMeans(d[, groups == centroid.group])

    unique.cols <- brewer.pal(8, 'Dark2')[seq(levels(groups))]
    col.cols <- as.character(plyr::mapvalues(groups,
                                             from = levels(groups),
                                             to = unique.cols))

    breaks <- c(seq(-sqrt(2), -0.1, length = 100),
                seq(-0.1, -0.1, length = 1),
                seq(0.1, sqrt(2), length = 100))
    colors <- bluered(length(breaks) - 1)

    ## cluster samples by pearson distance and average linkage
    hc <- as.dendrogram(
        agnes(as.dist(1 - cor(d, method = 'pearson')), method = 'average')
        )

    ## Determine column order
    require(dendextend)
    corr <- -apply(d, 2, cor, y = centroid, method = 'pearson')
    sample.order <- colnames(d)[rev(order(corr))]
    hc <- rotate(hc, sample.order)

    # heatmap.2(d, Colv = hc, Rowv = FALSE, dendrogram = 'col',
    #           scale = 'none', breaks = breaks, density.info = 'none',
    #           trace = 'none', col = colors, keysize = 1,
    #           ColSideColors = col.cols, labCol = groups, ...)

    list(centroid = centroid, dendrogram = hc)
}

SampleCorrelations <- function(eSet, centroid)
{
    ## Find common subset of genes to calculate correlation
    features <- intersect(featureNames(eSet), names(centroid))
    eSet <- eSet[features, ]
    centroid <- centroid[features]

    ## Calculate sample-wise correlation of dataset to centroid
    d <- t(scale(t(data.matrix(exprs(eSet)))))
    corr <- apply(d, 2, cor.test, y = centroid, method = 'pearson')
    corr.est <- sapply(corr, '[[', 'estimate')
    corr.pval <- p.adjust(sapply(corr, '[[', 'p.value'), method = 'fdr')

    ## Minimum significant positive correlation
    min.pos.corr <- min(corr.est[corr.est > 0 & corr.pval < 0.05])

    ## t-test to test for difference in correlation between groups
    t.res <- t.test(corr.est ~ pData(eSet)$receptor_status)$p.value

    ## Perform Fisher's Exact Test for coefficients >= 0.1 and 0.3
    tn.sample <- pData(eSet)$receptor_status == 'TN'
    fisher.res <- sapply(c(0.1, 0.3), function(coef) {
        corr.filter <- corr.est >= coef & corr.pval < 0.05

        M <- sum(tn.sample)
        N <- sum(!tn.sample)
        k <- sum(corr.filter)
        r <- sum(tn.sample & corr.filter)
        r2 <- sum(!tn.sample & corr.filter)

        list(pval = Fisher(M, N, k, r),
             TN = r,
             TN.perc = r / M,
             RP = r2,
             RP.perc = r2 / N,
             corr.TN.perc = r / k,
             corr.RP.perc = r2 / k)
    })

    ## Order samples by correlation to centroid
    pData(eSet)$corr.est <- corr.est
    eSet <- eSet[, order(pData(eSet)$corr.est)]

    # Draw plots
    require(ggplot2)
    plot <-
        ggplot(data = pData(eSet),
               aes(x = seq(corr.est), y = corr.est, fill = receptor_status)) +
        geom_bar(stat = 'identity', width = 0.7) +
        scale_fill_brewer(palette = 'Dark2') +
        theme_bw() +
        theme(legend.position = 'none',
              axis.line = element_line(size = 0.7),
              axis.line.x = element_blank(),
              panel.border = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              axis.title = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank())

    print(plot)
}

DiffEq <- function(eSet, by, ...) {
    ## Make a design matrix
    attr <- if (length(by) == 1) Biobase::pData(eSet)[, by] else by
    design <- model.matrix(~0 + attr)
    colnames(design) <- gsub('attr', '', colnames(design))

    ## Differential expression
    fit <- limma::lmFit(Biobase::exprs(eSet), design = design)
    contrasts <- limma::makeContrasts(..., levels = design)
    fit2 <- limma::contrasts.fit(fit, contrasts = contrasts)
    fit2 <- limma::eBayes(fit2)

    ## Construct list of results
    sapply(colnames(contrasts), function(contrast) {
        res <- data.frame(
            gene = rownames(fit2),
            logFC = fit2$coefficients[, contrast],
            t = fit2$t[, contrast],
            pval = fit2$p.value[, contrast],
            fdr = p.adjust(fit2$p.value[, contrast], 'fdr'),
            stringsAsFactors = FALSE
            )
        dplyr::arrange(res, pval)
    }, simplify = FALSE)
}

require(Biobase)
setClass('ExpressionSet2', contains = 'ExpressionSet',
    representation = representation(
        timeCol = 'character',
        eventCol = 'character'
        )
    )
setGeneric('ExpressionSet2',
    function(eSet, timeCol = '', eventCol = '') {
        standardGeneric('ExpressionSet2')
    }
    )
setMethod('ExpressionSet2',
    signature = signature(eSet = 'ExpressionSet'),
    function(eSet, timeCol, eventCol) {
        new('ExpressionSet2',
            assayData = Biobase::assayData(eSet),
            phenoData = Biobase::phenoData(eSet),
            featureData = Biobase::featureData(eSet),
            experimentData = Biobase::experimentData(eSet),
            annotation = Biobase::annotation(eSet),
            protocolData = Biobase::protocolData(eSet),
            timeCol = timeCol,
            eventCol = eventCol
            )
    }
    )

setGeneric('Time', function(eSet) standardGeneric('Time'))
setMethod('Time',
    signature = signature(eSet = 'ExpressionSet2'),
    function(eSet) {
        if (eSet@timeCol == '')
            NULL
        else
            Biobase::pData(eSet)[, eSet@timeCol]
    }
    )

setGeneric('Event', function(eSet) standardGeneric('Event'))
setMethod('Event',
    signature = signature(eSet = 'ExpressionSet2'),
    function(eSet) {
        if (eSet@eventCol == '')
            NULL
        else
            Biobase::pData(eSet)[, eSet@eventCol]
    }
    )

GroupedUnivariateAnalysis <- function(eSet, genes) {
    .genes <- intersect(featureNames(eSet), genes)
    plyr::rbind.fill(lapply(.genes, function(g) {
        cbind(Gene = g,
              UnivariateAnalysis(exprs(eSet)[g, ],
                                 pData(eSet)$receptor_status,
                                 Time(eSet),
                                 Event(eSet)))
    }))
}

UnivariateAnalysis <- function(exprs, group, time, event) {
    require(survival)

    valid <- !is.na(group)
    .exprs <- scale(exprs[valid])
    .group <- group[valid]
    .time <- time[valid]
    .event <- event[valid]

    SubsetAnalysis <- function(subset) {
        ..exprs <- .exprs[subset]
        ..time <- .time[subset]
        ..event <- .event[subset]

        reg <- coxph(Surv(..time, ..event) ~ ..exprs)
        summ <- summary(reg)
        conf.int <- summ$conf.int
        p <- summ$logtest['pvalue']

        data.frame(N = length(..exprs),
                   HazardRatio = sprintf('%.3f (%.3f - %.3f)',
                                         conf.int[, 1],
                                         conf.int[, 3],
                                         conf.int[, 4]),
                   LogRank.P = p,
                   stringsAsFactors = FALSE)
    }

    subsets <- cbind(
        rep(T, length(.exprs)),
        sapply(levels(.group), function(g) .group == g))

    res <- apply(subsets, 2, SubsetAnalysis)
    res <- plyr::rbind.fill(res)
    res <- cbind(Subset = c('All', levels(.group)), res)

    res
}

DiscretizeSamplesByQuantile <- function(exprs, time, event, .quantile) {
    rank <- round(.quantile * length(exprs) / 100)
    groups <- DiscretizeSamplesByRank(exprs, rank)
    reg <- coxph(Surv(time, event) ~ groups)

    list(P = summary(reg)$sctest['pvalue'],
         quantile = .quantile,
         groups = groups)
}

DiscretizeSamplesByRank <- function(exprs, rank) {
    .exprs <- exprs[order(exprs)]
    cutoff <- .exprs[rank]

    DiscretizeSamples(exprs, cutoff)
}

DiscretizeSamples <- function(exprs, cutoff) {
    ifelse(exprs < cutoff, 1, 2)
}

GetOptimalQuantile <- function(gene.exprs, time, event, min.quantile = 0.1) {
    start <- ceiling(min.quantile * length(gene.exprs))
    end <- floor((1 - min.quantile) * length(gene.exprs))

    p.dist <- lapply(start:end, function(rank) {
        groups <- DiscretizeSamplesByRank(gene.exprs, rank)
        reg <- coxph(Surv(time, event) ~ groups)

        list(P = summary(reg)$sctest['pvalue'],
             quantile = rank / length(gene.exprs) * 100,
             groups = groups)
    })

    optimal <- which.min(sapply(p.dist, function(x) x$P))
    p.dist[[optimal]]
}

DrawKMPlot <- function(eSet, gene, title.prefix, filename,
                       .quantile = NULL, ...) {
    require(survival)

    temp <- if (is.null(.quantile)) {
        GetOptimalQuantile(exprs(eSet)[gene, ], Time(eSet), Event(eSet))
    } else {
        DiscretizeSamplesByQuantile(exprs(eSet)[gene, ], Time(eSet),
                                    Event(eSet), .quantile)
    }
    fit <- survfit(Surv(Time(eSet), Event(eSet)) ~ temp$groups)

    pdf(filename)
    col <- c('black', 'red')
    title <- sprintf('%s (n = %i)\nLog-rank p = %.3e',
                     title.prefix, length(temp$groups), temp$P)
    par(mar = c(5, 6, 4, 2))
    plot(fit, col = col, main = title, cex.lab = 2,
         cex.axis = 1.7, cex.main = 1.3, las = 1, lwd = 3, font.lab = 2, ...)
    legend('bottomright',
           legend = c(
                sprintf('< %.2f percentile (%i)', temp$quantile,
                        sum(temp$groups == 1)),
                sprintf('>= %.2f percentile (%i)', temp$quantile,
                        sum(temp$groups == 2))),
           col = col,
           lty = c(1, 1),
           lwd = c(5, 5),
           cex = 1.7)
    dev.off()

    return(fit)
}


##
## Main
##

## Source datasets
library(brca)

data(tcga.brca)
data(ispy1)
data(gse25066)
data(yau)
data(chin)

## Stratify by receptor status
tcga.brca <- StratifyByReceptorStatus(tcga.brca,
                                      er = 'ER_Status_nature2012',
                                      pr = 'PR_Status_nature2012',
                                      her2 = 'HER2_Final_Status_nature2012',
                                      pos = 'Positive',
                                      neg = 'Negative')

ispy1 <- StratifyByReceptorStatus(ispy1,
                                  er = 'er',
                                  pr = 'pr',
                                  her2 = 'her2',
                                  pos = 1,
                                  neg = 0)

gse25066 <- StratifyByReceptorStatus(gse25066,
                                     er = 'er',
                                     pr = 'pr',
                                     her2 = 'her2',
                                     pos = 'P',
                                     neg = 'N')

yau <- StratifyByReceptorStatus(yau,
                                er = 'er',
                                pr = 'er',
                                her2 = 'her2',
                                pos = 1,
                                neg = 0)

pData(chin)$erb_ihc <- plyr::mapvalues(pData(chin)$erb_ihc,
                                       from = 0:1,
                                       to = c('Negative', 'Positive'))
chin <- StratifyByReceptorStatus(chin,
                                 er = 'er',
                                 pr = 'pr',
                                 her2 = 'erb_ihc',
                                 pos = 'Positive',
                                 neg = 'Negative')

## Count samples
sample.breakdown <- lapply(
    list(tcga = tcga.brca,
         ispy1 = ispy1,
         gse25066 = gse25066,
         yau = yau,
         chin = chin),
    function(eSet) {
        breakdown = table(pData(eSet)$receptor_status)
        rbind(as.matrix(breakdown), sum(breakdown))
    }
    )

## TCGA differential expression analysis
res <- DiffEq(tcga.brca, by = 'receptor_status', TN - RP)[[1]]

require(dplyr)

## Format and write out TCGA DE data
res2 <-
    merge(res, pData(featureData(tcga.brca)),
          by.x = 'gene', by.y = 'entrez', all.x = TRUE) %>%
    dplyr::select(symbol, logFC, t, pval, fdr) %>%
    arrange(pval)
write.csv(res2, file = 'tcga_receptor_status_de.csv', row.names = FALSE)

## Extract genes involved in fatty acid metabolic process
gene.set <- GeneSet()
goi <-
    filter(res, gene %in% gene.set$genes)  %>%
    arrange(desc(logFC))

## Format and write out DE data
goi2 <-
    merge(goi, pData(featureData(tcga.brca)),
          by.x = 'gene', by.y = 'entrez', all.x = TRUE) %>%
    dplyr::select(symbol, logFC, t, pval, fdr) %>%
    arrange(pval)
write.csv(goi2, file = 'fat_met_genes_TN_vs_RP.csv', row.names = FALSE)

## Gene set enrichment analysis
gsa.res <- GSA(
    list(tcga = res,
         ispy1 = DiffEq(ispy1, by = 'receptor_status', TN - RP)[[1]],
         gse25066 = DiffEq(gse25066, by = 'receptor_status', TN - RP)[[1]],
         yau = DiffEq(yau, by = 'receptor_status', TN - RP)[[1]],
         chin = DiffEq(chin, by = 'receptor_status', TN - RP)[[1]]),
    gsc = gene.set$gsc
    )

# Draw heatmap for genes of interest
jpeg('fat_met_genes_TN_vs_RP.jpg', width = 10, height = 10,
     units = 'in', res = 300)
heat <- DrawHeatmap(
    exprs(tcga.brca)[goi$gene, ],
    pData(tcga.brca)$receptor_status
    )
dev.off()

clades <- cutree(heat$dendrogram, 2)

jpeg('tcga_correlation.jpg', width = 18, height = 4, units = 'in', res = 600)
SampleCorrelations(tcga.brca, centroid = heat$centroid)
dev.off()

##
## ISPY1 dataset
##

## Correlation to TCGA TN centroid
jpeg('ispy1_correlation.jpg', width = 18, height = 4, units = 'in', res = 600)
SampleCorrelations(ispy1, centroid = heat$centroid)
dev.off()

##
## GSE25066 dataset
##

## Correlation to TCGA TN centroid
jpeg('gse25066_correlation.jpg', width = 18, height = 4, units = 'in',
     res = 600)
SampleCorrelations(gse25066, centroid = heat$centroid)
dev.off()

##
## YAU dataset
##

## Correlation to TCGA TN centroid
jpeg('YAU_correlation.jpg', width = 18, height = 4, units = 'in', res = 600)
SampleCorrelations(yau, centroid = heat$centroid)
dev.off()

##
## CHIN dataset
##

## Correlation to TCGA TN centroid
jpeg('chin_correlation.jpg', width = 18, height = 4, units = 'in', res = 600)
SampleCorrelations(chin, centroid = heat$centroid)
dev.off()

##
## Univariate Analyses
##

datasets <- list(
    ispy1     = ExpressionSet2(ispy1,
                               timeCol = 'rfs.t',
                               eventCol = 'rfs.e'),
    gse25066  = ExpressionSet2(gse25066,
                               timeCol = 'drfs_t',
                               eventCol = 'drfs_e'),
    yau       = ExpressionSet2(yau,
                               timeCol = 't_dmfs',
                               eventCol = 'e_dmfs'),
    chin      = ExpressionSet2(chin,
                               timeCol = 'disease_time',
                               eventCol = 'disease_binary')
    )
univariate <- lapply(datasets, GroupedUnivariateAnalysis, genes = goi$gene)
univariate <- plyr::rbind.fill(mapply(
    function(d, res) cbind(res, Dataset = d),
    names(univariate),
    univariate,
    SIMPLIFY = FALSE,
    USE.NAMES = FALSE
    ))

## Merge in symbols
gene.symbols <- read.table(
    '~/Documents/tcga/TCGA_BRCA_exp_HiSeqV2-2014-05-02/gene_map',
    sep = '|',
    as.is = TRUE,
    fill = TRUE)
univariate <-
    merge(univariate, gene.symbols, by.x = 'Gene', by.y = 'V2',
          all.x = TRUE) %>%
    dplyr::select(Gene = V1, Subset, N, HazardRatio, LogRank.P, Dataset) %>%
    arrange(Dataset, Gene, Subset)
write.csv(univariate, file = 'univariate_analysis.csv', row.names = FALSE)

##
## KM Plots
##

a <- DrawKMPlot(
    datasets$gse25066,
    '32',
    'Pooled Neoadjuvant Chemotherapy Treated\nAll Samples',
    'gse25066_ACACB_all.pdf',
    xlab = 'Time (years)')

b <- DrawKMPlot(
    datasets$gse25066[, pData(datasets$gse25066)$receptor_status == 'RP'],
    '32',
    'Pooled Neoadjuvant Chemotherapy Treated\nReceptor-Positive',
    'gse25066_ACACB_RP.pdf',
    xlab = 'Time (years)')

d <- DrawKMPlot(
    datasets$gse25066[, pData(datasets$gse25066)$receptor_status == 'TN'],
    '32',
    'Pooled Neoadjuvant Chemotherapy Treated\nTriple-Negative',
    'gse25066_ACACB_TN.pdf',
    xlab = 'Time (years)')
