setwd('~/Documents/roman/mmtv_metabolites_de_analysis_061114')

KEGGPathways <- function()
{
    require(KEGGREST)

    pathways <- keggList('pathway', 'mmu')
    pathway.ids <- gsub('path:', '', names(pathways))
    pathway.names <- sprintf('%s (%s)',
                             gsub('(.*) - .+$', '\\1', pathways),
                             pathway.ids)

    pathway.members <- mapply(function(p, n) {
        tryCatch({
            temp <- keggGet(p)[[1]][['COMPOUND']]
            stopifnot(!is.null(temp))
            data.frame(compound = names(temp),
                       set = n,
                       stringsAsFactors = FALSE)
        }, error = function(err) {
            data.frame()
        })
    }, pathway.ids, pathway.names, SIMPLIFY = FALSE)

    piano::loadGSC(as.data.frame(data.table::rbindlist(pathway.members)))
}

DiffEq <- function(exprs, groups, desc, ...)
{
    require(limma)

    ## Make a design matrix
    design <- model.matrix(~0 + groups)
    colnames(design) <- gsub('groups', '', colnames(design))

    ## Differential expression
    fit <- lmFit(exprs, design = design)
    contrasts <- makeContrasts(..., levels = design)
    fit2 <- contrasts.fit(fit, contrasts = contrasts)
    fit2 <- eBayes(fit2)

    ## Construct list of results
    sapply(colnames(contrasts), function(contrast) {
        res <- data.frame(
            feature = rownames(fit2),
            logFC = fit2$coefficients[, contrast],
            t = fit2$t[, contrast],
            pval = fit2$p.value[, contrast],
            fdr = p.adjust(fit2$p.value[, contrast], 'fdr'),
            desc = desc,
            stringsAsFactors = FALSE
            )
        dplyr::arrange(res, pval)
    }, simplify = FALSE)
}

VolcanoPlot <- function(features, logfc, pval)
{
    groups <- mapply(function(a, b) {
        ifelse(b < 0.05,
               ifelse(abs(a) >= 2,
                      'A',
                      ifelse(abs(a) >= 1,
                             'B',
                             'C')),
               'D')
    }, logfc, pval)

    ## legend labels
    counts <- table(groups)
    labels <- c(sprintf('Log2 Fold Change >= 2\nRaw p-value < 0.05\n(N = %i)',
                        counts['A']),
                sprintf('Log2 Fold Change >= 1\nRaw p-value < 0.05\n(N = %i)',
                        counts['B']),
                sprintf('Log2 Fold Change < 1\nRaw p-value < 0.05\n(N = %i)',
                        counts['C']),
                sprintf('Raw p-value >= 0.05\n(N = %i)',
                        counts['D']))

    require(ggplot2)
    df <- data.frame(features, logfc, pval, groups, stringsAsFactors = FALSE)
    plot <-
        ggplot(data = df,
               aes(x = logfc, y = -log10(pval), fill = factor(groups),
                   size = abs(logfc))) +
        geom_vline(xintercept = c(-2, 2), linetype = 'longdash',
                   color = 'grey') +
        geom_vline(xintercept = c(-1, 1), linetype = 'longdash',
                   color = 'grey') +
        geom_hline(yintercept = -log10(0.05), linetype = 'longdash',
                   color = 'grey') +
        geom_point(alpha = 0.5, pch = 21, color = 'black') +
        scale_fill_brewer(name = 'Legend', labels = labels, palette = 'Set1') +
        scale_size_continuous(breaks = NULL) +
        xlab('Log2 Fold Change') + ylab('-log10(Raw p-value)') +
        theme_bw() +
        theme(axis.title.x = element_text(size = 12, face = 'bold'),
              axis.text.x = element_text(size = 12),
              axis.title.y = element_text(size = 12, face = 'bold'),
              axis.text.y = element_text(size = 12),
              legend.title = element_text(size = 12),
              legend.text = element_text(size = 10),
              legend.key.size = grid::unit(1.5, 'cm'))

    print(plot)
}

##
## Main
##

## Read in raw metabolite expression data and log2 transform
raw.files <- c(polar = 'polar_metabolites.csv',
               lipids = 'lipid_metabolites.csv')
dat <- lapply(raw.files, function(source) {
    raw <- read.csv(source, as.is = TRUE, fill = TRUE, check.names = FALSE,
                    row.names = 1)

    list(exprs = t(log2(data.matrix(raw[, -1]))),
         groups = raw[, 1])
})

## Convert polar metabolite names to KEGG IDs
kegg.ids <- read.csv('polar_metabolites_name_map.csv')
rownames(dat$polar$exprs) <- kegg.ids$KEGG

## Differential expression analysis
de <- mapply(function(d, name) {
    DiffEq(d$exprs, d$groups, name, Tumor - Control)[[1]]
}, dat, names(dat), SIMPLIFY = FALSE)

## KEGG pathway enrichment analysis for polar metabolites
require(piano)
stats <- de$polar$t
names(stats) <- de$polar$feature
polar.pa <- GSAsummaryTable(
    runGSA(stats, geneSetStat = 'maxmean', gsc = KEGGPathways(), nPerm = 1e6)
    )
polar.pa <- polar.pa[order(polar.pa[, 4]), ]
write.csv(polar.pa, file = 'polar_metabolites_pathway_analysis.csv',
          row.names = FALSE)

## Merge in raw expression values into DE dataframes
merged <- mapply(function(de, dat) {
    dplyr::arrange(merge(de, dat$exprs, by.x = 'feature', by.y = 0), pval)
}, de, dat, SIMPLIFY = FALSE)

## Convert KEGG IDs in DE result to human-readable names
merged$polar$feature <- kegg.ids[match(de$polar$feature, kegg.ids$KEGG),
                                 'Match']

## Write out DE data
merged <- as.data.frame(data.table::rbindlist(merged))
write.csv(merged, file = 'all_metabolites_de_res.csv', row.names = FALSE)

## Draw combined volcano plot
pdf('combined_volcano_plot.pdf')
VolcanoPlot(merged$feature, merged$logFC, merged$pval)
dev.off()
