setwd('~/Documents/roman/acute_eto_study_metabolomics_010415')

ProcessDataset <- function(dat)
{
    temp <- read.csv(dat$filename, as.is = TRUE, fill = TRUE,
                     check.names = FALSE, row.names = 1)
    samples <- rownames(temp)

    ## Remove metabolites with no values
    bad.mets <- sapply(temp, function(met) all(is.na(met)))
    temp <- temp[, !bad.mets]

    ## Impute missing values with minimum value for that metabolite
    temp <- as.data.frame(lapply(temp, Hmisc::impute, min))

    ## Normalize values by normalization control
    ctrl <- which(names(temp) == dat$norm.control)
    temp2 <- as.data.frame(apply(temp, 1, function(row) {
        .row <- as.numeric(row)
        .row / .row[ctrl]
    }))
    rownames(temp2) <- gsub('(.+?)\\.Results', '\\1', names(temp))
    names(temp2) <- samples

    ## Remove normalization control from expression values
    temp2 <- temp2[-ctrl, ]

    ## Log2 transform expression values for subsequent analysis using limma
    list(exprs = log2(temp2),
         groups = plyr::mapvalues(
            gsub('A.+?-(.+)', '\\1', samples),
            from = c('C', '40', '60'),
            to = c('control', 'low', 'high'), warn_missing = FALSE
            ))
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
        res <- dplyr::arrange(res, pval)
        names(res)[-1] <- paste(contrast, names(res)[-1], sep = '.')

        res
    }, simplify = FALSE)
}

Boxplots <- function(dat, features, center = 'control')
{
    # Center by median of control group
    med <- apply(dat$exprs, 1, function(row) {
        tapply(row, dat$groups, FUN = 'median')[center]
    })
    dat$exprs <- dat$exprs - med

    require(reshape2)
    .dat <- melt(
        data.frame(feature = rownames(dat$exprs), dat$exprs,
                   stringsAsFactors = FALSE),
        id.vars = 'feature'
        )

    require(dplyr)
    .dat <- transmute(.dat,
        feature,
        sample = variable,
        group = rep(dat$groups, each = nrow(dat$exprs)),
        expr = value
        )
    .dat <- filter(.dat, feature %in% features)

    require(ggplot2)
    update_geom_defaults('point', list(colour = NULL))
    plot <-
        ggplot(data = .dat,
               aes(x = factor(feature), y = expr, color = group)) +
        geom_boxplot(position = position_dodge(width = 0.7), width = 0.7,
                     outlier.colour = NULL) +
        scale_color_brewer(palette = 'Set1') +
        coord_flip() +
        ylab('Median-centered Log2 Abundance') +
        theme_bw() +
        theme(axis.title.x = element_text(size = 12, face = 'bold',
                                          vjust = 0),
              axis.text.x = element_text(size = 12),
              axis.title.y = element_blank(),
              axis.text.y = element_text(size = 12, angle = 30, hjust = 1,
                                         vjust = 0),
              axis.ticks.y = element_blank(),
              legend.position = 'none')
    update_geom_defaults('point', list(colour = 'black'))

    print(plot)
}


##
## Main
##

## Polar metabolites

## Read in raw metabolite expression data and log2 transform
raw.files <- list(
    pos = list(filename = 'acuteetostudy_polar_pos.csv',
               norm.control = 'D3.N15.Serine.Results'),
    neg = list(filename = 'acuteetostudy_polar_neg.csv',
               norm.control = 'd3.serine.IS.Results')
    )
dat <- lapply(raw.files, ProcessDataset)

## Differential expression analysis
de <- mapply(function(d, name) {
    res <- DiffEq(d$exprs, d$groups, name,
                  low = low - control,
                  high = high - control,
                  dose = high - low)
}, dat, names(dat), SIMPLIFY = FALSE)

## Merge contrast data
de <- lapply(de, function(x) {
    Reduce(function(x, y) merge(x, y, by = 'feature', all = TRUE), x)
})

## Merge in raw expression values into DE dataframes
merged <- mapply(function(a, b) {
    merge(a, b$exprs, by.x = 'feature', by.y = 0)
}, de, dat, SIMPLIFY = FALSE)

## Write out DE data
mapply(write.csv, merged, paste0(names(merged), '_metabolites_de_res.csv'),
       MoreArgs = list(row.names = FALSE))

## Draw boxplots
features <- c(
    'Acetyl.CoA',
    'alpha.ketoglutarate',
    'ATP',
    'Citrate',
    'fumarate',
    'GTP',
    'isocitrate',
    'malate',
    'NAD',
    'NADH',
    'oxaloacetate',
    'succinate',
    'transaconitate'
    )
pdf('tca_cycle_boxplots.pdf')
dat$neg$groups <- factor(dat$neg$groups, levels = c('high', 'low', 'control'))
Boxplots(dat$neg, features)
dev.off()

## Non-polar metabolites

## Read in raw metabolite expression data and log2 transform
raw.files <- list(
    pos = list(filename = 'acuteetostudy_nonpolar_pos.csv',
               norm.control = 'C12.MAGE.Results'),
    neg = list(filename = 'acuteetostudy_nonpolar_neg.csv',
               norm.control = 'PDA.Results')
    )
dat <- lapply(raw.files, ProcessDataset)

## Differential expression analysis
de <- mapply(function(d, name) {
    res <- DiffEq(d$exprs, d$groups, name,
                  low = low - control,
                  high = high - control,
                  dose = high - low)
}, dat, names(dat), SIMPLIFY = FALSE)

## Merge contrast data
de <- lapply(de, function(x) {
    Reduce(function(x, y) merge(x, y, by = 'feature', all = TRUE), x)
})

## Merge in raw expression values into DE dataframes
merged <- mapply(function(a, b) {
    merge(a, b$exprs, by.x = 'feature', by.y = 0)
}, de, dat, SIMPLIFY = FALSE)

## Write out DE data
mapply(write.csv, merged, paste0(names(merged), '_np_metabolites_de_res.csv'),
       MoreArgs = list(row.names = FALSE))
