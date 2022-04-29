#' Plot the Standard Deviation of Split Frequencies over the course of an MCMC
#' with windows growing in size in the style of MrBayes.
#'
#' This function takes two or more rwty.chain ojects and returns a plot of
#' ASDSF at number.points evenly spaced iterations.
#' Each ASDSF is calculated on the most recent window.size of samples
#' The solid line with points shows the Average Standard Deviation of Split Frequences at the current generation
#'
#' @param chains A list of rwty.chain objects.
#' @param burnin The number of trees to eliminate as burnin. Defaults to zero.
#' @param window.size The percentage of preceding samples to use for each estimate.
#' @param window.number The number of ASDSFs estimates to plot at evenly spaced number of iterations after burnin and including the final sample.
#' @param min.freq The minimum frequency for a node to be used for calculating ASDSF.
#' @param log.y Controls whether they Y axis is plotted on a log scale or not.  Which scale is more useful depends largely on the amount of disagreement between your chains.  Attempting to make an asdsf plot with a log Y axis for chains that include standard deviations of zero will result in warning messages.
#'
#' @return output A cumulative plot of ASDSF across all chains
#'
#' @keywords MCMC, phylogenetics, ASDSF, cumulative
#'
#' @export makeplot.asdsf.mb
#' @examples
#' \dontrun{
#'     data(fungus)
#'     p <- makeplot.asdsf.mb(fungus, burnin = 20)
#'     p
#' }

makeplot.asdsf.mb <- function(chains,
                              burnin = 0,
                              window.size = 0.75,
                              window.number = 10,
                              min.freq = 0.0,
                              log.y = TRUE) {

    message("Creating ASDSF.mb plot")

    chains <- check.chains(chains)
    # labels <- names(chains)
    slide.freq.list <- slide.freq.mb(chains, burnin, window.size, window.number)
    dat <- get.asdsfs.mb(slide.freq.list, min.freq)

    asdsf.plot <- ggplot(dat, aes(x = as.numeric(as.character(Generation)))) +
        geom_line(aes(y = ASDSF)) +
        geom_point(aes(y = ASDSF)) +
        geom_hline(aes(yintercept = 0.01), alpha = 0.25, linetype = "dashed") +
        theme(legend.position = "none") +
        theme_light() +
        xlab("Generation") +
        ylab("ASDSF") +
        xlim(0, NA) +
        ggtitle("Average Standard Deviation of Split Frequencies")

    if (log.y == TRUE) {
        asdsf.plot <- asdsf.plot + scale_y_log10()
    }

    return(list("asdsf.plot" = asdsf.plot))
}

# Calculating expanding sliding window frequencies
slide.freq.mb <- function(chains, burnin, window.size, window.number) {

    chains <- check.chains(chains)
    trees <- lapply(chains, function(x) x[['trees']])

    if (length(trees[[1]]) - burnin < window.number) {
        stop(
            "burnin is too large to make at least window.number points for the plot, quitting. Try setting a smaller burnin and/or window.number"
        )
    }

    slide.freq.list <- lapply(
        trees,
        get.slide.freq.table.mb,
        burnin = burnin,
        window.size = window.size,
        window.number = window.number,
        gens.per.tree = chains[[1]]$gens.per.tree
    )
    return(slide.freq.list)
}

# Get trees on each window
get.slide.freq.table.mb <- function(tree.list,
                                    burnin,
                                    window.size,
                                    window.number,
                                    gens.per.tree = 1) {
    # Instead of splitting trees disjointly, we create a list such that entry j
    # includes the most recent window.size proportion of trees up to the
    # iteration corresponding to the final entry in that window

    tree.list <- tree.list[(burnin + 1):length(tree.list)]

    n.tree <- length(tree.list)
    i.win <- seq_len(window.number)
    u <- round(n.tree * i.win / window.number)
    l <- round(u * (1 - window.size)) + 1

    # overlapping sliding windows of trees
    tree.windows <- purrr::map2(l, u, ~ tree.list[seq.int(.x, .y)])

    # clade frequencies for each list of trees
    clade.freq.list <- purrr::map(
        i.win,
        \(i) clade.freq(tree.windows[[i]], start = 1, end = u[i] - l[i] + 1)
    )

    # rename the list to the last tree index of each window
    names(clade.freq.list) <- prettyNum(
        (burnin + u) * gens.per.tree,
        sci = TRUE
    )

    # this is the table of frequencies in each window
    slide.freq.table <- clade.freq.list[[1]]
    colnames(slide.freq.table)[-1] <- names(clade.freq.list)[1]
    for (i in 2:length(clade.freq.list)) {
        slide.freq.table <- merge(
            slide.freq.table,
            clade.freq.list[[i]],
            by = "cladenames",
            all = TRUE
        )
        colnames(slide.freq.table)[
            which(colnames(slide.freq.table) == "cladefreqs")
        ] <- names(clade.freq.list)[i]
    }

    # tidy things up
    slide.freq.table[is.na(slide.freq.table)] <- 0.0
    rownames(slide.freq.table) <- slide.freq.table$cladenames
    slide.freq.table <- slide.freq.table[, -1]

    # # calculate sd and mean of cumulative frequency and mean
    # thissd <- apply(slide.freq.table, 1, sd)
    # thismean <- apply(slide.freq.table, 1, mean)
    # thisess <- apply(slide.freq.table, 1, effectiveSize)
    #
    # slide.freq.table$sd <- thissd
    # slide.freq.table$mean <- thismean
    # slide.freq.table$ess <- thisess
    #
    # # Sorting by sd, since these are usually the most interesting clades
    # slide.freq.table <- slide.freq.table[
    #     order(slide.freq.table$sd, decreasing = TRUE),
    # ]

    # Building a new table that contains parsed clade names
    translation.table <- cbind(
        as.numeric(as.factor(rownames(slide.freq.table))),
        as.character(rownames(slide.freq.table)),
        parse.slide.clades(rownames(slide.freq.table), tree.list)
    )
    colnames(translation.table) <- c("Clade number", "Tip numbers", "Tip names")

    # Setting slide.freq.table to the same names as the translation table
    rownames(slide.freq.table) <- as.numeric(
        as.factor(rownames(slide.freq.table))
    )

    output <- list(
        "slide.table" = slide.freq.table,
        "translation" = translation.table
    )
    class(output) <- "rwty.slide"
    return(output)
}


# Calculate ASDSF from sliding window split frequencies
get.asdsfs.mb <- function(slide.freq.list, min.freq) {

    x <- slide.freq.list

    # set initial values
    sets <- length(x)
    slide.wins <- ncol(x[[1]]$slide.table)

    # # remove mean and sd etc. columns
    # slide.wins <- length(
    #     setdiff(names(x[[1]]$slide.table), c("mean", "sd", "ess", "wcsf"))
    # )

    # use to label plot
    all_SDs <- list()

    # populate clade frequency tables, one for each window
    for (i in seq_len(slide.wins)) {
        winTable <- NULL

        # attach tip names to clade frequencies for the first chain (necessary
        # because clade name/numbers may not match among chains)
        winTable <- data.frame(
            cbind(x[[1]]$translation[, 3], x[[1]]$slide.table[, i]),
            stringsAsFactors = FALSE
        )

        colnames(winTable) <- c("Tip names", names(x)[1])
        class(winTable[, 2]) <- "numeric"
        # Populate the rest of the table, one chain at a time
        for (j in 2:sets) {
            thisTable <- data.frame(
                cbind(x[[j]]$translation[, 3], x[[j]]$slide.table[, i]),
                stringsAsFactors = FALSE
            )
            colnames(thisTable) <- c("Tip names", names(x)[j])
            winTable <- merge(winTable, thisTable, by = "Tip names", all = TRUE)
            winTable[is.na(winTable)] <- 0
            class(winTable[, j + 1]) <- "numeric"
        }
        # calculate ASDSF across all chains for each window
        # filter for clades below min freq then save ASDSF for that window
        winTable <- winTable[, -1]
        winTable <- winTable[
            apply(winTable, MARGIN = 1, function(x) any(x > min.freq)),
        ]
        all_SDs[[i]] <- apply(winTable, 1, sd, na.rm = TRUE)
    }

    # now that we have collected SDSF for all windows, return a data frame
    generations <- as.numeric(
        as.character(names(x[[1]]$slide.table[1:slide.wins]))
    )
    d <- data.frame(
        split.frequency = unlist(all_SDs),
        Generation = rep(
            generations[seq_along(all_SDs)],
            times = sapply(all_SDs, length)
        )
    )

    dat <- ddply(d, .(Generation), summarize, ASDSF = mean(split.frequency))
    return(dat)
}
