test_that("compare with makeplot.asdsf", {
    # asdsf.mb using 100% window sizes should compare produce the same output
    # as asdsf provided the window indices line up

    # burnin = 2 so 8 samples, 4 evenly spaced window endpoints
    obs1 <- makeplot.asdsf.mb(cl, 2, 1, 4)$asdsf.plot$data$ASDSF
    exp1 <- makeplot.asdsf(cl, 2, 2)$asdsf.plot$data$ASDSF
    expect_equal(obs1, exp1)

    # burnin = 1 so 9 samples, 9 evenly spaced window endpoints
    obs2 <- makeplot.asdsf.mb(cl, 1, 1, 9)$asdsf.plot$data$ASDSF
    exp2 <- makeplot.asdsf(cl, 1, 1)$asdsf.plot$data$ASDSF
    expect_equal(obs2, exp2)
})

test_that("check get.slide.freq.table.mb", {
    # burnin = 2, lookback = 75%, 2 evenly spaced window endpoints
    burnin <- 2
    window.size <- 0.75
    window.number <- 2
    gens.per.tree <- 1

    obs1 <- get.slide.freq.table.mb(
        c1$trees,
        burnin,
        window.size,
        window.number,
        gens.per.tree
    )
    obs2 <- get.slide.freq.table.mb(
        c2$trees,
        burnin,
        window.size,
        window.number,
        gens.per.tree
    )

    expect_equal(obs1$translation[-1, 2:3], obs2$translation[, 2:3])

    n.trees <- length(c1$trees) - burnin
    u <- round(seq_len(window.number) * n.trees / window.number)
    l <- round(u * (1 - window.size)) + 1

    # clade indicators for each tree
    # (1) 12, (2) 123, (3) 12345, (4) 125, (5) 134, (6) 145
    w1 <- c(0, 1, 1, 0, 0, 1)
    w2 <- c(0, 0, 1, 1, 1, 0)
    w3 <- c(1, 0, 1, 1, 0, 0)
    w <- cbind(w1, w2, w3)

    # tree counts for each interval, v{chain}_{interval}
    v1_1 <- c(1, 1, 1)
    v1_2 <- c(3, 1, 2)
    v1 <- cbind(v1_1, v1_2)

    v2_1 <- c(2, 1, 0)
    v2_2 <- c(4, 2, 0)
    v2 <- cbind(v2_1, v2_2)

    exp1 <- data.frame(w %*% scale(v1, FALSE, colSums(v1)))
    rownames(exp1) <- sapply(seq_len(nrow(exp1)), as.character)
    colnames(exp1) <- prettyNum((u + burnin) * gens.per.tree, sci = TRUE)
    expect_equal(obs1$slide.table, exp1)

    # split (1) isn't observed so we have to remove it
    m2 <- w %*% scale(v2, FALSE, colSums(v2))
    exp2 <- data.frame(m2[rowSums(m2) > 0, ])
    rownames(exp2) <- sapply(seq_len(nrow(exp2)), as.character)
    colnames(exp2) <- prettyNum((u + burnin) * gens.per.tree, sci = TRUE)
    expect_equal(obs2$slide.table, exp2)
})

test_that("check get.asdsfs.mb", {
    # burnin = 1, lookback = 2/3, 3 evenly spaced window endpoints
    burnin <- 1
    window.size <- 2/3
    window.number <- 3
    gens.per.tree <- 1

    obs_asdsf <- makeplot.asdsf.mb(
        cl,
        burnin,
        window.size,
        window.number
    )$asdsf.plot$data$ASDSF

    # check slide frequencies
    obs_sf1 <- get.slide.freq.table.mb(
        c1$trees,
        burnin,
        window.size,
        window.number,
        gens.per.tree
    )
    obs_sf2 <- get.slide.freq.table.mb(
        c2$trees,
        burnin,
        window.size,
        window.number,
        gens.per.tree
    )

    n.trees <- length(c1$trees) - burnin
    u <- round(seq_len(window.number) * n.trees / window.number)
    l <- round(u * (1 - window.size)) + 1

    # clade indicators for each tree
    # (1) 12, (2) 123, (3) 12345, (4) 125, (5) 134, (6) 145
    w1 <- c(0, 1, 1, 0, 0, 1)
    w2 <- c(0, 0, 1, 1, 1, 0)
    w3 <- c(1, 0, 1, 1, 0, 0)
    w <- cbind(w1, w2, w3)

    # tree counts for each interval, v{chain}_{interval}
    v1_1 <- c(1, 1, 0)
    v1_2 <- c(1, 2, 1)
    v1_3 <- c(3, 1, 2)
    v1 <- cbind(v1_1, v1_2, v1_3)
    u1 <- scale(v1, FALSE, colSums(v1))

    v2_1 <- c(2, 0, 0)
    v2_2 <- c(2, 2, 0)
    v2_3 <- c(4, 2, 0)
    v2 <- cbind(v2_1, v2_2, v2_3)
    u2 <- scale(v2, FALSE, colSums(v2))

    # checking frequencies
    exp_st1 <- data.frame(w %*% u1)
    rownames(exp_st1) <- sapply(seq_len(nrow(exp_st1)), as.character)
    colnames(exp_st1) <- prettyNum((u + burnin) * gens.per.tree, sci = TRUE)
    expect_equal(obs_sf1$slide.table, exp_st1)

    m2 <- w %*% scale(v2, FALSE, colSums(v2))
    exp_st2 <- data.frame(m2[rowSums(m2) > 0, ])
    rownames(exp_st2) <- sapply(seq_len(nrow(exp_st2)), as.character)
    colnames(exp_st2) <- prettyNum((u + burnin) * gens.per.tree, sci = TRUE)
    expect_equal(obs_sf2$slide.table, exp_st2)

    # manually checking asdsf
    exp_asdsf <- double(length(u))
    sp1 <- obs_sf1$translation[, "Tip numbers"]
    sp2 <- obs_sf2$translation[, "Tip numbers"]
    for (i in seq_along(exp_asdsf)) {
        sf1 <- obs_sf1$slide.table[[i]]
        sf2 <- obs_sf2$slide.table[[i]]

        sp_i <- unique(c(sp1[sf1 > 0], sp2[sf2 > 0]))
        asdsf_ik <- double(length(sp_i))

        for (k in seq_along(sp_i)) {
            f1_k <- sf1[sp1 == sp_i[k]]
            f2_k <- sf2[sp2 == sp_i[k]]

            g1_k <- ifelse(length(f1_k) == 0, 0, f1_k)
            g2_k <- ifelse(length(f2_k) == 0, 0, f2_k)

            asdsf_ik[k] <- sd(c(g1_k, g2_k))
        }
        exp_asdsf[i] <- mean(asdsf_ik)
    }
    expect_equal(obs_asdsf, exp_asdsf)
})
