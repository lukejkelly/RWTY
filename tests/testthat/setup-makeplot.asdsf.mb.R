make_rwty_chain <- function(trees) {
    f <- list(trees = trees, ptable = NULL, gens.per.tree = 1)
    class(f) <- "rwty.chain"
    return(f)
}

# abc and ade
t1 <- read.tree(text = "(a:2,(b:1,c:1):1,(d:1,e:1):1);", keep.multi = TRUE)
# abe and acd
t2 <- read.tree(text = "((d:1,c:1):1,(b:1,e:1):1,a:2);", keep.multi = TRUE)
# abe and ab
t3 <- read.tree(text = "((d:1,c:1):1,(a:1,b:1):1,e:2);", keep.multi = TRUE)

# par(mfrow = c(1, 3))
# purrr::walk(c(t1, t2, t3), ~ plot.phylo(., "unrooted"))

# 10 samples each
b1 <- c(t1, t2, t1, t2, t3, t1, t2, t3, t1, t1)
b2 <- c(t2, t3, t1, t1, t1, t2, t2, t1, t1, t1)

c1 <- make_rwty_chain(b1)
c2 <- make_rwty_chain(b2)
cl <- list(c1, c2)

# sufficiently long chains for debugging with makeplot.asdsf
d1 <- make_rwty_chain(c(b1, b1))
d2 <- make_rwty_chain(c(b2, b2))
dl <- list(d1, d2)
