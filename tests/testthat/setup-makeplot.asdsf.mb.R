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

c1 <- make_rwty_chain(c(t1, t2, t1, t2, t3, t1, t2, t3, t1, t1))
c2 <- make_rwty_chain(c(t2, t3, t1, t1, t1, t2, t2, t1, t1, t1))
cl <- list(c1, c2)

#
# #
# burnin = 2
# window.size = 0.5
# window.number = 4
#
#
#
# debug(makeplot.asdsf.mb)
# makeplot.asdsf.mb(list(f1, f2), 0, 0.5, 2)
