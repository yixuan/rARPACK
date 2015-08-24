library(Matrix)
library(rARPACK)
library(microbenchmark)
library(dplyr)
library(ggplot2)

m = 200
n = 100
k = 10

set.seed(123)
x100 = matrix(rnorm(n^2), n)
x100[sample(n^2, floor(n^2 / 2))] = 0
xsp100 = as(x100, "dgCMatrix")
y100 = crossprod(x100)
m100 = matrix(rnorm(m*n), m)
mt100 = t(m100)

m = 2000
n = 1000
x1000 = matrix(rnorm(n^2), n)
x1000[sample(n^2, floor(n^2 / 2))] = 0
xsp1000 = as(x1000, "dgCMatrix")
y1000 = crossprod(x1000)
m1000 = matrix(rnorm(m*n), m)
mt1000 = t(m1000)

nrun = 50

packageVersion("rARPACK")

res_eigs = microbenchmark(
    "eigs[100]"             = eigs(x100, k),
    "eigs[100, sparse]"     = eigs(xsp100, k),
    "eigs[100, symmetric]"  = eigs_sym(y100, k),
    "eigs[1000]"            = eigs(x1000, k),
    "eigs[1000, sparse]"    = eigs(xsp1000, k),
    "eigs[1000, symmetric]" = eigs_sym(y1000, k),
    times = nrun,
    control = list(order = "inorder")
)
write.csv(res_eigs, "eigs_0.8-0.csv")

res_svds = microbenchmark(
    "svds[100]"             = svds(m100, k),
    "sleep1"                = Sys.sleep(0.01),
    "svds[100, transpose]"  = svds(mt100, k),
    "sleep2"                = Sys.sleep(0.01),
    "svds[1000]"            = svds(m1000, k),
    "sleep3"                = Sys.sleep(0.01),
    "svds[1000, transpose]" = svds(mt1000, k),
    times = nrun,
    control = list(order = "inorder")
)
write.csv(res_svds, "svds_0.8-0.csv")

## Install version 0.7-0 and rerun the benchmarks above
write.csv(res_eigs, "eigs_0.7-0.csv")
write.csv(res_svds, "svds_0.7-0.csv")



## Compare results
res_eigs_07 = read.csv("eigs_0.7-0.csv")
res_eigs_08 = read.csv("eigs_0.8-0.csv")
res_svds_07 = read.csv("svds_0.7-0.csv")
res_svds_08 = read.csv("svds_0.8-0.csv")

res_eigs = rbind(
    res_eigs_07 %>% mutate(version = "0.7-0"),
    res_eigs_08 %>% mutate(version = "0.8-0")
)

res_svds = rbind(
    res_svds_07 %>% mutate(version = "0.7-0"),
    res_svds_08 %>% mutate(version = "0.8-0")
)

## Visualize results
visualize = function(res)
{
    dat = res %>% filter(!grepl("sleep", expr)) %>%
        mutate(size = ifelse(grepl("1000", expr), "1000", "100")) %>%
        group_by(expr, version, size) %>%
        summarize(medtime = median(time) / 1e6)
    
    pdat = res %>% filter(!grepl("sleep", expr)) %>%
        mutate(size = ifelse(grepl("1000", expr), "1000", "100"))
    
    g = ggplot(dat, aes(x = expr, y = medtime)) +
        geom_bar(aes(color = version), fill = "white",
                 stat = "identity", position = "dodge") +
        geom_jitter(aes(x = expr, y = time / 1e6, color = version), pdat) +
        facet_wrap(~ size, ncol = 2, scales = "free") +
        scale_y_continuous("Elapsed time (milliseconds)") +
        scale_x_discrete("Functions") +
        theme_grey(base_size = 18)
    print(g)
}

visualize(res_eigs)
visualize(res_svds)
