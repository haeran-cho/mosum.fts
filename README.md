# mosum.fts
Multiple change point detection under a static factor model

Codes accompanying the paper "Moving sum procedure for multiple change point detection in large factor models" by Barigozzi, Cho and Trapani (2024).

## file description

### codes.R
Contains the main routine `mosum.fts` for the proposed MOSUM procedure for time series segmentation under a factor model.

## usage example

```
source('codes.R')

set.seed(1234)
dd <- duan_dgp(n = 1000, p = 500, type = c(2, 1, 4), dep = TRUE)
G0 <- round(dim(dd$x)[2]^(max(2 / 5, 1 - min(1, log(dim(dd$x)[1])/log(dim(dd$x)[2])))) * log(dim(dd$x)[2])^1.1)
out <- mosum.fts(x = dd$x, G = G0)
out$cpts
plot(out$stat, type = 'l')
abline(h = out$thr, col = 3)
abline(v = out$cpts, col = 2); abline(v = dd$k0, col = 4)
```
