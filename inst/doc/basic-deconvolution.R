## ------------------------------------------------------------------------
library('dtangle')
names(shen_orr_ex)

## ------------------------------------------------------------------------
truth = shen_orr_ex$annotation$mixture
head(truth)

## ----cache=FALSE---------------------------------------------------------
pure_samples <- lapply(1:3, function(i) {
    which(truth[, i] == 1)
})
names(pure_samples) = colnames(truth)
pure_samples

## ------------------------------------------------------------------------
Y <- shen_orr_ex$data$log
Y[1:4,1:4]

## ------------------------------------------------------------------------
marker_list = find_markers(Y,pure_samples,data_type="microarray-gene",marker_method='ratio')

## ------------------------------------------------------------------------
lapply(marker_list$L,head)

## ------------------------------------------------------------------------
lapply(marker_list$V,head)

## ------------------------------------------------------------------------
q = .1
quantiles = lapply(marker_list$V,function(x)quantile(x,1-q))
K = length(pure_samples)
n_choose = sapply(1:K,function(i){max(which(marker_list$V[[i]] > quantiles[[i]]))})
n_choose

## ------------------------------------------------------------------------
marks = marker_list$L
dc <- dtangle(Y,pure_samples,n_choose,data_type='microarray-gene',markers=marks)

## ------------------------------------------------------------------------
dc

## ----results='asis',fig.height=5,fig.width=5-----------------------------
phats <- dc$estimates
plot(truth,phats,xlab="Truth",ylab="Estimates",xlim=c(0,1),ylim=c(0,1))
abline(coef=c(0,1))

## ----results='asis',fig.height=5,fig.width=5-----------------------------
dc <- dtangle(Y,pure_samples,n_choose,gamma=.7,markers=marks)
phats <- dc$estimates
plot(truth,phats,xlab="Truth",ylab="Estimates",xlim=c(0,1),ylim=c(0,1))
abline(coef=c(0,1))

## ------------------------------------------------------------------------
get_gamma('microarray-probe')
get_gamma('microarray-gene')
get_gamma('rna-seq')

## ----results='asis',fig.height=5,fig.width=5-----------------------------
dc <- dtangle(Y,pure_samples,n_choose=5,data_type='microarray-gene',markers=marks)
phats <- dc$estimates
plot(truth,phats,xlab="Truth",ylab="Estimates",xlim=c(0,1),ylim=c(0,1))
abline(coef=c(0,1))

## ----results='asis',fig.height=5,fig.width=5-----------------------------
dc <- dtangle(Y,pure_samples,n_choose=c(5,6,7),data_type='microarray-gene',markers=marks)
phats <- dc$estimates
plot(truth,phats,xlab="Truth",ylab="Estimates",xlim=c(0,1),ylim=c(0,1))
abline(coef=c(0,1))

## ----results='asis',fig.height=5,fig.width=5-----------------------------
dc <- dtangle(Y, pure_samples,n_choose,data_type='microarray-gene',marker_method = 'ratio')
phats <- dc$estimate
plot(truth,phats,xlab="Truth",ylab="Estimates",xlim=c(0,1),ylim=c(0,1))
abline(coef=c(0,1))

dc2 <- dtangle(Y, pure_samples,n_choose,data_type='microarray-gene',marker_method = 'diff')
phats2 <- dc2$estimates
points(truth,phats2,col='blue')

dc3 <- dtangle(Y, pure_samples,n_choose,data_type='microarray-gene',marker_method = 'regression')
phats3 <- dc3$estimates
points(truth,phats3,col='red')

