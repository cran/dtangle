## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----cache=TRUE----------------------------------------------------------
library('dtangle')
names(shen_orr_ex)

## ----cache=TRUE----------------------------------------------------------
head(shen_orr_ex$annotation$mixture)

## ----cache=TRUE----------------------------------------------------------
Y <- shen_orr_ex$data$log
Y[1:4,1:4]

## ----cache=TRUE----------------------------------------------------------
library('dtangle')
data = shen_orr_ex$data$log
mixture_proportions = shen_orr_ex$annotation$mixture

## ----cache=TRUE----------------------------------------------------------
mixture_proportions

## ----cache=TRUE----------------------------------------------------------
pure_samples = list(Liver=c(1,2,3),Brain=c(4,5,6),Lung=c(7,8,9))

dt_out = dtangle(Y=data, pure_samples = pure_samples)

matplot(mixture_proportions,dt_out$estimates, xlim = c(0,1),ylim=c(0,1),xlab="Truth",ylab="Estimates")

## ----cache=TRUE----------------------------------------------------------
mixture_samples = data[-(1:9),]
reference_samples = data[1:9,]

dt_out = dtangle(Y=mixture_samples, reference=reference_samples,pure_samples = pure_samples)

mixture_mixture_proportions = mixture_proportions[-(1:9),]
matplot(mixture_mixture_proportions,dt_out$estimates, xlim = c(0,1),ylim=c(0,1),xlab="Truth",ylab="Estimates")

## ----cache=TRUE----------------------------------------------------------
ref_reduced = t(sapply(pure_samples,function(x)colMeans(reference_samples[x,,drop=FALSE])))

dt_out = dtangle(Y=mixture_samples, reference=ref_reduced)

matplot(mixture_mixture_proportions,dt_out$estimates, xlim = c(0,1),ylim=c(0,1),xlab="Truth",ylab="Estimates")

## ----cache=TRUE----------------------------------------------------------
dt_out = dtangle(Y=mixture_samples, references = ref_reduced)

## ----cache=TRUE----------------------------------------------------------
dt_out = dtangle(Y=mixture_samples, references = ref_reduced,marker_method = "diff")

## ----cache=TRUE----------------------------------------------------------
dt_out$n_markers

## ----cache=TRUE----------------------------------------------------------
dt_out = dtangle(Y=mixture_samples, references = ref_reduced,marker_method = "diff",n_markers=100)

dt_out$n_markers

## ----cache=TRUE----------------------------------------------------------
dt_out = dtangle(Y=mixture_samples, references = ref_reduced,marker_method = "diff",n_markers=c(100,150,50))

dt_out$n_markers

## ----cache=TRUE----------------------------------------------------------
dt_out = dtangle(Y=mixture_samples, references = ref_reduced,marker_method = "diff",n_markers=.075)

dt_out$n_markers

dt_out = dtangle(Y=mixture_samples, references = ref_reduced,marker_method = "diff",n_markers=c(.1,.15,.05))

dt_out$n_markers

## ----cache=TRUE----------------------------------------------------------
marker_genes = list(c(120,253,316),
                    c(180,429,14),
                    c(1,109,206))

dt_out = dtangle(Y=mixture_samples, references = ref_reduced,markers=marker_genes)
dt_out$n_markers

## ----cache=TRUE----------------------------------------------------------
mrkrs = find_markers(Y=mixture_samples, references = ref_reduced)
names(mrkrs)

## ----cache=TRUE----------------------------------------------------------
dt_out = dtangle(Y = mixture_samples,references = ref_reduced,markers=mrkrs,n_markers=.1)

## ----cache=TRUE----------------------------------------------------------
dt_out = dtangle(Y = mixture_samples,references = ref_reduced,gamma=.9)

## ----cache=TRUE----------------------------------------------------------
dt_out = dtangle(Y = mixture_samples,references = ref_reduced,data_type="microarray-gene")

## ----cache=TRUE----------------------------------------------------------
dtangle:::gma

## ----cache=TRUE----------------------------------------------------------
dt_out = dtangle(Y = mixture_samples,references = ref_reduced,summary_fn=median)
head(dt_out$estimates)
dt_out = dtangle(Y = mixture_samples,references = ref_reduced,summary_fn=mean)
head(dt_out$estimates)

