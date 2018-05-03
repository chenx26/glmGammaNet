SE.glmGammaNet = function (data, ..., d = 7, alpha.EN = 0.5, keep = 1, standardize = FALSE,
                           return.coeffs = FALSE, prewhiten = FALSE)
{
        if (prewhiten) {
                res.ar = ar(data, ...)
                ar.coeffs = res.ar$ar
                data = na.omit(res.ar$resid)
        }
        N = length(data)
        my.periodogram = myperiodogram(data, ..., keep = keep)
        my.freq = my.periodogram$freq
        my.periodogram = my.periodogram$spec
        nfreq = length(my.freq)
        x.mat = rep(1, length(my.freq))
        for (col.iter in 1:d) {
                x.mat = cbind(x.mat, my.freq^col.iter)
        }
        mean_vec = apply(x.mat[, -1], 2, mean)
        sd_vec = apply(x.mat[, -1], 2, sd)
        if (standardize) {
                for (i in 2:ncol(x.mat)) {
                        tmp = x.mat[, i]
                        x.mat[, i] = (tmp - mean(tmp))/sd(tmp)
                }
        }
        res = cv.glmGammaNet(x.mat, my.periodogram, alpha.EN = alpha.EN)
        res = res$x
        if (return.coeffs) {
                if (standardize) {
                        variance = exp(sum(res * c(1, -mean_vec/sd_vec)))/N
                        if (prewhiten)
                                variance = variance/(1 - sum(ar.coeffs))^2
                        coeffs = res
                        return(list(variance, coeffs))
                }
                variance = exp(res[1])/N
                if (prewhiten)
                        variance = variance/(1 - sum(ar.coeffs))^2
                coeffs = res
                return(list(variance, coeffs))
        }
        if (standardize) {
                variance = exp(sum(res * c(1, -mean_vec/sd_vec)))/N
                if (prewhiten)
                        variance = variance/(1 - sum(ar.coeffs))^2
                return(variance)
        }
        variance = exp(res[1])/N
        if (prewhiten)
                variance = variance/(1 - sum(ar.coeffs))^2
        return(variance)
}

