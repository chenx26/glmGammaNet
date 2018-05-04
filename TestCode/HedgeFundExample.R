rm(list=ls())
library(EstimatorStandardError)
library(glmGammaNet)
#--------- Load data
data(edhec)
colnames(edhec)=c("CA", "CTAG", "DIS", "EM",
                  "EMN", "ED", "FIA", "GM", "L/S",
                  "MA", "RV", "SS", "FoF")
nboot=500
return.coeffs = TRUE
d.GLM.EN = 12
alpha.ES = 0.05
alpha.EN = 0.5
#--------- Set random seed
seed = 1234
k_fold_iter = 1000
set.seed(seed)
#--------- Expected Shortfall estimates and their S.E.'s
se.ES=ES.SE(edhec, p = 1 - alpha.ES, method="historical",nsim = nboot,
            se.method = c("IFiid","IFcor"),
#            se.method = c("IFiid","IFcor","BOOTiid","BOOTcor"),
            standardize = FALSE, return.coeffs = return.coeffs, 
            k_fold_iter = k_fold_iter,
            d.GLM.EN = d.GLM.EN,
            alpha = alpha.EN)


tmp = se.ES$IFcor
SEs = sapply(tmp, function(x) x[[1]])
se.ES$IFcor = SEs

se.ES.df = printSE(se.ES, round.digit  = 3, valonly = TRUE)

ar1 = apply(edhec, 2,
            function(x) {
                    res = ar(x, order.max = 1)
                    if(length(res$ar) == 0)
                            return(0)
                    return(res$ar[1])
            }
)

IF.ar1 = apply(edhec, 2,
               function(x) {
                       res = ar(ES.IF(x), order.max = 1)
                       if(length(res$ar) == 0)
                               return(0)
                       return(res$ar[1])
               }
) 

PI.IFcor = (se.ES.df$IFcor / se.ES.df$IFiid - 1) * 100

# use newey-west to compute SE's

NW = NULL
for(i in 1:ncol(edhec)){
        data = coredata(edhec[,i])
        data.IF = ES.IF(data, alpha.ES = alpha.ES)
        NW = c(NW, sqrt(nse.nw(data.IF)))
}
PI.NW = (NW / se.ES.df$IFiid - 1) * 100

# use SE.glmGammaNet to compute the SE's

i=1
seCorIFGamma = list()
for(i in 1:ncol(edhec)){
        data = coredata(edhec[,i])
        data.IF = ES.IF(data, alpha.ES = alpha.ES)
        seCorIFGamma[[i]] = SE.glmGammaNet(data.IF, d = d.GLM.EN, alpha.EN = alpha.EN,
                                           nfolds = 100,
                                            # abs_tol = abs_tol,
                                            # standardize = TRUE,
                                            #               has_intercept = has_intercept,
                                            #               num_lambda = num_lambda,
                                            #               twosided = FALSE,
                                            
                                            #               keep = keep
                                            return.coeffs = TRUE
        )
}
SE.seCorIFGamma = sqrt(sapply(seCorIFGamma, function(x) x[[1]]))
PI.seCorIFGamma = (SE.seCorIFGamma / se.ES.df$IFiid - 1) * 100

se.ES.df

se.ES.df = cbind(se.ES.df, PI.IFcor, SE.seCorIFGamma, PI.seCorIFGamma)

# se.ES.df = cbind(se.ES.df, PI.IFcor)[c(1,2,3,6,4,5)]
# SE.glmGammaNet(rnorm(10), d = d.GLM.EN, alpha.EN = alpha.EN,
#                # abs_tol = abs_tol,
#                # standardize = TRUE,
# #               has_intercept = has_intercept,
# #               num_lambda = num_lambda,
# #               twosided = FALSE,
#                
# #               keep = keep
#                 return.coeffs = TRUE
# )

se.ES.df = cbind(se.ES.df, ar1, IF.ar1)

# SE.glmGammaNet(as.numeric(edhec[,1]), d = d.GLM.EN, alpha.EN = 0.5,
#                # abs_tol = abs_tol,
#                # standardize = TRUE,
#                #               has_intercept = has_intercept,
#                #               num_lambda = num_lambda,
#                #               twosided = FALSE,
#                
#                #               keep = keep
#                return.coeffs = TRUE
# )
# 
# res.gamma = apply(edhec, 2, function(x) SE.glmGammaNet(data, d = d.GLM.EN, alpha.EN = 0.5,
#                                                        # abs_tol = abs_tol,
#                                                        # standardize = TRUE,
#                                                        #               has_intercept = has_intercept,
#                                                        #               num_lambda = num_lambda,
#                                                        #               twosided = FALSE,
#                                                        
#                                                        #               keep = keep
#                                                        return.coeffs = TRUE
# )
# )

colnames(se.ES.df) = c("ES", "seIidIF", "seCorIF", "PI.seCorIF", "seCorIF.gamma", "PI.seCorIF.gamma", "ar1", "IF.ar1")
knitr::kable(se.ES.df, digits = c(3,3,3,0,3,0,2,2), align = 'c')


N = nrow(edhec)
abstol = 1e-3
coeffs = sapply(tmp, function(x) x[[2]])
coeffs = t(coeffs)
# cbind(se.ES.df[,3], sqrt(exp(coeffs[,1])/N))
coeffs = cbind(coeffs, apply(coeffs, 1, function(x) sum(abs(x) <= abstol)))
rownames(coeffs) = colnames(edhec)
colnames(coeffs) = c(paste0("beta", 0:(d.GLM.EN)), "zero.coeffs")
if(d.GLM.EN <= 7)
        knitr::kable(coeffs, digits = c(rep(3, (d.GLM.EN + 1)), 0))
if(d.GLM.EN > 7)
        knitr::kable(coeffs[,c(1:4, (ncol(coeffs) - 4):ncol(coeffs))], 
                     digits = c(rep(3, 8), 0),
                     align = 'c')

N = nrow(edhec)
abstol = 1e-3
coeffs = sapply(seCorIFGamma, function(x) x[[2]])
coeffs = t(coeffs)
# cbind(se.ES.df[,3], sqrt(exp(coeffs[,1])/N))
coeffs = cbind(coeffs, apply(coeffs, 1, function(x) sum(abs(x) <= abstol)))
rownames(coeffs) = colnames(edhec)
colnames(coeffs) = c(paste0("beta", 0:(d.GLM.EN)), "zero.coeffs")
if(d.GLM.EN <= 7)
        knitr::kable(coeffs, digits = c(rep(3, (d.GLM.EN + 1)), 0))
if(d.GLM.EN > 7)
        knitr::kable(coeffs[,c(1:4, (ncol(coeffs) - 4):ncol(coeffs))], 
                     digits = c(rep(3, 8), 0),
                     align = 'c')
