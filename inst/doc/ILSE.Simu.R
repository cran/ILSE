## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library("ILSE")

## ----eval = FALSE-------------------------------------------------------------
#    n <- 100
#    p <- 6
#    X <- MASS::mvrnorm(n, rep(0, p), cor.mat(p, rho=0.5))
#    beta0 <- rep(c(1,-1), times=3)
#    Y <- -2+ X %*% beta0 + rnorm(n, sd=1)

## ----eval = FALSE-------------------------------------------------------------
#  ilse1 <- ilse(Y~X)
#  print(ilse1)

## ----eval = FALSE-------------------------------------------------------------
#  dat <- data.frame(Y=Y, X=X)
#  ilse1 <- ilse(Y~., data=dat)
#  print(ilse1)
#  Coef(ilse1) # access the coefficients
#  Fitted.values(ilse1)
#  Residuals(ilse1)
#  

## ----eval = FALSE-------------------------------------------------------------
#  s1 <- summary(ilse1)
#  s1

## ----eval = FALSE-------------------------------------------------------------
#  mis_rate <- 0.3
#  set.seed(1)
#  na_id <- sample(1:(n*p), n*p*mis_rate)
#  Xmis <- X
#  Xmis[na_id] <- NA
#  ncomp <- sum(complete.cases(Xmis))
#  message("Number of complete cases is ", ncomp, '\n')

## ----eval = FALSE-------------------------------------------------------------
#  lm1 <- lm(Y~Xmis)
#  s_cc <- summary.lm(lm1)
#  s_cc

## ----eval = FALSE-------------------------------------------------------------
#  ilse2 <- ilse(Y~Xmis, data=NULL, infor_output=T)
#  print(ilse2)

## ----eval = FALSE-------------------------------------------------------------
#  s2 <- summary(ilse2, Nbt=20)
#  s2

## ----eval = FALSE-------------------------------------------------------------
#  fimllm <- fimlreg(Y~Xmis)
#  print(fimllm)
#  

## ----eval = FALSE-------------------------------------------------------------
#  s_fiml <- summary(fimllm, Nbt=20)
#  s_fiml

## ----eval = FALSE-------------------------------------------------------------
#  pMat <- cbind(CC=s_cc$coefficients[,4], ILSE=s2[,4], FIML=s_fiml[,4])
#  library(ggplot2)
#  df1 <- data.frame(Pval= as.vector(pMat[-1,]),
#                      Method =factor(rep(c('CC', "ILSE", "FIML"),each=p)),
#                      covariate= factor(rep(paste0("X", 1:p), times=3)))
#  ggplot(data=df1, aes(x=covariate, y=Pval, fill=Method)) + geom_bar(position = "dodge", stat="identity",width = 0.5) + geom_hline(yintercept = 0.05, color='red') + geom_hline(yintercept = 0.1, color='blue')

## -----------------------------------------------------------------------------
sessionInfo()

