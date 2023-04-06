cmdat = VGAM::coalminers
cmdat
cmdat2 = matrix(NA, sum(cmdat[, 1:4]), 3)
tel = 1
for(i in 1:9){
  for(j in 1:4){
    # cat(tel, i, j, cmdat[i,j], "\n")
    cmdat2[tel:(tel-1 + cmdat[i,j]), 1] = cmdat[i, 5]
    cmdat2[tel:(tel-1 + cmdat[i,j]), 2] = ifelse((j == 1)| (j == 2), 1, 0)
    cmdat2[tel:(tel-1 + cmdat[i,j]), 3] = ifelse((j == 1)| (j == 3), 1, 0)
    tel = tel + cmdat[i,j]
  }
}
colnames(cmdat2) = c("age", "B", "W")
X = cmdat2[, 1, drop = FALSE]
X = (X - 42)/5
XX = cbind(X, X^2)
Y = cmdat2[, 2:3]

# Extended Stereotype Model - analysis
source("~/surfdrive/LogitMDA/gsm/stereo.R")

outa = esm(X = XX[ , 1, drop = FALSE], Y = Y, S = 1, ord.z = 2) # only age
summary(outa)

# VGAM analysis
library(VGAM)
coalminers <- transform(coalminers, Age = (age - 42) / 5)
fit <- vglm(cbind(nBnW, nBW, BnW, BW) ~ Age, binom2.or(zero = NULL), data = coalminers)
summary(fit)
