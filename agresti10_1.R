rm(list=ls())
library(CatDataAnalysis)
# This is table 10.1 in Agresti 2013
data("table_10.1")
mydat = table_10.1
rm(table_10.1)

# Analysis as in Agresti Table 10.2 
out1 = glm(count ~ r*g + a + c + m, data = mydat, family = poisson()); summary(out1) # Model 1 in Table 9.2
out2 = glm(count ~ r*g + r*(a + c + m) + g*(a + c + m) + a*(c + m) + c*m, data = mydat, family = poisson()); summary(out2) # model 2 in Table 9.2
out3 = glm(count ~  
             r * g * a + 
             r * g * c + 
             r * g * m +
             r * a * c +
             r * a * m +
             r * c * m + 
             g * a * c + 
             g * a * m + 
             g * c * m + 
             a * c * m, data = mydat, family = poisson()); summary(out3) # model 3 in Table 9.2
out4a = glm(count ~ r*g + r*(a + c + m) + g*(a + c + m) + a*(m) + c*m, data = mydat, family = poisson()); summary(out4a) # model 4a in Table 9.2
out4b = glm(count ~ r*g + r*(a + c + m) + g*(a + c + m) + a*(c) + c*m, data = mydat, family = poisson()); summary(out4b) # model 4b in Table 9.2
out4c = glm(count ~ r*g + r*(a + c + m) + g*(a + c + m) + a*(c + m), data = mydat, family = poisson()); summary(out4c)   # model 4c in Table 9.2
out4d = glm(count ~ r*g + r*(a + c + m) + g*(c + m) + a*(c + m) + c*m, data = mydat, family = poisson()); summary(out4d) # ETC...
out4e = glm(count ~ r*g + r*(c + m) + g*(a + c + m) + a*(c + m) + c*m, data = mydat, family = poisson()); summary(out4e) # 
out4f = glm(count ~ r*g + r*(a + c + m) + g*(a + m) + a*(c + m) + c*m, data = mydat, family = poisson()); summary(out4f) # 
out4g = glm(count ~ r*g + r*(a + m) + g*(a + c + m) + a*(c + m) + c*m, data = mydat, family = poisson()); summary(out4g) # 
out4h = glm(count ~ r*g + r*(a + c + m) + g*(a + c) + a*(c + m) + c*m, data = mydat, family = poisson()); summary(out4h) # 
out4i = glm(count ~ r*g + r*(a + c) + g*(a + c + m) + a*(c + m) + c*m, data = mydat, family = poisson()); summary(out4i) # 
out5 = glm(count ~ r*g + r*(a + m) + g*(a + m) + a*(c + m) + c*m, data = mydat, family = poisson()); summary(out5) # 
out6 = glm(count ~ r*g + r*(a) + g*(a + m) + a*(c + m) + c*m, data = mydat, family = poisson()); summary(out6) # 
out7 = glm(count ~ r*g + r*(a) + g*(a) + a*(c + m) + c*m, data = mydat, family = poisson()); summary(out7) # 

# Models including all associations among predictors and all associations among responses
# out8 = glm(count ~ r*g + a*c*m + (r + g) * (a + c + m), data = mydat, family = poisson()); summary(out8) # main effects
# out9 = glm(count ~ r*g + a*c*m + (r + g) * (a * c + a * m + c * m), data = mydat, family = poisson()); summary(out9) # effects of predictors on associations
# 
# out8a = glm(count ~ r*g + a*c*m + (r + g) * (c + m), data = mydat, family = poisson()); summary(out8a) # main effects - 
# out8b = glm(count ~ r*g + a*c*m + (r + g) * (a + m), data = mydat, family = poisson()); summary(out8b) # main effects
# out8c = glm(count ~ r*g + a*c*m + (r + g) * (a + c), data = mydat, family = poisson()); summary(out8c) # main effects

###############################################################################
###############################################################################
# transform data to individual level data
###############################################################################
###############################################################################

mydat2 = as.matrix(mydat[-which(mydat$count == 0), ])
mydatlong = matrix(NA, sum(mydat$count), 5)
tel = 1
for(i in 1:nrow(mydat2)){
  # cat(tel, i, j, cmdat[i,j], "\n")
  mydatlong[tel:(tel-1 + mydat2[i,6]), ] = outer(rep(1, mydat2[i, 6]), mydat2[i, 1:5]) - 1
  tel = tel + mydat2[i,6]
}
rm(i, tel)
colnames(mydatlong) = colnames(mydat[, 1:5])
head(mydatlong)
X = mydatlong[, 4:5]
Y = mydatlong[, 1:3]
XX = cbind(X, X[,1] * X[,2])
rm(mydatlong, mydat2)

###############################################################################
###############################################################################
# analysis with multinomial canonical decomposition
###############################################################################
###############################################################################

source("~/surfdrive/LogitMDA/gsm/mcd.R")

# mutual independence model 
mcd.out1 = mcd(X = XX[, 1:2], Y, S = 0, ord.z = 1, ord.m = 1, trace = FALSE, dcrit = 1e-9); summary(mcd.out1)

# the following model is comparable to the second glm: homogeneous association model
mcd.out2 = mcd(X = XX[, 1:2], Y = Y, S = 2, ord.z = 1, ord.m = 2, trace = FALSE, dcrit = 1e-9); summary(mcd.out2)

# the following model is comparable to the third glm - all three factor terms
mcd.out3 # not fittable with mcd 

mcd.out4a = mcd(X = X[, 1:2], Y = Y, S = 2, ord.z = 1, W = mcd.out2$W[ , -4], trace = FALSE, dcrit = 1e-9); summary(mcd.out4a)
mcd.out4b = mcd(X = X[, 1:2], Y = Y, S = 2, ord.z = 1, W = mcd.out2$W[ , -5], trace = FALSE, dcrit = 1e-9); summary(mcd.out4b)
mcd.out4c = mcd(X = X[, 1:2], Y = Y, S = 2, ord.z = 1, W = mcd.out2$W[ , -6], trace = FALSE, dcrit = 1e-9); summary(mcd.out4c)

mcd.out4a$deviance - mcd.out2$deviance # similar as anova(out2, out4a)
mcd.out4b$deviance - mcd.out2$deviance # similar as anova(out2, out4b)
mcd.out4c$deviance - mcd.out2$deviance # similar as anova(out2, out4c)


mcd.out4d # not fittable with mcd
mcd.out4e # not fittable with mcd
mcd.out4f # not fittable with mcd
mcd.out4g # not fittable with mcd
mcd.out4h # not fittable with mcd
mcd.out4i # not fittable with mcd

mcd.out5 = mcd(X = XX[, 1:2], Y = Y, S = 2, Z = mcd.out2$Z[ , -2], ord.m = 2, trace = FALSE, dcrit = 1e-9); summary(mcd.out5)
mcd.out6# not fittable with mcd
mcd.out7 = mcd(X = XX[, 1:2], Y = Y, S = 1, Z = mcd.out2$Z[ , 1, drop = FALSE], ord.m = 2, trace = FALSE, dcrit = 1e-9); summary(mcd.out7)



