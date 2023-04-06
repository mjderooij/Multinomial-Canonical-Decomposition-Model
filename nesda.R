rm(list=ls())
load("/Users/rooijmjde/surfdrive/BinomialMapping/nesda1.Rdata")
rm(xval.out, out)
X[ , "Gender"] = ifelse(X[ , "Gender"] > 0, 1, 0)
source("~/surfdrive/LogitMDA/gsm/stereo.R")

# AICs needed - what is the effective dimension?
out5 = esm(X, Y, S = 5, ord.z = 3)
out4 = esm(X, Y, S = 4, ord.z = 3)
out3 = esm(X, Y, S = 3, ord.z = 3)
out2 = esm(X, Y, S = 2, ord.z = 3)
out1 = esm(X, Y, S = 1, ord.z = 3)

#c(out5$deviance, out4$deviance, out3$deviance, out2$deviance, out1$deviance)
round(c(out5$AIC, out4$AIC, out3$AIC, out2$AIC, out1$AIC), digits = 1)

out2a = esm(X, Y, S = 2, ord.z = 3)
out2b = esm(X, Y, S = 2, ord.z = 2)
out2c = esm(X, Y, S = 2, ord.z = 1)

#c(out2a$deviance, out2b$deviance, out2c$deviance)
round(c(out2a$AIC, out2b$AIC, out2c$AIC), digits = 1)

# in between ord = 1 and ord = 2 models
Z = out2b$Z
out = list()
AICs = rep(NA, 10)
for(c in 6:15){
  i = c - 5
  out[[i]] = esm(X, Y, S = 2, Z = Z[ , c(1,2,3,4,5,c)])
  AICs[i] = out[[i]]$AIC
}
round(AICs, digits = 1)
sort(AICs)
out2c$AIC
which(AICs < out2c$AIC) + 5
Z = Z[ , c(1,2,3,4,5,7)]

# out2c.1 = esm(X[ , -1], Y, S = 2, Z = Z)
# out2c.2 = esm(X[ , -2], Y, S = 2, Z = Z)
# out2c.3 = esm(X[ , -3], Y, S = 2, Z = Z)
out2c.4 = esm(X[ , -4], Y, S = 2, Z = Z)
out2c.5 = esm(X[ , -5], Y, S = 2, Z = Z)
out2c.6 = esm(X[ , -6], Y, S = 2, Z = Z)
out2c.7 = esm(X[ , -7], Y, S = 2, Z = Z)
out2c.8 = esm(X[ , -8], Y, S = 2, Z = Z)

# c(out2c$deviance, 
#   out2c.1$deviance, 
#   out2c.2$deviance, 
#   out2c.3$deviance, 
#   out2c.4$deviance, 
#   out2c.5$deviance, 
#   out2c.6$deviance, 
#   out2c.7$deviance, 
#   out2c.8$deviance) - out2c$deviance
round(c(out2c.4$AIC, out2c.5$AIC, out2c.6$AIC, out2c.7$AIC, out2c.8$AIC), digits = 1)
out2c.67 = esm(X[ , -c(6, 7)], Y, S = 2, Z = Z)
out2c.67$AIC

# residual dependencies
out2c67.4 = esm(X[ , -c(6, 7)], Y, S = 2, Z = Z, ord.m = 4)
out2c67.3 = esm(X[ , -c(6, 7)], Y, S = 2, Z = Z, ord.m = 3)
out2c67.2 = esm(X[ , -c(6, 7)], Y, S = 2, Z = Z, ord.m = 2)

round(c(out2c67.4$AIC, out2c67.3$AIC, out2c67.2$AIC), digits = 1)

outfinal = out2c67.2
summary(outfinal)



####################################################################################
####################################################################################
# checking downward compatibility
####################################################################################
####################################################################################

# D and G involving the association
# out = esm(X, Y[ , c(1,3)], S = 2, ord = 2)
# out$A
# outfinal$A[, c(1, 3, 6)]
# abs(out$A - outfinal$A[, c(1, 5, 6)])
# 
# # M and S involving no association
# out = esm(X, Y[ , c(2,4)], S = 2, ord = 1)
# out$A
# outfinal$A[, c(2, 4)]
# abs(out$A - outfinal$A[, c(2, 4)])
# 
# # M and G involving partial association
# out = esm(X, Y[ , c(2,3)], S = 2, ord = 1)
# out$A
# outfinal$A[, c(2, 3)]
# abs(out$A - outfinal$A[, c(2, 3)])
