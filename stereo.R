sm = function(X, y, S = 2, trace = TRUE, maxiter = 65536, dcrit = 1e-6){
  # stereotype model of Anderson
  # i.e. reduced rank multinomial logistic model
  library(nnet)
  # preparation
  n = nrow(X)
  P = ncol(X)
  G = class.ind(y)
  C = ncol(G)
  ones.C = matrix(1,C,1)
  Jc = diag(C) - 1/C
  ones.n = matrix(1,n,1)
  Jn = diag(n) - 1/n

  # scaling of predictor matrices
  Xoriginal = X
  X = scale(X)
  mx = attr(X, "scaled:center")
  sdx = attr(X, "scaled:scale")

  eig.out = eigen(t(X) %*% X)
  iRx = eig.out$vectors %*% diag(1/sqrt(eig.out$values)) %*% t(eig.out$vectors)
  
  # initialization
  m = colMeans(G)
  udv = svd(1/sqrt(n) * iRx %*% t(X) %*% G)
  B = iRx %*% matrix(udv$u[, 1:S], P, S) * sqrt(n)
  V = matrix(udv$v[, 1:S], C, S) %*% diag(udv$d[1:S], nrow = S, ncol = S) / sqrt(n)
  theta = ones.n %*% t(m) + X %*% B %*% t(V)
  Ghat = exp(theta) / rowSums(exp(theta))
  dev.old = -2 * sum(log(Ghat[G == 1]))

  # iteration
  iter = 0; dif = 1
  while(dif > dcrit){
    iter = iter + 1
    # update m
    Z = theta + 4 * (G - Ghat)
    m = colMeans(Z - X %*% B %*% t(V))

    # update B and V
    theta = ones.n %*% t(m) + X %*% B %*% t(V)
    Ghat = exp(theta) / rowSums(exp(theta))
    Z = theta + 4 * (G - Ghat) - ones.n %*% m
    udv = svd(iRx %*% t(X) %*% Z)
    B = iRx %*% matrix(udv$u[, 1:S], P, S) * sqrt(n)
    V = matrix(udv$v[, 1:S], C, S) %*% diag(udv$d[1:S], nrow = S, ncol = S) / sqrt(n)

    # deviance
    theta = ones.n %*% t(m) + X %*% B %*% t(V)
    Ghat = exp(theta) / rowSums(exp(theta))
    dev.new = -2 * sum(log(Ghat[G == 1]))
    
    # convergence
    dif = 2 * (dev.old - dev.new)/ 2 * (dev.old + dev.new)
    if ( dif < dcrit ) break
    if ( trace ) cat(iter, dev.old, dev.new, dif, "\n")
    if ( iter == maxiter ) warning("Maximum number of iterations reached - not converged (yet)")
    dev.old = dev.new
  } # end iteration

  # create output object
  results = list(
    Xoriginal = Xoriginal,
    X = X,
    mx = mx,
    sdx = sdx,
    G = G,
    m = m,
    B = B,
    U = X %*% B,
    V = V,
    Ghat = Ghat,
    deviance = dev.new,
    iter = iter
  )
  return(results)
}

# esm = function(X, Y, S = 2, Z = NULL, ord = 1, scale.x = FALSE, trace = TRUE, maxiter = 65536, dcrit = 1e-6){
#   # extended stereotype model for multivariate logistic analysis
#   # i.e. double constrained reduced rank multinomial logistic model
#   # preparation
#   cal = match.call()
#   n = nrow(X)
#   ones.n = matrix(1,n,1)
#   P = ncol(X)
#   R = ncol(Y)
#   
#   out = make.profile(Y = Y, ord = ord)
#   if(is.null(Z)){ Z = out$Z }
#   G = out$G
#   A = out$A
#   Q = ncol(Z)  
# 
#   # scaling of predictor matrices
#   Xoriginal = X
#   if(scale.x){
#     X = scale(X)
#     mx = attr(X, "scaled:center")
#     sdx = attr(X, "scaled:scale")
#   }
#   else{
#     X = X
#     mx = rep(0, P)
#     sdx = rep(1, P)
#   }
#   
#   if(P == 1){
#     iRx = 1/sqrt(t(X) %*% X)
#   }
#   else{
#     eig.out = eigen(t(X) %*% X)
#     iRx = eig.out$vectors %*% diag(1/sqrt(eig.out$values)) %*% t(eig.out$vectors)
#   }
#   
#   eig.out = eigen(t(Z) %*% Z)
#   iRz = eig.out$vectors %*% diag(1/sqrt(eig.out$values)) %*% t(eig.out$vectors)
# 
#   # initialization
#   m = colMeans(G)
#   udv = svd(1/sqrt(n) * iRx %*% t(X) %*% G %*% Z %*% iRz)
#   Bx = iRx %*% matrix(udv$u[, 1:S], P, S) * sqrt(n)
#   Bz = iRz %*% matrix(udv$v[, 1:S], Q, S) %*% diag(udv$d[1:S], nrow = S, ncol = S) / sqrt(n)
#   theta = ones.n %*% t(m) + X %*% Bx %*% t(Z %*% Bz)
#   Ghat = exp(theta) / rowSums(exp(theta))
#   dev.old = -2 * sum(log(Ghat[G == 1]))
#   
#   # iteration
#   iter = 0; dif = 1
#   while(dif > dcrit){
#     iter = iter + 1
#     # update m
#     ZZ = theta + 4 * (G - Ghat)
#     m = colMeans(ZZ - X %*% Bx %*% t(Z %*% Bz))
#     
#     # update B and V
#     theta = ones.n %*% t(m) + X %*% Bx %*% t(Z %*% Bz)
#     Ghat = exp(theta) / rowSums(exp(theta))
#     ZZ = theta + 4 * (G - Ghat) - ones.n %*% m
#     udv = svd(iRx %*% t(X) %*% ZZ %*% Z %*% iRz)
#     Bx = iRx %*% matrix(udv$u[, 1:S], P, S) * sqrt(n)
#     Bz = iRz %*% matrix(udv$v[, 1:S], Q, S) %*% diag(udv$d[1:S], nrow = S, ncol = S) / sqrt(n)
#     
#     # deviance
#     theta = ones.n %*% t(m) + X %*% Bx %*% t(Z %*% Bz)
#     Ghat = exp(theta) / rowSums(exp(theta))
#     dev.new = -2 * sum(log(Ghat[G == 1]))
#     
#     # convergence
#     dif = 2 * (dev.old - dev.new)/ ((dev.old + dev.new))
#     if ( trace ) cat(iter, dev.old, dev.new, dif, "\n")
#     if ( dif < dcrit ) break
#     if ( iter == maxiter ) warning("Maximum number of iterations reached - not converged (yet)")
#     dev.old = dev.new
#   } # end iteration
#   
#   BB = Bx %*% t(Bz)
#   
#   # decompose m with full factorial design
#   qrz = qr(model.matrix(formula(paste("~.^", R, sep = "")), data = A))
#   l = solve.qr(qrz, m)
#   l = l[l != 0]
#   mm = qr.fitted(qrz, m)
# 
#   idx = which(colnames(Z) %in% names(l))
#   l.in = l[idx]
#   l.out = l[-idx]
#   
#   # W = model.matrix(formula(paste("~.^", R, sep = "")), data = A)
#   # l = solve(t(W) %*% W, t(W) %*% m)
#   # mm = W %*% l
# 
#   # name giving
#   rownames(Bz) = colnames(Z)
#   rownames(Bx) = colnames(X)
#   rownames(BB) = colnames(X)
#   colnames(BB) = colnames(Z)
#   
#   npar = length(m) - 1 + (P + Q - S) * S
#   # create output object
#   results = list(
#     call = cal,
#     Xoriginal = Xoriginal,
#     X = X,
#     mx = mx,
#     sdx = sdx,
#     Y = Y,
#     pnames = out$Aprofile,
#     xnames = colnames(X),
#     znames = colnames(Z),
#     Z = Z,
#     G = G,
#     m = m,
#     mm = mm,
#     l = l,
#     #l.in = l.in,
#     #l.out = l.out,
#     Bx = Bx,
#     Bz = Bz,
#     A = BB,
#     U = X %*% Bx,
#     V = Z %*% Bz,
#     Ghat = Ghat,
#     deviance = dev.new,
#     df = npar,
#     AIC = dev.new + 2 * npar,
#     iter = iter,
#     svd = udv
#   )
#   class(results) <- "esm"
#   return(results)
# }

# summary.esm = function(object){
#   Q = ncol(object$A)
#   coef.matrix = rbind(object$l[2:(Q+1)] - object$l[1], 
#                       object$A)
#   rownames(coef.matrix)[1] = "1"
#   
#   cat("\n")
#   cat("Call:", "\n")
#   print(object$call) 
#   cat("\n")
#   cat("Fitted extended stereotype model", "\n")
#   cat("Residual deviance:", object$deviance, "\n")
#   cat("Number of fitted parameters", object$df, "\n" )
#   cat("AIC", object$AIC, "\n" )
#   cat("\n")
#   cat("Coefficients:", "\n")
#   print(coef.matrix)
#   cat("\n")
#   if(length(object$l) > (Q+1)){
#     cat("Residual Associations:", "\n")
#     print(object$l[-(1:(Q+1))] - object$l[1])
#   }
# }


esm = function(X, Y, S = 2, Z = NULL, W = NULL, ord.z = 1, ord.m = R, scale.x = FALSE, trace = TRUE, maxiter = 65536, dcrit = 1e-6){
  # extended stereotype model for multivariate logistic analysis
  # i.e. double constrained reduced rank multinomial logistic model
  # preparation
  cal = match.call()
  n = nrow(X)
  ones.n = matrix(1,n,1)
  P = ncol(X)
  R = ncol(Y)
  
  out = make.profile(Y = Y, ord = ord.z)
  if(is.null(Z)){ Z = out$Z }
  G = out$G
  A = out$A
  Q = ncol(Z)  

  if(is.null(W)){
    if(ord.m ==1){
      # W = model.matrix(~., data = A)[, -1]
      qrz = qr(model.matrix(~., data = A)[, -1])
      }
    else{
      # W = model.matrix(formula(paste("~.^", ord.m, sep = "")), data = A)[, -1]
      qrz = qr(model.matrix(formula(paste("~.^", ord.m, sep = "")), data = A)[ , -1])
    }
  }
  # iWWW = solve( t(W) %*% W + 1e-6 * diag(ncol(W)) ) %*% t(W)
  # bm = solve.qr(qrz, m)
  # m = qr.fitted(qrz, m)
  
  # scaling of predictor matrices
  Xoriginal = X
  if(scale.x){
    X = scale(X)
    mx = attr(X, "scaled:center")
    sdx = attr(X, "scaled:scale")
  }
  else{
    X = X
    mx = rep(0, P)
    sdx = rep(1, P)
  }
  
  if(P == 1){
    iRx = 1/sqrt(t(X) %*% X)
  }
  else{
    eig.out = eigen(t(X) %*% X)
    iRx = eig.out$vectors %*% diag(1/sqrt(eig.out$values)) %*% t(eig.out$vectors)
  }
  
  if(ncol(Z) == 1){
    iRz = 1/sqrt(t(Z) %*% Z)
  }
  else{
    eig.out = eigen(t(Z) %*% Z)
    iRz = eig.out$vectors %*% diag(1/sqrt(eig.out$values)) %*% t(eig.out$vectors)
  }
  
  # initialization
  m = colMeans(G)
  # bm =  iWWW %*% m
  # m = W %*% bm
  bm = solve.qr(qrz, m)
  m = qr.fitted(qrz, m)
  #
  udv = svd(1/sqrt(n) * iRx %*% t(X) %*% G %*% Z %*% iRz)
  Bx = iRx %*% matrix(udv$u[, 1:S], P, S) * sqrt(n)
  Bz = iRz %*% matrix(udv$v[, 1:S], Q, S) %*% diag(udv$d[1:S], nrow = S, ncol = S) / sqrt(n)
  theta = ones.n %*% t(m) + X %*% Bx %*% t(Z %*% Bz)
  Ghat = exp(theta) / rowSums(exp(theta))
  dev.old = -2 * sum(log(Ghat[G == 1]))
  
  # iteration
  iter = 0; dif = 1
  while(dif > dcrit){
    iter = iter + 1
    # update m
    ZZ = theta + 4 * (G - Ghat)
    m = colMeans(ZZ - X %*% Bx %*% t(Z %*% Bz))
    # bm = iWWW %*% m
    # m = W %*% bm
    bm = solve.qr(qrz, m)
    m = qr.fitted(qrz, m)

    # update Bx (U) and Bz (V)
    theta = ones.n %*% t(m) + X %*% Bx %*% t(Z %*% Bz)
    Ghat = exp(theta) / rowSums(exp(theta))
    ZZ = theta + 4 * (G - Ghat) - ones.n %*% t(m)
    udv = svd(iRx %*% t(X) %*% ZZ %*% Z %*% iRz)
    Bx = iRx %*% matrix(udv$u[, 1:S], P, S) * sqrt(n)
    Bz = iRz %*% matrix(udv$v[, 1:S], Q, S) %*% diag(udv$d[1:S], nrow = S, ncol = S) / sqrt(n)
    
    # deviance
    theta = ones.n %*% t(m) + X %*% Bx %*% t(Z %*% Bz)
    Ghat = exp(theta) / rowSums(exp(theta))
    dev.new = -2 * sum(log(Ghat[G == 1]))
    
    # convergence
    dif = 2 * (dev.old - dev.new)/ ((dev.old + dev.new))
    if ( trace ) cat(iter, dev.old, dev.new, dif, "\n")
    if ( dif < dcrit ) break
    if ( iter == maxiter ) warning("Maximum number of iterations reached - not converged (yet)")
    dev.old = dev.new
  } # end iteration
  
  BB = Bx %*% t(Bz)
  
  # name giving
  rownames(Bz) = colnames(Z)
  rownames(Bx) = colnames(X)
  rownames(BB) = colnames(X)
  colnames(BB) = colnames(Z)
  rownames(bm) = colnames(W)
  
  npar = length(bm) + (P + Q - S) * S
  # create output object
  results = list(
    call = cal,
    Xoriginal = Xoriginal,
    X = X,
    mx = mx,
    sdx = sdx,
    Y = Y,
    pnames = out$Aprofile,
    xnames = colnames(X),
    znames = colnames(Z),
    Z = Z,
    W = W,
    G = G,
    m = m,
    bm = bm,
    Bx = Bx,
    Bz = Bz,
    A = BB,
    U = X %*% Bx,
    V = Z %*% Bz,
    Ghat = Ghat,
    deviance = dev.new,
    df = npar,
    AIC = dev.new + 2 * npar,
    iter = iter,
    svd = udv
  )
  class(results) <- "esm"
  return(results)
}

boot.esm = function(object){
  # performs a bootstrap for a model fitted with the esm() function
}

summary.esm = function(object){
  Q = ncol(object$A)
  coef.matrix = rbind(object$bm[1:Q],object$A)
  rownames(coef.matrix)[1] = "1"
  
  resids = object$bm[-(1:Q)]
  #names(resids) = rownames(object$bm)[-(1:Q)]

  cat("\n")
  cat("Call:", "\n")
  print(object$call) 
  cat("\n")
  cat("Fitted extended stereotype model", "\n")
  cat("Residual deviance:", object$deviance, "\n")
  cat("Number of fitted parameters", object$df, "\n" )
  cat("AIC", object$AIC, "\n" )
  cat("\n")
  cat("Coefficients:", "\n")
  print(round(coef.matrix, digits = 2))
  cat("\n")
  if(length(object$bm) > (Q)){
    cat("Residual Associations:", "\n")
    print(round(resids, digits = 2))
  }
}

predict.esm = function(object, newX){
  # predict function for extended stereotype model
  m = object$m
  Bx = object$Bx
  Bz = object$Bz
  Z = object$Z
  #
  X = as.matrix(scale(newX, center = object$mx, scale = object$sdx))
  #
  theta = matrix(1, nrow(X), 1) %*% t(m) + X %*% Bx %*% t(Z %*% Bz)
  Ghat = exp(theta) / rowSums(exp(theta))
  #
  theta2 = X %*% Bx %*% t(Bz) # effects on log-odds/log-odds ratio?
  #
  output = list(
    Ghat = Ghat,
    theta2 = theta2
  )
  return(output)
}

plot.esm = function(object, dims = c(1,2), ycol = "darkgreen", xcol = "lightskyblue", ocol = "grey"){
  # plot the results of an extended stereotype model in a 2D biplot
  # object: output object from esm()
  # dims: selction of dimensions for plotting (default first two)
  #
  library(ggplot2)
  ######################################################
  # retrieve information from object
  ######################################################
  Y = object$Y
  
  U = as.data.frame(object$U[ , dims])
  N = nrow(U)
  colnames(U) = c("dim1", "dim2")
  
  V = object$V[ , dims]
  VV = as.data.frame(V)
  K = nrow(V)
  colnames(VV) = c("dim1", "dim2")
  rownames(VV) = object$pnames
  
  ######################################################
  # retrieve information for response variables variable axes
  ######################################################
  
  Z = object$Z
  Q = ncol(Z)
  Bz = object$Bz[ , dims]
  
  # for solid line
  MCz1 <- data.frame(labs=character(),
                     varz = integer(),
                     dim1 = double(),
                     dim2 = double(), stringsAsFactors=FALSE)
  # for markers
  MCz2 <- data.frame(labs=character(),
                     varz = integer(),
                     dim1 = double(),
                     dim2 = double(), stringsAsFactors=FALSE)
  
  ll = 0
  for(q in 1:Q){
    b = matrix(Bz[q , ], 2, 1)
    markers1 = matrix(c(-1,1), 2, 1)
    markerscoord1 = outer(markers1, b) # markers1 %*% t(b %*% solve(t(b) %*% b))
    MCz1[(ll + 1): (ll + 2), 1] = paste0(c("min", "max"), q)
    MCz1[(ll + 1): (ll + 2), 2] = q
    MCz1[(ll + 1): (ll + 2), 3:4] = markerscoord1
    ll = ll + 2
  } # loop q
  
  ######################################################
  # retrieve information for predictor variables variable axes
  ######################################################
  X = object$X
  P = ncol(X)
  Bx = object$Bx[ , dims]
  Xo = object$Xoriginal
  
  # for solid line
  MCx1 <- data.frame(labs=character(),
                     varx = integer(),
                     dim1 = double(),
                     dim2 = double(), stringsAsFactors=FALSE)
  # for markers
  MCx2 <- data.frame(labs=character(),
                     varx = integer(),
                     dim1 = double(),
                     dim2 = double(), stringsAsFactors=FALSE)
  
  ll = 0
  lll = 0
  for(p in 1:P){
    b = matrix(Bx[p , ], 2, 1)
    # solid line
    minx = min(Xo[, p])
    maxx = max(Xo[, p])
    m.x1 = c(minx,maxx)
    markers1 = matrix((m.x1 - object$mx[p])/object$sdx[p], 2, 1)
    markerscoord1 = outer(markers1, b) # markers1 %*% t(b %*% solve(t(b) %*% b))
    MCx1[(ll + 1): (ll + 2), 1] = paste0(c("min", "max"), p)
    MCx1[(ll + 1): (ll + 2), 2] = p
    MCx1[(ll + 1): (ll + 2), 3:4] = markerscoord1
    ll = ll + 2
    # markers
    m.x2 = pretty(Xo[, p])
    m.x2 = m.x2[which(m.x2 > minx & m.x2 < maxx)]
    l.m = length(m.x2)
    markers2 = matrix((m.x2 - object$mx[p])/object$sdx[p], l.m, 1)
    markerscoord2 = outer(markers2, b) # markers2 %*% t(b %*% solve(t(b) %*% b))
    MCx2[(lll + 1): (lll + l.m), 1] = paste(m.x2)
    MCx2[(lll + 1): (lll + l.m), 2] = p
    MCx2[(lll + 1): (lll + l.m), 3:4] = markerscoord2
    lll = lll + l.m
  } # loop p

  ######################################################
  # plotting - objects
  ######################################################
  plt = ggplot() 
  plt = plt + 
    geom_point(data = U, aes_string(x = 'dim1', y = 'dim2'), colour = ocol, alpha = 0.5) +
    xlab(paste("Dimension", dims[1])) +
    ylab(paste("Dimension", dims[2]))

  ######################################################
  # plotting - profile points
  ######################################################
  plt = plt + geom_point(data = VV, aes_string(x = 'dim1', y = 'dim2'), colour = ycol, size = 2) + 
    xlab(paste("Dimension", dims[1])) +
    ylab(paste("Dimension", dims[2])) # + 
    # geom_text(label = rownames(VV)) # werkt nog niet

  
  ######################################################
  # variable axes with ticks and markers for predictors
  ######################################################
  plt = plt + geom_abline(intercept = 0, slope = Bx[,2]/Bx[,1], colour = xcol, linetype = 3) +
      geom_line(data = MCx1, aes_string(x = 'dim1', y = 'dim2', group = 'varx'), col = xcol, linewidth = 1.5) +
      geom_point(data = MCx2, aes_string(x = 'dim1', y = 'dim2'), col = xcol) +
      geom_text(data = MCx2, aes_string(x = 'dim1', y = 'dim2', label = 'labs'), nudge_y = -0.08, size = 1.5)

  ######################################################
  # variable axes for responses
  ######################################################
  plt = plt + geom_abline(intercept = 0, slope = Bz[,2]/Bz[,1], colour = ycol, linetype = 3) + 
    geom_line(data = MCz1, aes_string(x = 'dim1', y = 'dim2', group = 'varz'), col = ycol, linewidth = 1.5) 

  ######################################################
  # variable labels
  ######################################################
  a = max(abs(c(ggplot_build(plt)$layout$panel_scales_x[[1]]$range$range, ggplot_build(plt)$layout$panel_scales_y[[1]]$range$range))) + 0.01
  
  idx2 = apply(abs(Bz), 1, which.max)
  t2 = s2 = rep(NA,Q)
  for(rr in 1:Q){
    t2[rr] = (a * 1.1)/(abs(Bz[rr,idx2[rr]])) * Bz[rr,-idx2[rr]]
    s2[rr] = sign(Bz[rr,idx2[rr]])
  }
  CC2 = cbind(idx2, t2, s2)
  
  idx1 = apply(abs(Bx), 1, which.max)
  t1 = s1 = rep(NA,P)
  for(pp in 1:P){
    t1[pp] = (a * 1.1)/(abs(Bx[pp,idx1[pp]])) * Bx[pp,-idx1[pp]]
    s1[pp] = sign(Bx[pp,idx1[pp]])
    CC1 = cbind(idx1, t1, s1)
  }
  
  CC = rbind(CC1, CC2)
  rownames(CC) = c(object$xnames, object$znames) 
  colnames(CC) = c("idx", "t", "s")
  
  bottom = which(CC[, "idx"] == 2 & CC[, "s"] == -1)
  top =  which(CC[, "idx"] == 2 & CC[, "s"] == 1)
  right = which(CC[, "idx"] == 1 & CC[, "s"] == 1)
  left = which(CC[, "idx"] == 1 & CC[, "s"] == -1)
  
  
  if(length(CC[top, "t"])==0){
    plt <- plt + scale_x_continuous(limits = c(-a,a), breaks = CC[bottom, "t"], labels = rownames(CC)[bottom],
                                    sec.axis = sec_axis(trans ~ ., breaks = 0, labels = ""))
  }else{
    plt <- plt + scale_x_continuous(limits = c(-a,a), breaks = CC[bottom, "t"], labels = rownames(CC)[bottom],
                                    sec.axis = sec_axis(trans ~ ., breaks = CC[top, "t"], labels = rownames(CC)[top]))
  }
  
  
  if(length(CC[right, "t"])==0){
    plt = plt + scale_y_continuous(limits = c(-a,a), breaks = CC[left, "t"], labels = rownames(CC)[left],
                                   sec.axis = sec_axis(trans ~ ., breaks = 0, labels = ""))
  }else{
    plt = plt + scale_y_continuous(limits = c(-a,a), breaks = CC[left, "t"], labels = rownames(CC)[left],
                                   sec.axis = sec_axis(trans ~ ., breaks = CC[right, "t"], labels = rownames(CC)[right]))
  }
  
  
  plt = plt + theme(axis.line = element_line(colour = "black"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank())
  plt = plt + coord_fixed()
  
  #suppressWarnings(print(plt))
  
  return(plt)
}

make.profile = function(Y, ord){
  library(nnet)
  n = nrow(Y)
  profile = matrix(NA, n, 1)
  for(i in 1:n){profile[i, 1] = paste(Y[i,], collapse = "")}
  G = class.ind(profile)
  C = ncol(G)
  
  A = unique(Y)
  Aprofile = matrix(NA, nrow(A),1)
  AA = ifelse(A == 0, "-", "+")
  for(i in 1:nrow(AA)){Aprofile[i, 1] = paste(AA[i,], sep="", collapse = "")}
  Ai = t(class.ind(Aprofile))
  aidx = rep(NA, nrow(A)); for(j in 1:nrow(A)){aidx[j] = which(Ai[j,] == 1, arr.ind = TRUE)}
  A = A[aidx,]
  A = (2 * A - 1)/2
  A = as.data.frame(A)
  rownames(A) = colnames(G)
  colnames(A) = colnames(Y)
  if (ord == 1){
    Z = model.matrix(~ ., data = A)[, -1]
  }
  else if (ord == 2){
    Z = model.matrix(~ .^2, data = A)[, -1]
  }
  else if (ord == 3){
    Z = model.matrix(~ .^3, data = A)[, -1]
  }
  else if(ord == R){
    Z = diag(nrow(A))
  }
  #
  output = list(
    Z = Z,
    G = G, 
    A = A,
    Aprofile = Aprofile)
  return(output)
}
