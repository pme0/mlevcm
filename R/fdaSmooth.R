#' @title  ???
#'
#' @description  ???.
#'
#' @param  ???
#'
#' @export
#'
fdaSmooth <- function(X_rough, wgts=1, wvls=350:2500){

  # Bspline_dim <- floor(0.5*ncol(X_rough))
  # t_range <- 350:2500
  # NN <- nrow(X_rough);   PP <- ncol(X_rough);
  # breaks <- seq(t_range[1], t_range[length(t_range)], len = Bspline_dim-2)
  # B    <- bsplineS(t_range, breaks, norder=4, returnMatrix=T)
  # B_d2 <- bsplineS(t_range, breaks, norder=4, returnMatrix=T, nderiv=2) # 2nd derivative
  # Pen <- t(B_d2) %*% B_d2  # penmatalty matrix
  # ll <- 10;   lams <- 10^seq(-7, -2, len=ll);   gcv <- matrix(0,NN,ll)
  # W <- rep(1,P) #1/apply(X, 2, var)
  # BB <- diag(sqrt(W)) %*% B
  # BtB <- t(BB) %*% BB
  # tB <- t(BB)

  NN <- nrow(X_rough)
  PP <- ncol(X_rough)
  ww <- if(length(wgts)==1){ rep(wgts,PP) }else{ wgts }
  X_smooth <- matrix(0,NN,PP);   lam_opt <- rep(NA,NN)

  for(i in 1:NN){
    XX <- X_rough[i,] #sqrt(W) * X_rough[i,]
    # find optimal lambda
    # for(l in 1:ll){
    #   H <- BB %*% solve(BtB + PP * lams[l] * Pen) %*% tB
    #   I_minus_H <- diag(PP) - H
    #   num <- (1/PP) * norm(I_minus_H %*% XX, "F")^2
    #   den <- ((1/PP) * sum(diag(I_minus_H)))^2
    #   gcv[i,l] <- num / den
    # }
    # plot(lams, gcv[i,], type="l", lwd=2); abline(v=lam_opt, col="red"); legend("topright", paste0("lam=",lam_opt), bty="n", cex=0.8)
    #lam_opt[i] <- lams[which.min(gcv[i,])]
    #X_smooth[i,] <- drop(BB %*% solve(BtB + PP * lam_opt[i] * Pen) %*% tB %*% X_rough[i,]) # plot(X_rough[i,],t="l"); lines(X_smooth[i,],col=2)
    aux <- smooth.spline(x = wvls, w = ww, y = X_rough[i,], cv=F)
    X_smooth[i,] <- predict(aux)$y
    lam_opt[i] <- aux$lambda
    #print(paste0(i, " out of ", NN))
  }

  return(X_smooth)
}

# X_smooth <- fdaSmooth(X_rough, wvls)
#
# i <- 2
# X_mult <- X_rough[i,] %*% t(t(B) / colSums(B))
# plot(t_range, X_rough[i,], type="l")
# lines(breaks, X_mult[-c(1,ncol(X_mult))], col="green")
# lines(t_range, X_smooth[i,], col="red")

#write.table(X_smooth, file="nirs/data/datasets/experiments_2017_SouthKensington/A_/X_smooth__A_.txt", sep=",", col.names=F, row.names=F)


