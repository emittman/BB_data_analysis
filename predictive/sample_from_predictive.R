#predictive K-M

sample_glfp <- function(n,p,m1,s1,m2,s2){
  f1 <- rweibull(n, 1/s1, exp(m1))
  f2 <- rweibull(n, 1/s2, exp(m2))
  defective <- sapply(1:n, function(i) sample(c(TRUE,FALSE), 1, prob=c(p[i], 1-p[i])))
  f1 <- ifelse(defective, pmin(f1,f2), f2)
  return(pmin(f1,f2))
}

sample_glfp_rep <- function(n,p,m1,s1,m2,s2,left_tr=0, depth=0){
  if(depth>10) return(rep(0, n))
  f1 <- rweibull(n, 1/s1, exp(m1))
  f2 <- rweibull(n, 1/s2, exp(m2))
  defective <- sample(c(TRUE,FALSE), n, prob=c(p, 1-p), replace=TRUE)
  f1 <- ifelse(defective, pmin(f1,f2), f2)
  f <- pmin(f1,f2)
  d <- depth + 1
  if(!(all(f>left_tr))){
    id <- which(f<=left_tr)
    print(cat("proportion trimmed: ", signif(length(id)/n, 2), "\n"))
    f_sub <- sample_glfp_rep(length(id),p,m1,s1,m2,s2,left_tr[id], d)
    f[id] <- f_sub
  }
  return(f)
}
