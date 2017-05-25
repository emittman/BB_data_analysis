#predictive K-M

sample_glfp <- function(n,p,m1,s1,m2,s2){
  f1 <- rweibull(n, 1/s1, exp(m1))
  f2 <- rweibull(n, 1/s2, exp(m2))
  defective <- sapply(1:n, function(i) sample(c(TRUE,FALSE), 1, prob=c(p[i], 1-p[i])))
  f1 <- ifelse(defective, pmin(f1,f2), f2)
  return(pmin(f1,f2))
}