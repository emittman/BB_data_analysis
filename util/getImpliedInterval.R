getImpliedInterval <- function(distr, a=0, b=1, prob=.95, alpha=1-prob){
  callingFun <- paste("q",distr, sep="")
  qtls <- c(alpha/2, 1-alpha/2)
  eval(call(callingFun, qtls, a, b))
}
