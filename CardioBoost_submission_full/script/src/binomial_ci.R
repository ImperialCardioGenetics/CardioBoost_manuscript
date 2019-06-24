#binomial Interval to estimate the CI of PR AUC
#http://pages.cs.wisc.edu/~boyd/aucpr_final.pdf

#prauc -  the point estimate of prauc
#n -  the number of positive class
#alpha corresponding to Confidence Interval: [alpha/2, 1-alpha/2]
binomial_ci<-function(prauc,n,alpha){
  z=qnorm(1-alpha/2)
  se=z*sqrt(prauc*(1-prauc)/n)
  lower=prauc-se
  upper=prauc+se
  c(lower,upper)
}
