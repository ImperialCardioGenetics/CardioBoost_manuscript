#permutation test: Test whether the stat fun computed on two datasets are significant different 
## About permutation test: compute the sampling distribution of any test statistic under null hypothesis
## If the null hypothesis is true, then changing the exposure should not have effect on the outcome.
## Null hypothesis: the two datasets are same. -> we can combine them to creat one whole dataset (i.e., permute the classifier labels)
## creat the distribution of the test stat under null hypothesis: randomly choose two datasets of equal size and then calculate the stat
## and estimate the prob that observed stat is larger than the proposed one. Return that prob as p-value
## Alternative hypothesis: the two datasets are different
# Reference: 
# [1]http://faculty.washington.edu/kenrice/sisg/SISG-08-06.pdf
# [2]http://genomicsclass.github.io/book/pages/permutation_tests.html


permutation<-function(data1,data2,fun,n=2000){
  d.obs<-abs(fun(data1)-fun(data2))
  tail.prob<-0
  for(i in 1:n){
  data<-sample_n(rbind(data1,data2),size=nrow(data1)+nrow(data2))
  v1<-data[1:nrow(data1),]
  v2<-data[nrow(data1)+1:nrow(data2),]
  d.sim<-abs(fun(v1)-fun(v2))
  if(d.sim>=d.obs){
    tail.prob<-tail.prob+1
  }
  }
  tail.prob <- round((tail.prob+1)/(n+1),7)
  
  return(tail.prob)
  
}

