## ------------------------------------------------------------------------
patientID <- c(1,2,3,4)
age <- c(25,34,28,52)
diabetes <- c("Type1","Type2","Type1","Type1")
status <- c("Poor","Improved","Excellent","Poor")
patientdata <- data.frame(patientID,age,diabetes,status)
patientdata
table(patientdata$diabetes,patientdata$status)

## ------------------------------------------------------------------------
mpg <- c(39.4,36.4,33.8,31.5,29.6,27.8,26.3,24.9)
wt <- c(1897,1665,1559,1550,1490,1500,1480,1380)
lm(mpg~wt)
lmfit <- lm(mpg~wt)
summary(lmfit)
plot(lmfit)

## ------------------------------------------------------------------------
x <- c(0,1,2,3,4)
p <- c(.1,.2,.2,.2,.3)
cp <- cumsum(p)
m <- 1000
r <- numeric(m);
#find all indexes of x we want under the condition we design: F(x[i-1])<u<=F(x[i])??then x=x[i]
r <- x[findInterval(runif(m),cp)+1]
ct <- as.vector(table(r))
ct
#caculate the probabilty in the random sample we generate
tp <- ct/sum(ct)
ratio <- tp/p
data.frame(x,p,tp,ratio)
# #the relative frequency table comparing the empirical with the theoretical probabilities.

## ------------------------------------------------------------------------
set.seed(1234)
r2 <- sample(c(0,1,2,3,4),m,replace=TRUE,prob=c(.1,.2,.2,.2,.3))
ct2 <- as.vector(table(r2))
tp2 <- ct2/sum(ct2)
ratio <- tp2/p
data.frame(x,p,tp2,ratio)
#R sample function

## ------------------------------------------------------------------------
n <- 1000
j<-k<-0
y <- numeric(n)
while (k < n) 
  {
    u <- runif(1)
    j <- j + 1
    x <- runif(1) #random variate from g
    if (x^2 * (1-x) > u) 
      {
        #we accept x
        k <- k + 1
        y[k] <- x
      }
   }
j
# #(experiments) for n random numbers

## ----echo=TRUE-----------------------------------------------------------
hist(y, prob = TRUE, main = expression(f(x)==12*x^2*(1-x)))
x <- seq(0, 1, .001)
lines(x, 12*x^2*(1-x),col="red")
# #Histogram of the sample with the theoretical $Beta(3,2)$ density superimposed,which shows that the random sample generated fits the theoretical model $Beta(3,2)$ well. 

## ----fig.height=3--------------------------------------------------------
library(evir)
n <- 1000; r <- 4; b <- 2
lambda <- rgamma(n, r, b) #lambda is random
#generate the mixture  
x <- rexp(n, lambda)
#compare the mixture with the GPD
hist(x,prob=TRUE,col="grey",main="Density of GPD(0.25,0,0.5)")
y <- seq(0, 100, .001)
lines(y,1/2^4*(1/2+1/4*y)^(-5),col="red") #the density function of generalized Pareto distribution:1/2^4*(1/2+1/4*x)^(-5)

# #The figure above shows that the random model generated is good especially when x extends to right further.
  

## ------------------------------------------------------------------------
#the theoretical model: generalized Pareto distribution
y1 <- rgpd(n,xi=1/r,mu=0,beta=b/r)
par(mfrow=c(1,2))
hist(x,main="Histogram of Simulation",col="pink")
hist(y1,main="Histogram of GPD",col="light green")
# #The comparing histograms shows that the simulation is good.

## ------------------------------------------------------------------------
MCofBeta <- function(n,t,a,b){
  x <- runif(n, 0, t)    #generate random x
  theta.hat <- mean(x^(a-1)*(1-x)^(b-1))*t
}
print(c(MCofBeta(1e4,1,3,3),beta(3,3)))   #compare the MC estimation and the theoretical value
##the result shows that the MC estimate works well

## ------------------------------------------------------------------------
x <- y <- ratio<-c(.1,.2,.3,.4,.5,.6,.7,.8,.9)
for (i in 1:9){
  F[i] <- MCofBeta(1e4,x[i],3,3)/MCofBeta(1e4,1,3,3)  #generate F(x) based on estimation B(x;a,b)
  y[i]<-pbeta(x[i],3,3)      #theretical value
  ratio[i] <- F[i]/y[i]
}
  data.frame(F,y,ratio)     #compare
##figure above shows that the MC estimate works well,the ratio is almost equal to 1.

## ------------------------------------------------------------------------
MC.R <- function(x,sigma,n,antithetic) {
n <- numeric(n)
u <- runif(n/2)       #generate random u ~ U(0,1)

if (!antithetic) 
  v <- runif(n/2)     #if not antithentic then generate the other half of n from U(0,1)
else
  v <- 1 - u          #if antithentic then generate the other half of n as 1-u
u <- c(u, v)          #combine the whole u 

cdf <- numeric(length(x))
for (i in 1:length(x)) {
  g <- x[i]^2*u/(sigma^2)*exp(-(u *x[i])^2/(2*sigma^2))  #the estimate of x multipling density of Rayleigh for each x
  cdf[i] <- mean(g)   #the estimate of cdf
}
cdf
}

## ----echo=TRUE-----------------------------------------------------------
set.seed(1234)
sigma <- 1;n <- 1e4
x <- R <- seq(.1, 2.5,length=5) #value the x
for(i in 1:length(x))
  R[i] <- 1-exp(-x[i]^2/(2*sigma^2))   #The theoretical cdf of Rayleigh

  MC1 <- MC.R(x,sigma,n,antithetic = FALSE)  #(x1+x2)/2 with independent x1,x2
  MC2 <- MC.R(x,sigma,n,antithetic = TRUE)   #(x+x`)/2 with negatively correlated x,x`
  
print(round(rbind(x, MC1, MC2, R), 5))
##compare the estimate with or without antithetic variables and the theoretical value,which shows that both MC1 and MC2 fit well while MC2 fits better.

## ------------------------------------------------------------------------
set.seed(1234)
m <- 1e2;sigma <- 1;x <- 6;n <- 1e4
MC1 <- MC2 <- numeric(m)
for (i in 1:m) {
MC1[i] <- MC.R(x,sigma,n,antithetic = FALSE)
MC2[i] <- MC.R(x,sigma,n,antithetic = TRUE)
}
var1 <- var(MC1)
var2 <- var(MC2)
ratio <- (var1-var2)/var1   #the percent of reduction
print(c(var1,var2,ratio))

## ------------------------------------------------------------------------
 x <- seq(1,7,.01)
w<-2
g <- exp(-x^2/2)*x^2/sqrt(2*pi) 
f1 <- 1/(x^2)
f2 <- x*exp((1-x^2)/4)/sqrt(2*pi)
gs <- c(expression(g(x)==e^{-x^2/2}*x^2/sqrt(2*pi)),expression(f[1](x)==1/(x^2)),expression(f[2](x)==x*e^{(1-x^2)/4}/sqrt(2*pi)))
    #for color change lty to col
    par(mfrow=c(1,2))
    #figure (a)
    plot(x, g, type = "l", ylab = "", ylim = c(0,2),lwd = w,col=1,main='(A)')
    lines(x, f1, lty = 2, lwd = w,col=2)
    lines(x, f2, lty = 3, lwd = w,col=3)
    legend("topright", legend = gs,
           lty = 1:3, lwd = w, inset = 0.02,col=1:3)

    #figure (b)
   plot(x, g/f1,  type = "l",ylim = c(0,3.2),ylab = "",lty = 2, lwd = w,col=2)
    lines(x, g/f2, lty = 3, lwd = w,col=3)
    legend("topright", legend = gs[-1],
           lty = 2:3, lwd = w, inset = 0.02,col=2:3)
##The figures above show that f2 is more similiar to g according to figure(A), and it is also more aclinic than f1 in figure(B).f2 is the better estimate than f1.

## ------------------------------------------------------------------------
g <- function(x) {x^2*exp(-x^2/2)/sqrt(2*pi)*(x>1)}
m <- 10000   #size
theta.hat <- se <- numeric(2)
#using inverse transform method to generate random sample f1
u <- runif(m) 
x <- 1/(1-u)
fg <- g(x)*x^2
theta.hat[1] <- mean(fg)
se[1] <- sd(fg)           #sigma for f1

#using inverse transform method to generate random sample f2
u <- runif(m)
x <- (4*(1-log((1-u)*sqrt(2*pi)/2)))^0.5
fg <- g(x)/(x*exp((1-x^2)/4)/sqrt(2*pi))
theta.hat[2] <- mean(fg)
se[2] <- sd(fg)
rbind(theta.hat,se)
##As we have inferred in the figures, the data also prensents that f2 is the better function to get closer to g with a smaller variance.

## ------------------------------------------------------------------------
g <- function(x) {x^2*exp(-x^2/2)/sqrt(2*pi)*(x>1)}
m <- 1e7                  #size
u <- runif(m)             #generate the random u ~U(0,1)
x <- (4*(1-log((1-u)*sqrt(2*pi)/2)))^0.5       #using inverse transform method to generate random sample f
fg <- g(x)/(x*exp((1-x^2)/4)/sqrt(2*pi))
theta.hat <- mean(fg)     #the estimate of E(g(x)/f(x)),also the estimate of integral
theta.hat

## ------------------------------------------------------------------------
GiniEstimation <-  function(cdf,mu,n,m,a,b)  # cdf for the supposed distribution,mu for theoritical expectation,n for size of random numbers,m for repeating times, a and b for the better xlim for histogram(u can choose a better one when u observed its distribution in histogram with random a b )
{
Gcdf<-numeric(m)
set.seed(1234)
v <- numeric(n)

for(i in 1:m)  
  {
  x<-sort(cdf(n))     #sort the x
  
  for(j in 1:n)   
    {
    v[j] <- (2*j-n-1)*x[j]
  } 
  s <- sum(v)          
  Gcdf[i]<-(1/(n^2*mu))*s   #caculate the G for every i
}
mean<-mean(Gcdf)
median<-quantile(Gcdf,0.5)
deciles<-quantile(Gcdf,probs=seq(0,1,0.1))
print(mean)
print(median)
print(deciles)
hist(Gcdf,freq=F,xlim=c(a,b),main="Gini ratio of supposed distribution")  #density histograms of the replicates
}



## ------------------------------------------------------------------------
GiniEstimation(rlnorm,exp(0.5),1e3,1e2,0.45,0.65)
##The mean is much approaching to median as the data shows and we can also roughly conclude from the histogram.

## ------------------------------------------------------------------------
GiniEstimation(runif,0.5,1e3,1e2,0.32,0.35)
##The mean is much approaching to median as the data shows and we can also roughly conclude from the histogram.

## ------------------------------------------------------------------------
n <- m <- 1e3
mu <- 100
Gcdf<-numeric(m)
set.seed(1234)
v <- numeric(n)

for(i in 1:m)  
  {
  x<-sort(rbinom(n,1000,0.1))     #sort the x
  
  for(j in 1:n)   
    {
    v[j] <- (2*j-n-1)*x[j]
  } 
  s <- sum(v)          
  Gcdf[i]<-(1/(n^2*mu))*s   #caculate the G for every i
}
mean<-mean(Gcdf)
median<-quantile(Gcdf,0.5)
deciles<-quantile(Gcdf,probs=seq(0,1,0.1))
print(mean)
print(median)
print(deciles)
hist(Gcdf,freq=F,xlim=c(0.048,0.058),main="Gini ratio of Bernoulli(0.1)")  #density histograms of the replicates
##The mean is much approaching to median as the data shows and we can also roughly conclude from the histogram.

## ------------------------------------------------------------------------
n<-1000              #generate 1000 random numbers
m<-100               #generate 100 gini ratio every times and repeate 100 times
alpha <- 0.05
UCL <- LCL <- numeric(m)
G <- function(b)      #the real value of gamma,b for sigma
  {
 
  return(2*pnorm(b/sqrt(2))-1)
}
CI.sigma<-function(x,n){                     #function to construct confident interval for sigma
  s2<-(1/(n-1))*sum((x-mean(x))^2)
  LCL<- (n-1)*s2/qchisq(1-alpha/2,n-1)
  UCL<-(n-1)*s2/qchisq(alpha/2,n-1)
  c(LCL,UCL)
}
ciG<-cbind(numeric(m),numeric(m))
C<-CR<-numeric(m)
e <- G(1)
set.seed(1234)
for(j in 1:m){
  for(i in 1:m){

  x<-rlnorm(n)
  ci<-CI.sigma(log(x),n)
  ciG[i,1]<-G(ci[1])
  ciG[i,2]<-G(ci[2])
  if(e>ciG[i,1] && e<ciG[i,2])     # for every ith ci,judge whether real value e in the ci
    C[i]<-1
  else
    C[i]<-0
   }

CR[j]<-mean(C)
}
mean(CR)
##The coverage rate of estimation is a little bigger than 0.95.

## ------------------------------------------------------------------------
library(MASS)
n <- 20
m <- 1000
rho.s <- c(seq(-0.8,-0.1,0.1),seq(0.1,0.8,0.1))   #alternatives
M <- length(rho.s)
power.s <- numeric(M)
set.seed(1234)
for (i in 1:M) {
rho <-rho.s[i]
pvalues.s <- replicate(m, expr = {       #simulate under alternative mu1
x <- mvrnorm(n, c(0,0), matrix(c(1,rho,rho,1),2,2))
cortest.s <- cor.test(x[,1],x[,2],
alternative = "two.sided",method="spearman")
cortest.s$p.value } )
power.s[i] <- mean(pvalues.s<= .05)
}
power.s

## ------------------------------------------------------------------------
n <- 20
m <- 1000
rho.k <- c(seq(-0.8,-0.1,0.1),seq(0.1,0.8,0.1))    #alternatives
M <- length(rho.k)
power.k <- numeric(M)
set.seed(1234)
for (i in 1:M) {
rho <-rho.k[i]
pvalues.k <- replicate(m, expr = {       #simulate under alternative mu1
x <- mvrnorm(n, c(0,0), matrix(c(1,rho,rho,1),2,2))
cortest.k <- cor.test(x[,1],x[,2],
alternative = "two.sided",method="kendall")
cortest.k$p.value } )
power.k[i] <- mean(pvalues.k<= .05)
}
power.k

## ------------------------------------------------------------------------
n <- 20
m <- 1000
rho.p <- c(seq(-0.8,-0.1,0.1),seq(0.1,0.8,0.1))       #alternatives
M <- length(rho.p)
power.p <- numeric(M)
set.seed(1234)
for (i in 1:M) {
rho <-rho.p[i]
pvalues.p <- replicate(m, expr = {     #simulate under alternative mu1

x <- mvrnorm(n, c(0,0), matrix(c(1,rho,rho,1),2,2))
cortest.p <- cor.test(x[,1],x[,2],
alternative = "two.sided",method="pearson")
cortest.p$p.value } )
power.p[i] <- mean(pvalues.p<= .05)
}
power.p

## ------------------------------------------------------------------------
power.rho<-rbind(power.s,power.k,power.p)
colnames(power.rho)<-c(seq(-0.8,-0.1,0.1),seq(0.1,0.8,0.1))
power.rho

## ------------------------------------------------------------------------
n <- 20
m <- 1000
set.seed(1234)
pvalues.s1 <- replicate(m, expr = {   #simulate under alternative mu1
x <- rt(n,2)
y<-numeric(n)
for(i in 1:n){
  if(x[i]<0) 
    y[i]<-runif(1,-1,0)
  else       
    y[i]<-runif(1,0,1)
}
cortest.s1 <- cor.test(x,y,
alternative = "two.sided",method="spearman")
cortest.s1$p.value } )
power.s1 <- mean(pvalues.s1<= .05)
power.s1
##nice power

## ------------------------------------------------------------------------
n <- 20
m <- 1000
set.seed(431)
pvalues.p1 <- replicate(m, expr = {
#simulate under alternative mu1
x <- rt(n,2)
y<-numeric(n)
for(i in 1:n){
  if(x[i]<0) y[i]<-runif(1,-1,0)
  else       y[i]<-runif(1,0,1)
}
cortest.p1 <- cor.test(x,y,
alternative = "two.sided",method="pearson")
cortest.p1$p.value } )
power.p1 <- mean(pvalues.p1<= .05)
power.p1

## ------------------------------------------------------------------------
#declare the varibles
library(bootstrap)
n <-nrow(law)   #sample size
theta.jack <- numeric(n)

## ------------------------------------------------------------------------
theta.hat <- cor(law[,1], law[,2]) #law[,1] is LSAT, law[,2] is GPA
b.cor <- function(x,i)cor(x[i,1],x[i,2])
for(i in 1:n)
  {
    theta.jack[i] <- b.cor(law,(1:n)[-i])  #leave one out
}
bias.jack <- (n-1)*(mean(theta.jack)-theta.hat)     #the eatimate bias 
se.jack <- sqrt((n-1)*mean((theta.jack-theta.hat)^2))   #the eatimate standard error
round(c(original=theta.hat,bias=bias.jack,se=se.jack),3)
##Jackknife estimate of the bias and the standard error of the correlation between LSAT and GPA in law data is shown above. 

## ------------------------------------------------------------------------
#declare the varibles
library(boot)
n <- 100 #sample size
m <- 100# times for replicates
ci.norm<-ci.basic<-ci.perc<-ci.bca<-matrix(NA,m,2)

## ------------------------------------------------------------------------

boot.mean <- function(x,i) mean(x$hours[i])
de <- boot(aircondit,statistic=boot.mean,R = 999)
ci <- boot.ci(de,type=c("norm","basic","perc","bca")) #the confidence interval
ci
##The four kinds of CIs have been shown above,and they differ greatly,it is may beacause of the small fixed sample size which is not appproching the normal distribution. 

## ------------------------------------------------------------------------
set.seed(4)
mu <- mean(aircondit$hours)#the estimate of mean time
lambda <- 1/mu
tboot.mean <- function(x,i) mean(x[i])
for (i in 1:m) {
  U<-runif(n);R<--log(1-U)/lambda #generate random sample with expotential distribution
  de <- boot(data=R,statistic=tboot.mean, R = 999)
  ci <- boot.ci(de,type=c("norm","basic","perc","bca"))
  ci.norm[i,]<-ci$norm[2:3];ci.basic[i,]<-ci$basic[4:5]
  ci.perc[i,]<-ci$percent[4:5];ci.bca[i,]<-ci$bca[4:5]
}
cat('norm =',mean(ci.norm[,1]<=mu & ci.norm[,2]>=mu),
    'basic =',mean(ci.basic[,1]<=mu & ci.basic[,2]>=mu),
    'perc =',mean(ci.perc[,1]<=mu & ci.perc[,2]>=mu),
    'BCa =',mean(ci.bca[,1]<=mu & ci.bca[,2]>=mu))
  
##The result above showS that the percentile method has the best effect in expotential distribution with lambda equals the mle of lambda from aircondit data

## ------------------------------------------------------------------------
#declare the varibles
library(bootstrap)
n <- nrow(scor)
m <- ncol(scor)
covmatrix <- matrix(NA,m,m)# for covariance matrix
eigenvalues <- lambda <- numeric(m) #for eigenvalues,sorted eigenvalues(decreasing)


## ------------------------------------------------------------------------
covmatrix <- cov(scor) #covariance matrix
eigenvalues <- eigen(covmatrix)$values #eigenvalues
lambda <- sort(eigenvalues) #sort the eigenvalues
theta.hat <- lambda[1]/sum(lambda) #the empirical theta 
b.theta <- function(x,i) 
             {
               covmatrix <- cov(x[i,])
               eigenvalues <- eigen(covmatrix)$values
               lambda <- sort(eigenvalues)
               theta.hat <- lambda[1]/sum(lambda)
             }
for(i in 1:n)
  {
    theta.jack[i] <- b.theta(scor,(1:n)[-i])   #leave one out
  }
bias.jack <- (n-1)*(mean(theta.jack)-theta.hat)     #the eatimate bias 
se.jack <- sqrt((n-1)*mean((theta.jack-theta.hat)^2))   #the eatimate standard error
round(c(original=theta.hat,bias=bias.jack,se=se.jack),3)
##Jackknife estimate of the bias and the standard error of the theta.hat is shown above. 

## ------------------------------------------------------------------------
#declare the varibles
library(DAAG)
attach(ironslag)
n <- length(magnetic) #in DAAG ironslag
e1 <- e2 <- e3 <- e4 <- numeric(n) # for the prediction error

## ------------------------------------------------------------------------
# for n-fold cross validation
# fit models on leave-two-out samples
for (k in 1:round(n/2))
  {
y <- magnetic[-(2*k-1)] #delete the x_2k-1
y <- y[-2*k]          #delete the x_2k
x <- chemical[-2*k-1]
x <- x[-2*k]

J1 <- lm(y ~ x) #Linear model
yhat1.1 <- J1$coef[1] + J1$coef[2] * chemical[2*k-1]
e1[2*k-1] <- magnetic[2*k-1] - yhat1.1
yhat1.2 <- J1$coef[1] + J1$coef[2] * chemical[2*k]
e1[2*k] <- magnetic[2*k] - yhat1.2

J2 <- lm(y ~ x + I(x^2)) #Quardratic model
yhat2.1 <- J2$coef[1] + J2$coef[2] * chemical[2*k-1] +
J2$coef[3] * chemical[2*k-1]^2
e2[2*k-1] <- magnetic[2*k-1] - yhat2.1
yhat2.2 <- J2$coef[1] + J2$coef[2] * chemical[2*k] +
J2$coef[3] * chemical[2*k]^2
e2[2*k] <- magnetic[2*k] - yhat2.2

J3 <- lm(log(y) ~ x) #Exponential model
logyhat3.1 <- J3$coef[1] + J3$coef[2] * chemical[2*k-1]
yhat3.1 <- exp(logyhat3.1)
e3[2*k-1] <- magnetic[2*k-1] - yhat3.1
logyhat3.2 <- J3$coef[1] + J3$coef[2] * chemical[2*k]
yhat3.2 <- exp(logyhat3.2)
e3[2*k] <- magnetic[2*k] - yhat3.2

J4 <- lm(log(y) ~ log(x)) #Log-Log model
logyhat4.1 <- J4$coef[1] + J4$coef[2] * log(chemical[2*k-1])
yhat4.1 <- exp(logyhat4.1)
e4[2*k-1] <- magnetic[2*k-1] - yhat4.1
logyhat4.2 <- J4$coef[1] + J4$coef[2] * log(chemical[2*k])
yhat4.2 <- exp(logyhat4.2)
e4[2*k] <- magnetic[2*k] - yhat4.2
}
c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))
##The result above shows that the J2 Model 2, the quadratic model,would be the best fit for the data according to the prediction error criterion.


## ------------------------------------------------------------------------
J2
##The related quadratic model is shown above

## ------------------------------------------------------------------------
#declare the varibles
attach(chickwts)    #for data
x <- sort(as.vector(weight[feed == "soybean"]))
y <- sort(as.vector(weight[feed == "linseed"]))
detach(chickwts)

z <- c(x,y) #pool the sample
m <- length(x);n <- length(y);N <- length(z)
S <- numeric(N) #for Cram¡äer-von Mises statistic
R<- 999 #number of replicates
W <- k <- numeric(R) #W for storage of replicates k for permutation

## ------------------------------------------------------------------------
#caculate the Cram¡äer-von Mises statistic
cvm <- function(x1,y1){
Fn <- ecdf(x1);Gm <- ecdf(y1);z1 <- c(x1,y1)
for (i in 1:N) 
  {
    S[i] <- (Fn(z1[i])-Gm(z1[i]))^2 
  }
  m*n/(m+n)^2*sum(S)
  }
W0 <- cvm(x,y)
for (i in 1:R) {
#generate indices k for the first sample
k <- sample(1:N, size = m, replace = FALSE)
x1 <- z[k]
y1 <- z[-k] #complement of x1
W[i] <- cvm(x1,y1)
}
p <- mean(c(W0, W) >= W0)#cacaulate the p value
p
## As the result shows above,it does not support the alternative hypothesis that distributions differ.

## ------------------------------------------------------------------------
#declare the varibles
library(RANN)
library(energy)
library(Ball)
library(boot)
m <- 103; k<-3; p<-2; mu <- 0.5; set.seed(12345)
n1 <- n2 <- 50; R<-999; n <- n1+n2; N = c(n1,n2)

## ------------------------------------------------------------------------
Tn <- function(z, ix, sizes,k) {
n1 <- sizes[1]
n2 <- sizes[2]
n <- n1 + n2
z <- z[ix, ]
o <- rep(0, NROW(z))
z <- as.data.frame(cbind(z, o))
NN <- nn2(z, k=k+1)
block1 <- NN$nn.idx[1:n1,-1]
block2 <- NN$nn.idx[(n1+1):n,-1]
i1 <- sum(block1 < n1 + .5)
i2 <- sum(block2 > n1 + .5)
return((i1 + i2) / (k * n))
}
eqdist.nn <- function(z,sizes,k){
  boot.obj <- boot(data=z,statistic=Tn,R=R,
  sim = "permutation", sizes = sizes,k=k)
  ts <- c(boot.obj$t0,boot.obj$t)
  p.value <- mean(ts>=ts[1])
  list(statistic=ts[1],p.value=p.value)
}
p.values <- matrix(NA,m,3)

#Unequal variances and equal expectations
for(i in 1:m){
  x <- matrix(rnorm(n1*p),ncol=p);
  y <- cbind(rnorm(n2),rnorm(n2,mean=mu));
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value  #NN
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value  #energy method
  p.values[i,3] <-
    bd.test(x=x,y=y,R=R,seed=i*12345)$p.value  #ball method
}
alpha <- 0.05                        #confidence level 
pow <- colMeans(p.values<alpha)
pow

## ------------------------------------------------------------------------
#Unequal variances and unequal expectations
mu <- 0.6
for(i in 1:m){
  x <- matrix(rnorm(n1*p),ncol=p);
  y <- cbind(rnorm(n2),rnorm(n2,mean=mu));
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value  #NN
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value  #energy method
  p.values[i,3] <-
    bd.test(x=x,y=y,R=R,seed=i*12345)$p.value  #ball method
}
alpha <- 0.05                        #confidence level 
pow <- colMeans(p.values<alpha)
pow

## ------------------------------------------------------------------------
#Non-normal distributions: t distribution with 1 df (heavy-tailed distribution), bimodel distribution (mixture of two normal distributions)
for(i in 1:m){
  x <- cbind(rt(n1,1),rnorm(n2, mean = sample(c(0, 1), size = n2, prob = c(0.5, 0.5), replace = T), sd =1 ))
  y <- cbind(rt(n2,1),rnorm(n2, mean = sample(c(1, 2), size = n2, prob = c(0.5, 0.5), replace = T), sd =1 ))
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value  #NN
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value  #energy method
  p.values[i,3] <-
    bd.test(x=x,y=y,R=R,seed=i*12345)$p.value  #ball method
}
alpha <- 0.05                        #confidence level 
pow <- colMeans(p.values<alpha)
pow

## ------------------------------------------------------------------------
#Unbalanced samples (say, 1 case versus 10 controls)
n1 <- 10;n2 <- 1e2;N = c(n1,n2)
for(i in 1:m){
  x <- matrix(rnorm(n1*p),ncol=p);
  y <- cbind(rnorm(n2),rnorm(n2,mean=mu));
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value  #NN
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value  #energy method
  p.values[i,3] <-
    bd.test(x=x,y=y,R=R,seed=i*12345)$p.value  #ball method
}
alpha <- 0.05                        #confidence level 
pow <- colMeans(p.values<alpha)
pow

## ------------------------------------------------------------------------
#declare the varibles
library(stats)
theta <- 1;eta <- 0
m <- 10000
x <- numeric(m)

## ------------------------------------------------------------------------
f <- function(x, theta,eta) {
if (theta > 0) 
return(1/(theta*pi*(1+((x-eta)/theta)^2)))
}

x[1] <- runif(1,min=-1,max=1)
k <- 0
u <- runif(m)
for (i in 2:m) {
xt <- x[i-1]
y <- xt+runif(1,min=-1,max=1) #proposal distribution
num <- f(y, theta,eta) * dnorm(xt, mean = y,sd=0.5)

den <- f(xt, theta,eta) * dnorm(y, mean = xt,sd=0.5)

if (u[i] <= num/den) x[i] <- y 
else {
x[i] <- xt
k <- k+1 #y is rejected
}
}
print(k)

## ------------------------------------------------------------------------
b <- 1001 #discard the burnin sample
y <- x[b:m]
a <- ppoints(100)
QC <- qcauchy(a) #quantiles of Cauchy
Q <- quantile(x, a)
qqplot(QC, Q, main="",xlim=c(-2,2),ylim=c(-2,2),xlab="Cauchy Quantiles", ylab="Sample Quantiles")

## ------------------------------------------------------------------------
hist(y, breaks="scott", main="", xlab="", freq=FALSE)
lines(QC, f(QC,theta,eta))

## ------------------------------------------------------------------------
#declare the varibles
gsize <- c(125,18,20,34)  #group size
m <- 5000
w <- .25 
burn <- 1000

## ------------------------------------------------------------------------
#the target density
prob <- function(y, gsize) {
  if (y < 0 || y >1)
    return (0)
  else
  return((1/2+y/4)^gsize[1] *((1-y)/4)^gsize[2]*((1-y)/4)^gsize[3]*(y/4)^gsize[4])
}

u <- runif(m)  #for accept/reject step
v <- runif(m, -w, w)  #proposal distribution
x[1] <- .25
for (i in 2:m) {
  y <- x[i-1] + v[i]
  if (u[i] <= prob(y, gsize) / prob(x[i-1], gsize))
    x[i] <- y 
  else
    x[i] <- x[i-1]
}
xtheta <- x[(burn+1):m]
theta.hat <- mean(xtheta)
theta.hat
##the estimate of posterior distribution of ¦È is shown above.

## ------------------------------------------------------------------------
#declare the varibles
gsize <- c(125,18,20,34)  #group size
m <- 5000
w <- .25 
burn <- 1000
x <- numeric(m)

## ------------------------------------------------------------------------
#the target density
prob <- function(y, gsize) {
  if (y < 0 || y >1)
    return (0)
  else
  return((1/2+y/4)^gsize[1] *((1-y)/4)^gsize[2]*((1-y)/4)^gsize[3]*(y/4)^gsize[4])
}

u <- runif(m)  #for accept/reject step
v <- runif(m, -w, w)  #proposal distribution
x[1] <- .25
for (i in 2:m) {
  y <- x[i-1] + v[i]
  if (u[i] <= prob(y, gsize) / prob(x[i-1], gsize))
    x[i] <- y 
  else
    x[i] <- x[i-1]
}
x
xtheta <- x[(burn+1):m]
theta.hat <- mean(xtheta)
theta.hat
##the estimate of posterior distribution of theta is shown above.

## ------------------------------------------------------------------------
#declare the varibles
q <- c(-222,-22,-2,0,4,44,444,Inf)
m <-length(q)
CDF <- numeric(m) #for CDF with different p
eta <- 0;theta <- 1 #the location parameter and the scale parameter of Cauchy density

## ------------------------------------------------------------------------
f <- function(x, eta, theta) {
1/(1+((x-eta)/theta)^2)/theta/pi
}

for (i in 1:m ){
  
CDF[i] <- integrate(f,lower=-Inf,upper=q[i],rel.tol=.Machine$double.eps^0.25,eta=eta,theta=theta)$value
}
pcauchy <- pcauchy(q,location = eta, scale = theta)
data.frame(q,pcauchy,CDF)
##As the result shows above, the caculated CDF is the same to the results of pcauchy function with different q.

## ------------------------------------------------------------------------
#declare the varibles
nA <- 28;nB <- 24;nAB <- 70;nOO <- 41 #observed data
L <- c(0.3,0.2) #initial values for p q 
tol <- .Machine$double.eps^0.5 # for caculation accuracy
n <- 10000 #max. number of iterations
EI <- M <- numeric(n)# for the values of "M-step"

## ------------------------------------------------------------------------
p <- L[1];q <- L[2]   #initial values
for (i in 1:n) {
  a<-p/(2-p-2*q)
  b<-q/(2-2*p-q)
  fp=nA*(1+a)+nAB
  fq=nB*(1+b)+nAB
  fpq=nA*(1-a)+nB*(1-b)+2*nOO
  
  M[i]<-fp*log(p)+fpq*log(1-p-q)+fq*log(q) #Record the    maximum likelihood values in M-steps
  
  x<-p; y<-q                #store ith values of p and q
  p<-fp/(fp+fq+fpq); q<-fq/(fp+fq+fpq) #update the parameters
  #k<-k+1     #the times until converge
  if (abs(p-x)/x<tol && abs(q-y)/y<tol) break  #control the iterations
  p.hat<-p; q.hat<-q
 
}
EI<-M[1:i]#the values of "M-step"
data.frame(p.hat,q.hat)
##the estimates of p q is shown above.


## ------------------------------------------------------------------------
i ##the times until converge

## ------------------------------------------------------------------------
-EI 

## ------------------------------------------------------------------------
x <- seq(1,i,1)
plot(x,EI)
##the maximum likelihood values in M-steps are  increasing

## ------------------------------------------------------------------------
formulas <- list(
mpg ~ disp,
mpg ~ I(1 / disp),
mpg ~ disp + wt,
mpg ~ I(1 / disp) + wt
)

## ------------------------------------------------------------------------
 out_formulas <- vector('list', length(formulas))
   for(i in seq_along(formulas)) {
       out_formulas[[i]] <- lm(formulas[[i]], data = mtcars)
   }
   out_formulas

## ------------------------------------------------------------------------
lapply(formulas, lm, data = mtcars)

## ------------------------------------------------------------------------
bootstraps <- lapply(1:10, function(i) {
rows <- sample(1:nrow(mtcars), rep = TRUE)
mtcars[rows, ]
})

## ------------------------------------------------------------------------
out_bootstrap <- vector('list', length(bootstraps))
   for(i in seq_along(bootstraps)) {
       out_bootstrap[[i]] <- lm(mpg~disp, data = bootstraps[[i]])
   }
   out_bootstrap

## ------------------------------------------------------------------------
 lapply(bootstraps, lm, formula = mpg~disp)

## ------------------------------------------------------------------------
rsq <- function(mod) summary(mod)$r.squared

## ------------------------------------------------------------------------
lapply(out_formulas, rsq)

## ------------------------------------------------------------------------
lapply(out_bootstrap, rsq)

## ------------------------------------------------------------------------
trials <- replicate(
100,
t.test(rpois(10, 10), rpois(7, 10)),
simplify = FALSE
)

## ------------------------------------------------------------------------
sapply(trials, function(mod) mod$p.value)

## ------------------------------------------------------------------------
sapply(trials, `[[`, 'p.value')

## ------------------------------------------------------------------------
 library(parallel)
   mcvMap <- function(f, FUN.VALUE , ...) {
       out <- mcMap(f, ...)
       vapply(out, identity, FUN.VALUE)
   }
   

## ------------------------------------------------------------------------
expected <- function(colsum, rowsum, total) {
  (colsum / total) * (rowsum / total) * total
}                      #for $n_{i??}n_{??j}/n$

chi_stat <- function(observed, expected) {
  ((observed - expected) ^ 2) / expected
}


chisq_test2 <- function(x, y) {#x for the first row, y for the second row
  total <- sum(x) + sum(y)
  rowsum_x <- sum(x)
  rowsum_y <- sum(y)
  chistat <- 0
  for (i in seq_along(x)) {
    colsum <- x[i] + y[i]
    expected_x <- expected(colsum, rowsum_x, total)
    expected_y <- expected(colsum, rowsum_y, total)
    chistat <- chistat + chi_stat(x[i], expected_x) 
    chistat <- chistat + chi_stat(y[i], expected_y)
  }
  chistat
}


## ------------------------------------------------------------------------
#examples of the use of the new chisq_test2

print(chisq_test2(seq(3, 8), seq(4, 9)))

print(chisq.test(seq(3, 8), seq(4, 9)))

##there is a doubt that the value of chisq from the new chisq_test is much different from the chisq.test function


## ------------------------------------------------------------------------
table2 <- function(x, y) {
  x_val <- unique(x)
  y_val <- unique(y)
  mat <- matrix(0L, length(x_val), length(y_val))
  for (i in seq_along(x)) {
    mat[which(x_val == x[[i]]), which(y_val == y[[i]])] <-
      mat[which(x_val == x[[i]]),  which(y_val == y[[i]])] + 1L
  }
  dimnames <- list(x_val, y_val)
  names(dimnames) <- as.character(as.list(match.call())[-1])  # R has names for dimnames... :/
  tab <- array(mat, dim = dim(mat), dimnames = dimnames)
  class(tab) <- "table"
  tab
}


## ------------------------------------------------------------------------
#examples of the use of the new table2()
a <- c(1, 2, 3)
identical(table(a, a), table2(a, a))


b <- c(2, 3, 4)
identical(table(a, b), table2(a, b))


c <- c(1, 2, 3, 1, 2, 3)
d <- c(2, 3, 4, 2, 3, 4)
identical(table(c, d), table2(c, d))


e <- c(1, 2, 2)
identical(table(a, e), table2(a, e))

identical(table(b, e), table2(b, e))

identical(table(e, e), table2(e, e))


f <- c(1, 1, 1)
identical(table(f, f), table2(f, f))

identical(table(e, f), table2(e, f))


g <- c(1, 4, 9)
identical(table(g, g), table2(g, g))

identical(table(g, f), table2(g, f))




## ------------------------------------------------------------------------
expected <- function(colsum, rowsum, total) {
  (colsum / total) * (rowsum / total) * total
}                      #for $n_{i??}n_{??j}/n$

chi_stat <- function(observed, expected) {
  ((observed - expected) ^ 2) / expected
}
x <- y <- NULL
z <- table2(x,y)
chisq_test3 <- function(z) {#x for the first row, y for the second row
  total <- sum(x) + sum(y)
  rowsum_x <- sum(x)
  rowsum_y <- sum(y)
  chistat <- 0
  for (i in seq_along(x)) {
    colsum <- x[i] + y[i]
    expected_x <- expected(colsum, rowsum_x, total)
    expected_y <- expected(colsum, rowsum_y, total)
    chistat <- chistat + chi_stat(x[i], expected_x) 
    chistat <- chistat + chi_stat(y[i], expected_y)
  }
  chistat
}


## ------------------------------------------------------------------------
x <- seq(3, 8)
y <- seq(4, 9)
system.time(print((chisq_test2(x, y))))

system.time(print(chisq_test3(table2(x, y))))
##we can see that the new table function dose speed up the new chisq_test function

