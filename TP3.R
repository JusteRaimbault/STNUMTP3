#comparison of estimators for position parameter

n=10000;
X=rnorm(n,mean=5,sd=2);
par(bg="cornsilk",lwd=2,col="darkblue")
hist(X,breaks=40,freq=F,col="cyan")
curve(dnorm(x,mean=5,sd=2),add=T)


location_estimator=function(U,theta){
  n=length(U)
  theta1=cumsum(U)/(1:n)
  theta2=1:n;theta3=1:n;
  for (i in 1:n){
    theta2[i]=median(U[1:i])
    theta3[i]=(min(U[1:i])+max(U[1:i]))/2
  }
  par(mfrow=c(1,3),bg="cornsilk",lwd=2,col="darkblue")
  plot(1:n,theta1,type="l",main="moyenne")
  abline(h=theta,col="darkred")
  plot(1:n,theta2,type="l",main="mediane")
  abline(h=theta,col="darkred")
  plot(1:n,theta3,type="l",main="Minmax")
  abline(h=theta,col="darkred")
  #return(matrix(c(theta1,theta2,theta3),n,3))
}

location_estimator(rnorm(100000,mean=5,sd=2),5)

compare_estimators=function(N,n){
  par(mfrow=c(1,1))
  X=matrix(rcauchy(N*n,location=5,scale=2),N,n)
  theta1=apply(X,1,mean)
  theta2=apply(X,1,median)
  #theta3=(apply(X,1,min)+apply(X,1,max))/2
  par(bg="cornsilk",lwd=2,col="darkblue")
  boxplot(theta1,theta2,col="cyan")
}

compare_estimators(200,1000)




check_student=function(N,n,mu,sigma)
{
  X=matrix(rnorm(N*n,mean=mu,sd=sigma),N,n)
  t=sqrt(n)*(apply(X,1,mean)-mu)/apply(X,1,sd)
  return(t)
}

N=20000; n=1000;
t=check_student(N,n,100,50)
par(bg="cornsilk",lwd=2,col="darkblue")
hist(t,breaks=100,freq=F,col="cyan")
curve(dt(x,n-1),add=T,lwd=2)

quantile(x=t,probs=c(0.025,0.975))

quants = function(n,Nrange,ind){
  res=c()
  for(N in Nrange){
    res = append(res,quantile(x=check_student(N,n,100,50),probs=c(0.025,0.975))[[1]])
  }
  return(res)
}

plot(x=seq(from=5000,to=50000,by=1000),y=quants(1000,seq(from=5000,to=50000,by=1000),1),type="l")

X=c(19.6, 19.9, 20.4, 19.8, 20.5, 21.0, 18.5, 19.7, 18.4, 19.4)
n=10

mean(X) - quantile(x=check_student(20000,10,mean(X),sd(X)),probs=c(0.025,0.975)) * sd(X)/sqrt(10)


#X=c(77.551, 45.195, 50.626, 39.878, 29.137, 57.321, 39.140, 66.776, 48.028, 42.325, 31.200, 38.632, 42.914, 60.969, 22.076, 52.446, 45.257, 42.626, 62.504, 22.684, 69.196, 42.383, 61.339, 45.803, 74.707, 33.048, 72.423, 43.670, 65.279, 42.714, 59.785, 101.742, 59.641, 44.749, 44.161, 58.488, 46.448, 25.280, 67.619, 66.846, 80.208, 98.492, 41.149, 40.395, 22.220, 34.628, 77.768, 48.161,48.909, 66.267)

X = c(-7.54,  82.51,  14.27,   3.96,  189.98,   17.20,  -20.07,   52.66,   93.47,
-33.57,  13.13,  -1.26,  12.69,   53.33,    2.85,   -7.25,   13.30,   -5.67,
-38.99,  24.24,   4.17,  12.30,   21.59,   -6.70,    1.24,   13.91,   30.24,
3.35,   6.45, -26.22,  72.65,   10.12,   -1.64,   21.49,  391.11,   26.53,
146.60,   2.11,   5.84,  14.25,    7.17,    4.96,   -9.55,    7.89,   -2.31,
91.11,   8.39,   6.23,  25.45,    9.36,  102.44,   -7.28,  -40.02,   -8.86,
14.11,   6.84, -11.15,  -6.67,  -84.82, -241.41,   -0.14,  -72.95,   21.09,
53.47,  -3.80, -10.64,  19.71,   45.89, -124.30,   -2.02,   -1.67,    7.81,
-9.76,   6.25,  16.68,   8.88,   32.14,    1.29,  -10.00,   -5.03,  -66.77,
12.85,  15.32,  31.27,   6.59,    3.92,    8.61,   15.38,   -1.34,   14.11,
10.53,   2.35, -94.19,  16.45,    2.97,   12.26,    4.15,   10.63,    5.47)

ll=function(a=0.5, sigma=1.5)
{
  if(a > 0 && sigma > 0)
    -sum(dcauchy(X, location=a, scale=sigma, log=TRUE))
  else
    NA
}

library(stats4)
fit = mle(ll)
summary(fit)

vcov(fit)
par(mfrow=c(2,1),bg="cornsilk",col="blue",lwd=2)
plot(profile(fit), absVal=FALSE)
confint(fit,level=0.95)


# Calcul des valeurs de la log-vraisemblance
K=80
x=(1:K)/4;  y=(1:K)/4;  z=c();
for (i in 1:length(x))
  for (j in 1:length(y))
    z=c(z,ll(x[i],y[j]))
# Transformation des valeurs calculées
z=matrix(z,length(x),length(y))
z=log(0.001+((z-min(z))/(max(z)-min(z))))
# Le contenu des 7 lignes suivantes peut etre utilisé comme une
# boîte noire
nrz <- nrow(z)
ncz <- ncol(z)
jet.colors <- colorRampPalette( c("blue", "green") )
nbcol <- 100
color <- jet.colors(nbcol)
zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
facetcol <- cut(zfacet, nbcol)
# Visualisation des résultats
par(bg="cornsilk",lwd=1,mfrow=c(1,1))
#image(x,y,z,col = cm.colors(50))
#contour(x,y,z,add=T,col="darkred")
persp(x, y, z,ticktype="detailed",expand=0.5,col=color[facetcol],shade=0.4)

X = rgamma(1000,shape=5,scale=3)
fit = mle(ll)
summary(fit)
confint(fit,level=0.95)

