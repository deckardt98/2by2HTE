setwd("/Users/deckard/Desktop/Fan Li Project/Factorial HTE/z-test with interaction")
#f(m)
f <- function(icc,m){
  return((icc*(m-1)+1))
}
#function approximating p0 to p3
approxp <- function(m, rho00, rho0, sigma0, b4, piz){
  #Input
  #m: cluster size
  #rho00: rho_0,0
  #rho0: rho_0
  #sigma0: sigma_0^@
  #b4: sigma_b_4^2
  #piz: pi_Z
  
  #output:
  #app: vector of (p0,p1,p2,p3)
  
  #sigma_epsilon^2
  ep <- sigma0*(1-rho0)
  #mean of m1i
  mb <- m*piz
  #sigma_Z^2
  sZ <- piz*(1-piz)
  
  #determinant
  D <- f(rho00,m)^2*b4^2+4*b4*rho00*ep*f(rho00,m)
  #roots
  x1 <- (f(rho00,m)*b4-sqrt(D))/(2*b4*rho00)
  x2 <- (f(rho00,m)*b4+sqrt(D))/(2*b4*rho00)
  #factors
  f1 <- ((x1+x2)*(x1+x2-3*mb)-x1*x2+3*mb^2)/((x1-mb)^3*(x2-mb)^3)
  f2 <- (x1+x2-2*mb)/((x1-mb)^2*(x2-mb)^2)-2*mb*f1
  f3 <- 3*mb^2*f1+(x1*x2-mb^2)/((x1-mb)^2*(x2-mb)^2)-
    2*mb*(x1*x2*(x1+x2-3*mb)+mb^3)/((x1-mb)^3*(x2-mb)^3)
  #third moment for m1i
  tm <- m*(m-1)*(m-2)*piz^3+3*m*(m-1)*piz^2+m*piz
  if (b4==0){
    app0 <- 1/(ep*f(rho00,m))
    app1 <- mb*app0
    app2 <- (m*sZ+mb^2)*app0
    app3 <- tm*app0
  }else {
    #approximated p0
    app0 <- -1/(b4*rho00)*(f1*(m*sZ+mb^2)+f2*mb+f3)
    #approximated p1
    app1 <- -1/(b4*rho00)*(f1*tm+f2*(m*sZ+mb^2)+f3*mb)
    #approximated p2
    app2 <- -1/(b4*rho00)+f(rho00,m)/rho00*app1+ep*f(rho00,m)/(b4*rho00)*app0
    #approximated p3
    app3 <- -mb/(b4*rho00)+f(rho00,m)/rho00*app2+ep*f(rho00,m)/(b4*rho00)*app1
  }
  #vector output
  app <- c(app0,app1,app2,app3)
  return(app)
}

#function calculating lambda11 lambda12 lambda22
callamb <- function(m, rho00, rho0, sigma0, b4, piz){
  #Input
  #m: cluster size
  #rho00: rho_0,0
  #rho0: rho_0
  #sigma0: sigma_0^@
  #b4: sigma_b_4^2
  #piz: pi_Z
  
  #output:
  #estlamb lambda11,lambda12,lambda22
  
  #sigma_epsilon^2
  ep <- sigma0*(1-rho0)
  esp <- approxp(m, rho00, rho0, sigma0, b4, piz)
  p0 <- esp[1]
  p1 <- esp[2]
  p2 <- esp[3]
  p3 <- esp[4]
  
  #calculate lambda11 lambda12 lambda22
  lamb11 <- (m+b4*(2*m*rho00-f(rho00,m))*p2-rho00*m^2*b4*p1-rho00*m^2*ep*p0)/(ep)
  lamb12 <- (m*piz+rho00*b4*p3-f(rho00,m)*b4*p2-m*rho00*ep*p1)/(ep)
  lamb22 <- (m*piz+rho00*b4*p3-(f(rho00,m)*b4+rho00*ep)*p2)/(ep)
  estlamb <- c(lamb11,lamb12,lamb22)
  return(estlamb)
}

#function calculating the number of cluster needed 
nc <- function(type, delta, m, rho00, rho0, sigma0, b4, piz, pix){
  #Input
  #type: test type 1=cex 2=cez 3=mex 4=mez 5=ie
  #delta: effect size
  #m: cluster size
  #rho00: rho_0,0
  #rho0: rho_0
  #sigma0: sigma_0^2
  #b4: sigma_b_4^2
  #piz: pi_Z
  #pix: pi_X
  
  #output:
  #n: the number of cluster 
  
  #type i error rate
  alpha <- 0.05
  #1-power
  mu <- 0.2
  #sigma_Z^2
  sz <- piz*(1-piz)
  #calculate lambda
  lamb <- callamb(m, rho00, rho0, sigma0, b4, piz)
  l11 <- lamb[1]
  l12 <- lamb[2]
  l22 <- lamb[3]
  if (type==1){
    #variance part 1
    vp1 <- (sigma0*((m-1)*rho0+1)*((m-1)*rho0*(1-piz)+1-rho0))/(m*(1-pix)*(1-piz)*((m-2)*rho0+1))
    #variance part2
    vp2 <- l22/(pix*(l11*l22-l12^2))
  }else if (type==2){
    #variance part 1
    vp1 <- (sigma0*((m-1)*rho0+1)*(1-rho0))/(m*(1-pix)*sz*((m-2)*rho0+1))
    #variance part2
    vp2 <- 0
  }else if (type==3){
    #variance part 1
    vp1 <- (sigma0*((m-1)*rho0+1))/(m*(1-pix))
    #variance part2
    vp2 <- (piz^2*l11-2*piz*l12+l22)/(pix*(l11*l22-l12^2))
  } else if (type==4){
    #variance part 1
    vp1 <- (sigma0*(1-rho0)*((m-1)*rho0+1)*(1-pix))/(m*sz*((m-2)*rho0+1))
    #variance part2
    vp2 <- (pix*l11)/(l11*l22-l12^2)
  }else if (type==5){
    #variance part 1
    vp1 <- (sigma0*(1-rho0)*((m-1)*rho0+1))/(m*(1-pix)*sz*((m-2)*rho0+1))
    #variance part2
    vp2 <- l11/(pix*(l11*l22-l12^2))
  }
  n <- (qnorm(1-alpha/2)+qnorm(1-mu))^2/delta^2*(vp1+vp2)
  return(n)
}

#function calculating the number of cluster needed for reduced model
nc_tian <-  function(type, delta, m, rho0, sigma0, piz, pix){
  #Input
  #type: test type 1=cex 2=cez 3=mex 4=mez 5=ie
  #delta: effect size
  #m: cluster size
  #rho0: rho_0
  #sigma0: sigma_0^2
  #piz: pi_Z
  #pix: pi_X
  
  #output:
  #n: the number of cluster 
  
  #type i error rate
  alpha <- 0.05
  #1-power
  mu <- 0.2
  #sigma_X^2
  sx <- pix*(1-pix)
  #sigma_Z^2
  sz <- piz*(1-piz)
  #variance 
  if (type==1){
    v <- sigma0*(1+(m-1)*rho0)*(rho0*(m-1)*(1-piz)+(1-rho0))/m/(1+(m-2)*rho0)/(1-piz)/sx
  }else if (type==2){
    v <- (sigma0*((m-1)*rho0+1)*(1-rho0))/(m*(1-pix)*sz*((m-2)*rho0+1))
  }else if (type==3){
    v <- sigma0*((m-1)*rho0+1)/m/sx
  } else if (type==4){
    v <- sigma0*(1-rho0)*(rho0*(m-1)+1)/m/sz/(rho0*(m-2)+1)
  }else if (type==5){
    v <- sigma0*(1-rho0)*(rho0*(m-1)+1)/m/sz/sx/(rho0*(m-2)+1)
  }
  n <- (qnorm(1-alpha/2)+qnorm(1-mu))^2/delta^2*v
  return(n)
}

calnctian <- function(z){
  nctian <- c()
  for (i in 1:nrow(z)){
    temp <- ceiling(nc_tian(z[i,1], z[i,11], z[i,10], z[i,7], 1, 0.5, 0.5))
    if (temp%%2==1){temp <- temp+1}
    nctian[i] <- temp
  }
  return(nctian)
}
z1 <- read.csv("z-test for CEX updated.csv")
z1 <- z1[,-c(1)]
z2 <- read.csv("z-test for CEZ updated.csv")
z2 <- z2[,-c(1)]
z3 <- read.csv("z-test for MEX updated.csv")
z3 <- z3[,-c(1)]
z4 <- read.csv("z-test for MEZ updated.csv")
z4 <- z4[,-c(1)]
z5 <- read.csv("z-test for IE updated.csv")
z5 <- z5[,-c(1)]




z1$nc.tian <- calnctian(z1)
z2$nc.tian <- calnctian(z2)
z3$nc.tian <- calnctian(z3)
z4$nc.tian <- calnctian(z4)
z5$nc.tian <- calnctian(z5)

z1 <- z1[,c(1,12,13,2,3,4,5,7,8,9,10,11)]
#write.csv(z1,"z-test CEX FINAL.csv")
z2 <- z2[,c(1,12,13,2,3,4,5,7,8,9,10,11)]
#write.csv(z2,"z-test CEZ FINAL.csv")
z3 <- z3[,c(1,12,13,2,3,4,5,7,8,9,10,11)]
#write.csv(z3,"z-test MEX FINAL.csv")
z4 <- z4[,c(1,12,13,2,3,5,4,7,8,9,10,11)]
#write.csv(z4,"z-test MEZ FINAL.csv")
z5 <- z5[,c(1,12,13,2,3,5,4,7,8,9,10,11)]
#write.csv(z5,"z-test IE FINAL.csv")

library(latex2exp)
library(ggplot2)
library(gridExtra)

plot <- function(z,label,xl,yl){
  require(latex2exp)
  require(ggplot2)
  z$m <- factor(z$m)
  z$rho0 <- factor(z$rho0)
  p1 <- ggplot(z, aes(x=number.of.cluster, y=nc.tian, shape=m, color=m)) + 
    geom_point(size=2.5) + 
    geom_abline(intercept = 0, slope = 1, lty=2, alpha=0.5) +
    xlab(xl) +
    ylab(yl) +
    ggtitle(TeX(label)) +
    xlim(0,205) +
    ylim(0,205) +
    theme_bw(base_size = 24) +
    scale_color_manual(labels = c("20", "50", "100"), values = c("navy", "orange", "red")) 
  return(p1)
}
p1 <- plot(z1,"(A): Testing $CE_X$",TeX("$n_{CE_X}$"),TeX("$n_{CE_X}^{*}$"))
p1
p2 <- plot(z2,"(B): Testing $CE_Z$",TeX("$n_{CE_Z}$"),TeX("$n_{CE_Z}^{*}$"))
p3 <- plot(z3,"(C): Testing $ME_X$",TeX("$n_{ME_X}$"),TeX("$n_{ME_X}^{*}$"))
p4 <- plot(z4,"(D): Testing $ME_Z$ or $IE$",TeX("$n_{ME_Z}$ or $n_{IE}$"),TeX("$n_{ME_Z}^{*}$ or $n_{IE}^{*}$"))
#p5 <- plot(z5,"(E): Test for $IE$")

#grid.arrange(p1,p2,p3,p4,p5,nrow=3,ncol=2)
grid.arrange(p1,p2,p3,p4,nrow=2,ncol=2)
#function calculating the number of cluster needed without interaction
nc_wi <- function(type, delta, m, rho0, rho1, sigma0, piz, pix){
  #Input
  #type: test type 1=ex 2=ez
  #delta: effect size
  #m: cluster size
  #rho0: rho_0
  #rho1: rho_1
  #sigma0: sigma_0^2
  #piz: pi_Z
  #pix: pi_X
  
  #output:
  #n: the number of cluster 
  
  #type i error rate
  alpha <- 0.05
  #1-power
  mu <- 0.2
  #sigma_Z^2
  sz <- piz*(1-piz)
  #sigma_1^2
  sigma1 <- sigma0*(1-rho0)/(1-rho1)
  
  if (type==1){
    #variance
    v <- sigma0*f(rho0,m)/(1-pix)/m + sigma1*f(rho1,m)/pix/m
  }else if (type==2){
    #variance
    v <- 1/m/sz/((1-pix)*f(rho0,(m-1))/sigma0/(1-rho0)/f(rho0,m)+
                   pix*f(rho1,(m-1))/sigma1/(1-rho1)/f(rho1,m))
  }
  n <- (qnorm(1-alpha/2)+qnorm(1-mu))^2/delta^2*v
  return(n)
}

calncwi <- function(z,type){
  ncwi <- c()
  for (i in 1:nrow(z)){
    temp <- ceiling(nc_wi(type, z[i,11], z[i,10], z[i,7], z[i,8],1, 0.5, 0.5))
    if (temp%%2==1){temp <- temp+1}
    ncwi[i] <- temp
  }
  return(ncwi)
}
z1 <- read.csv("z-test for CEX updated.csv")
z1 <- z1[,-c(1)]
z2 <- read.csv("z-test for CEZ updated.csv")
z2 <- z2[,-c(1)]
z3 <- read.csv("z-test for MEX updated.csv")
z3 <- z3[,-c(1)]
z4 <- read.csv("z-test for MEZ updated.csv")
z4 <- z4[,-c(1)]
z5 <- read.csv("z-test for IE updated.csv")
z5 <- z5[,-c(1)]




z1$nc.wi <- calncwi(z1,1)
z2$nc.wi <- calncwi(z2,2)
z3$nc.wi <- calncwi(z3,1)
z4$nc.wi <- calncwi(z4,2)

library(latex2exp)
library(ggplot2)
library(gridExtra)

plot <- function(z,xlab,ylab,title){
  require(latex2exp)
  require(ggplot2)
  z$m <- factor(z$m)
  z$rho0 <- factor(z$rho0)
  p1 <- ggplot(z, aes(x=number.of.cluster, y=nc.wi, shape=m, color=m)) + 
    geom_point(size=2.5) + 
    geom_abline(intercept = 0, slope = 1, lty=2, alpha=0.5) +
    xlab(TeX(xlab)) +
    ylab(TeX(ylab)) +
    ggtitle(TeX(title)) +
    xlim(0,205) +
    ylim(0,205) +
    theme_bw(base_size = 24) +
    scale_color_manual(labels = c("20", "50", "100"), values = c("navy", "orange", "red")) 
  return(p1)
}

p1 <- plot(z1,"$n_{CE_X}$","$n_{E_X}$","(A): Testing $CE_X$")
p1
p2 <- plot(z2,"$n_{CE_Z}$","$n_{E_Z}$","(B): Testing $CE_Z$")
p3 <- plot(z3,"$n_{ME_X}$","$n_{E_X}$","(C): Testing $ME_X$")
p4 <- plot(z4,"$n_{ME_Z}$","$n_{E_Z}$","(D): Testing $ME_Z$")

grid.arrange(p1,p2,p3,p4,nrow=2,ncol=2)


#ratio of our model vs. reduced model 
ratio_our_reduce <- function(type, delta, m, rho00, rho0, sigma0, b4, piz, pix){
  ratio <- nc(type, delta, m, rho00, rho0, sigma0, b4, piz, pix)/nc_tian(type, delta, m, rho0, sigma0, piz, pix)
  return(ratio)
}

#plot heat map 
phm <- function(type,delta,rho0,m,label,maxl){
  require(latex2exp)
  x <- seq(rho0,0.2,by=0.001)
  y <- seq(rho0,0.2,by=0.001)
  z <- matrix(data= NA, nrow=length(x),ncol=length(y))
  for (i in 1:length(x)){
    for (j in 1:length(y)){
      if (x[i]>y[j]){
        rho11 <- x[i]
        rho00 <- y[j]
      } else {
        rho00 <- x[i]
        rho11 <- y[j]
      }
      b4 <- (1-rho0)*(rho11-rho00)/(1-rho00)/(1-rho11)
      z[i,j] <- ratio_our_reduce(type, delta, m, rho00, rho0, 1, b4, 0.5, 0.5)
    }
  }
  
  filled.contour(x=x,
                 y=y,
                 z=z,
                 levels = seq(0,maxl,length.out = 100),
                 plot.title = title(TeX(label), xlab=TeX("$\\rho_{0,0}$"),
                                    ylab=TeX("$\\rho_{1,1}$")),
                 plot.axes = {
                   axis(1, at=seq(rho0,0.20,by=0.01));
                   axis(2, at=seq(rho0,0.20,by=0.01));
                   contour(x,y,z,nlevels =20,add=TRUE,labcex = 1,col="black")
                 },
                 key.title = {par(cex.main=1.2);title(main="ratio")}
                 #key.axes = axis(4, at=seq(0.5,1.3,by=0.1))
  )
}
pix <- piz <- 0.5
type <- 5
delta1 <- 0.30
delta2 <- 0.25
delta3 <- 0.30
delta4 <- 0.18
delta5 <- 0.30
m <- c(20,50,100)
sigma0 <- 1
rho0 <- c(rep(0.01,10), rep(0.10,3))
rho00 <- c(c(rep(0.01,4), rep(0.05,3), rep(0.10,2), 0.20),
           c(rep(0.10,2), rep(0.20,1)))
rho11 <- c(c(0.01, 0.05, 0.10, 0.20, 0.05, 0.10, 0.20, 0.10, 0.20, 0.20),
           c(0.10,0.20,0.20))
dftemp <- NULL
dff <- NULL
for (k in 1:3){

    for (j in 1:13){
      b4 <- (1-rho0[j])*(rho11[j]-rho00[j])/(1-rho00[j])/(1-rho11[j])
      n1 <- ceiling(nc(type, delta5, m[k], rho00[j], rho0[j], sigma0, b4, piz, pix))
      print(n1)
      n1 <- n1+n1%%2
      n2 <- ceiling(nc_tian(type, delta5, m[k], rho0[j], sigma0, piz, pix))
      n2 <- n2+n2%%2
      n3 <- ceiling(nc_wi(2, delta4, m[k], rho0[j], rho00[j], sigma0, piz, pix))
      n3 <- n3+n3%%2
      dftemp <- rbind(dftemp,c(rho0[j],rho00[j],rho11[j],n1,n2,n3))
    }
}
dftemp
dff <- cbind(dftemp[1:13,],dftemp[14:26,c(4,5,6)],dftemp[27:39,c(4,5,6)])
dff
setwd("/Users/deckard/Desktop/Fan Li Project/Project 2")
write.csv(dff, file="DataExample_MEZ.csv")

nc(type, delta, m, rho00, rho0, sigma0, b4, piz, pix)
ratio_our_reduce(type, delta1, m[3], 0)



phm(4,delta4,0.01,m[1],"(a): m=20",4)
phm(4,delta4,0.01,m[2],"(b): m=50",4)
phm(4,delta4,0.01,m[3],"(c): m=100",4)

nc_tian(1, delta1, m[1], 0.01, 1, 0.5, 0.5)
nc_tian(1, delta1, m[2], 0.01, 1, 0.5, 0.5)
nc_tian(1, delta1, m[3], 0.01, 1, 0.5, 0.5)

nc_tian(2, delta2, m[1], 0.01, 1, 0.5, 0.5)
nc_tian(2, delta2, m[2], 0.01, 1, 0.5, 0.5)
nc_tian(2, delta2, m[3], 0.01, 1, 0.5, 0.5)

nc_tian(3, delta3, m[1], 0.01, 1, 0.5, 0.5)
nc_tian(3, delta3, m[2], 0.01, 1, 0.5, 0.5)
nc_tian(3, delta3, m[3], 0.01, 1, 0.5, 0.5)

nc_tian(4, delta4, m[1], 0.01, 1, 0.5, 0.5)
nc_tian(4, delta4, m[2], 0.01, 1, 0.5, 0.5)
nc_tian(4, delta4, m[3], 0.01, 1, 0.5, 0.5)

nc_tian(5, delta5, m[1], 0.01, 1, 0.5, 0.5)
nc_tian(5, delta5, m[2], 0.01, 1, 0.5, 0.5)
nc_tian(5, delta5, m[3], 0.01, 1, 0.5, 0.5)
