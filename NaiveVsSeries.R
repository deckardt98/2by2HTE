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
  #approximated p0
  app0 <- -1/(b4*rho00)*(f1*(m*sZ+mb^2)+f2*mb+f3)
  #third moment for m1i
  tm <- m*(m-1)*(m-2)*piz^3+3*m*(m-1)*piz^2+m*piz
  #approximated p1
  app1 <- -1/(b4*rho00)*(f1*tm+f2*(m*sZ+mb^2)+f3*mb)
  #approximated p2
  app2 <- -1/(b4*rho00)+f(rho00,m)/rho00*app1+ep*f(rho00,m)/(b4*rho00)*app0
  #approximated p3
  app3 <- -mb/(b4*rho00)+f(rho00,m)/rho00*app2+ep*f(rho00,m)/(b4*rho00)*app1
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


# function calculating p0-p3 using naive approach 
approxp_na <- function(m, rho00, rho0, sigma0, b4, piz){
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
  
  #approximated p0
  app0 <- 1/(-b4*rho00*mb^2+f(rho00,m)*b4*mb+ep*f(rho00,m))
  #third moment for m1i
  tm <- m*(m-1)*(m-2)*piz^3+3*m*(m-1)*piz^2+m*piz
  #approximated p1
  app1 <- mb*app0
  #approximated p2
  app2 <- -1/(b4*rho00)+f(rho00,m)/rho00*app1+ep*f(rho00,m)/(b4*rho00)*app0
  #approximated p3
  app3 <- -mb/(b4*rho00)+f(rho00,m)/rho00*app2+ep*f(rho00,m)/(b4*rho00)*app1
  #vector output
  app <- c(app0,app1,app2,app3)
  return(app)
  
}

#function calculating lambda11 lambda12 lambda22 via naive approach
callamb_na <- function(m, rho00, rho0, sigma0, b4, piz){
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
  esp <- approxp_na(m, rho00, rho0, sigma0, b4, piz)
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

#function calculating the number of cluster needed via naive approach
nc_na <- function(type, delta, m, rho00, rho0, sigma0, b4, piz, pix){
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
  lamb <- callamb_na(m, rho00, rho0, sigma0, b4, piz)
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

#function calculating the ratio between sample sizes by Taylor approximation and naive approach
ratio <- function(type, delta, m, rho00, rho0, sigma0, b4, piz, pix){
  nct <- nc(type, delta, m, rho00, rho0, sigma0, b4, piz, pix)
  ncn <- nc_na(type, delta, m, rho00, rho0, sigma0, b4, piz, pix)
  ra <- nct/ncn
  return(ra)
}

piz <- 0.5
pix <- 0.5
m <- c(20,50,100)
rho0 <- c(0.01,0.05)
rho00 <- c(0.02,0.06)
rho11 <- c(0.05,0.10)
b4 <- c()
for (i in 1:length(rho0)){
  b4[i] <- (1-rho0[i])*(rho11[i]-rho00[i])/(1-rho00[i])/(1-rho11[i])
}
delta <- 0.2


#ratio(1, delta1[1], m[1], rho00[1], rho0[1], sigma0[1], b4[1], 0.01, 0.01)
#plot heat map
phm <- function(type,delta,rho0,rho00,b4,m,label){
  
  require(latex2exp)
  x <- seq(0.1,0.9,by=0.005)
  y <- seq(0.1,0.9,by=0.005)
  z <- matrix(data= NA, nrow=length(x),ncol=length(y))
  for (i in 1:length(x)){
    for (j in 1:length(y)){
      z[i,j] <- ratio(type, delta, m, rho00, rho0, 1, b4, y[j], x[i])
    }
  }
  contour(x=x,
          y=y,
          z=z,
          main = TeX(label),
          xlab = TeX("$\\pi_X$"),
          ylab = TeX("$\\pi_Z$"),
          col = "blue",
          lwd = 2,
          cex.main = 2.0,
          cex.lab = 1.3,
          cex = 5)
}
par(mfrow = c(2, 2))
phm(1,delta,rho0[1],rho00[1],b4[1],m[1],"(a): $m$=20, $\\alpha,\\rho_0,\\rho_1=0.01,0.02,0.05$")
phm(1,delta,rho0[1],rho00[1],b4[1],m[2],"(b): $m$=100, $\\alpha,\\rho_0,\\rho_1=0.01,0.02,0.05$")
phm(1,delta,rho0[2],rho00[2],b4[2],m[1],"(c): $m$=20, $\\alpha,\\rho_0,\\rho_1=0.05,0.06,0.10$")
phm(1,delta,rho0[2],rho00[2],b4[2],m[2],"(d): $m$=100, $\\alpha,\\rho_0,\\rho_1=0.05,0.06,0.10$")

phm(3,delta,rho0[1],rho00[1],b4[1],m[1],"(a): $m$=20, $\\alpha,\\rho_0,\\rho_1=0.01,0.02,0.05$")
phm(3,delta,rho0[1],rho00[1],b4[1],m[2],"(b): $m$=100, $\\alpha,\\rho_0,\\rho_1=0.01,0.02,0.05$")
phm(3,delta,rho0[2],rho00[2],b4[2],m[1],"(c): $m$=20, $\\alpha,\\rho_0,\\rho_1=0.05,0.06,0.10$")
phm(3,delta,rho0[2],rho00[2],b4[2],m[2],"(d): $m$=100, $\\alpha,\\rho_0,\\rho_1=0.05,0.06,0.10$")

phm(4,delta,rho0[1],rho00[1],b4[1],m[1],"(a): $m$=20, $\\alpha,\\rho_0,\\rho_1=0.01,0.02,0.05$")
phm(4,delta,rho0[1],rho00[1],b4[1],m[2],"(b): $m$=100, $\\alpha,\\rho_0,\\rho_1=0.01,0.02,0.05$")
phm(4,delta,rho0[2],rho00[2],b4[2],m[1],"(c): $m$=20, $\\alpha,\\rho_0,\\rho_1=0.05,0.06,0.10$")
phm(4,delta,rho0[2],rho00[2],b4[2],m[2],"(d): $m$=100, $\\alpha,\\rho_0,\\rho_1=0.05,0.06,0.10$")

phm(5,delta,rho0[1],rho00[1],b4[1],m[1],"(a): $m$=20, $\\alpha,\\rho_0,\\rho_1=0.01,0.02,0.05$")
phm(5,delta,rho0[1],rho00[1],b4[1],m[2],"(b): $m$=100, $\\alpha,\\rho_0,\\rho_1=0.01,0.02,0.05$")
phm(5,delta,rho0[2],rho00[2],b4[2],m[1],"(c): $m$=20, $\\alpha,\\rho_0,\\rho_1=0.05,0.06,0.10$")
phm(5,delta,rho0[2],rho00[2],b4[2],m[2],"(d): $m$=100, $\\alpha,\\rho_0,\\rho_1=0.05,0.06,0.10$")

m <- c(20,50,100)
rho0 <- c(0.01,0.01,0.01,0.01,0.05,0.05,0.10,0.10)
rho00 <- c(0.02,0.02,0.04,0.04,0.06,0.10,0.12,0.13)
rho11 <- c(0.05,0.12,0.05,0.12,0.08,0.13,0.14,0.15)
b4 <- c()
for (i in 1:length(rho0)){
  b4[i] <- (1-rho0[i])*(rho11[i]-rho00[i])/(1-rho00[i])/(1-rho11[i])
}
delta1 <- c(0.2,0.35)
delta2 <- c(0.15,0.3)
delta3 <- c(0.2,0.35)
delta4 <- c(0.1,0.15)
delta5 <- c(0.2,0.3)

#function generating simulation results
nc_na_whole <- function(type, delta, m, rho00, rho0, sigma0, b4){
  result <- c()
  for (j in 1:length(delta)){
    for (k in 1:length(m)){
      for (l in 1:length(rho00)){
        ncna <- ceiling(nc_na(type, delta[j], m[k], rho00[l], rho0[l], sigma0, b4[l], 0.5, 0.5))
        if (ncna%%2==1){ncna <- ncna+1}
        result <- c(result,ncna)
      }
    }
  }
  return(result)
}


