#f(m)
f <- function(icc,m){
  return((icc*(m-1)+1))
}
#function approximating p0 to p3 via naive approach
approxp_naive <- function(m, rho00, rho0, sigma0, b4, piz){
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
  
  #approximated p0
  app0 <- 1/(ep*f(rho00,m)+mb*b4*f(rho00,(m-mb)))
  #approximated p1
  app1 <- app0*mb
  #approximated p2
  app2 <- -1/(b4*rho00)+f(rho00,m)/rho00*app1+ep*f(rho00,m)/(b4*rho00)*app0
  #approximated p3
  app3 <- -mb/(b4*rho00)+f(rho00,m)/rho00*app2+ep*f(rho00,m)/(b4*rho00)*app1
  #vector output
  app <- c(app0,app1,app2,app3)
  return(app)
}

#function calculating lambda11 lambda12 lambda22 via naive approach
callamb_naive <- function(m, rho00, rho0, sigma0, b4, piz){
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
  esp <- approxp_naive(m, rho00, rho0, sigma0, b4, piz)
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
nc_naive <- function(type, delta, m, rho00, rho0, sigma0, b4, piz, pix){
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
  lamb <- callamb_naive(m, rho00, rho0, sigma0, b4, piz)
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

B <- 5000
m <- c(20,50,100)
rho0 <- c(0.01,0.01,0.01,0.01,0.10,0.10,0.10,0.10)
rho00 <- c(0.02,0.02,0.04,0.04,0.12,0.13,0.12,0.13)
rho11 <- c(0.05,0.12,0.05,0.12,0.14,0.15,0.14,0.15)
b4 <- c()
for (i in 1:length(rho0)){
  b4[i] <- (1-rho0[i])*(rho11[i]-rho00[i])/(1-rho00[i])/(1-rho11[i])
}
delta1 <- c(0.2,0.4)
delta2 <- c(0.15,0.3)
delta3 <- c(0.2,0.4)
delta4 <- c(0.1,0.15)
delta5 <- c(0.2,0.3)
