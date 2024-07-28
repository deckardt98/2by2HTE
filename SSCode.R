#####################Sample size calculation####################################

#Eigenvalue of the correlation matrix for the basic model, i.e., lambda(m) 
f <- function(icc,m){
  #Input:
  #icc: alpha or rho0
  #m: cluster size
  
  #Output: eigenvalue of the correlation matrix
  return((icc*(m-1)+1))
}

#function calculating approximations for p0-p3
approxp <- function(m, rho0, alpha, sigma0, b4, piz){
  #Input
  #m: cluster size
  #rho0
  #alpha
  #sigma0: sigma_0^2
  #b4: sigma_{b_4}^2
  #piz: pi_Z
  
  #Output:
  #app: a vector of (p0,p1,p2,p3)
  
  #sigma_epsilon^2
  ep <- sigma0*(1-alpha)
  #mean of m1i
  mb <- m*piz
  #sigma_Z^2
  sZ <- piz*(1-piz)
  
  #determinant
  D <- f(rho0,m)^2*b4^2+4*b4*rho0*ep*f(rho0,m)
  #roots
  x1 <- (f(rho0,m)*b4-sqrt(D))/(2*b4*rho0)
  x2 <- (f(rho0,m)*b4+sqrt(D))/(2*b4*rho0)
  #factors
  f1 <- ((x1+x2)*(x1+x2-3*mb)-x1*x2+3*mb^2)/((x1-mb)^3*(x2-mb)^3)
  f2 <- (x1+x2-2*mb)/((x1-mb)^2*(x2-mb)^2)-2*mb*f1
  f3 <- 3*mb^2*f1+(x1*x2-mb^2)/((x1-mb)^2*(x2-mb)^2)-
    2*mb*(x1*x2*(x1+x2-3*mb)+mb^3)/((x1-mb)^3*(x2-mb)^3)
  #approximated p0
  app0 <- -1/(b4*rho0)*(f1*(m*sZ+mb^2)+f2*mb+f3)
  #third moment for m1i
  tm <- m*(m-1)*(m-2)*piz^3+3*m*(m-1)*piz^2+m*piz
  #approximated p1
  app1 <- -1/(b4*rho0)*(f1*tm+f2*(m*sZ+mb^2)+f3*mb)
  #approximated p2
  app2 <- -1/(b4*rho0)+f(rho0,m)/rho0*app1+ep*f(rho0,m)/(b4*rho0)*app0
  #approximated p3
  app3 <- -mb/(b4*rho0)+f(rho0,m)/rho0*app2+ep*f(rho0,m)/(b4*rho0)*app1
  
  #output
  app <- c(app0,app1,app2,app3)
  
  return(app)
}

#Function calculating lambda11 lambda12 lambda22
callamb <- function(m, rho0, alpha, sigma0, b4, piz){
  #Input
  #m: cluster size
  #rho0 
  #alpha
  #sigma0: sigma_0^@
  #b4: sigma_b_4^2
  #piz: pi_Z
  
  #output:
  #estlamb lambda11,lambda12,lambda22
  
  #sigma_epsilon^2
  ep <- sigma0*(1-alpha)
  esp <- approxp(m, rho0, alpha, sigma0, b4, piz)
  p0 <- esp[1]
  p1 <- esp[2]
  p2 <- esp[3]
  p3 <- esp[4]
  
  #calculate lambda11 lambda12 lambda22
  lamb11 <- (m+b4*(2*m*rho0-f(rho0,m))*p2-rho0*m^2*b4*p1-rho0*m^2*ep*p0)/(ep)
  lamb12 <- (m*piz+rho0*b4*p3-f(rho0,m)*b4*p2-m*rho0*ep*p1)/(ep)
  lamb22 <- (m*piz+rho0*b4*p3-(f(rho0,m)*b4+rho0*ep)*p2)/(ep)
  
  #output:
  estlamb <- c(lamb11,lamb12,lamb22)
  return(estlamb)
}

#Function calculating the number of cluster  
nc <- function(type, delta, m, rho0, alpha, sigma0, b4, piz, pix){
  #Input
  #type: test type 1=cex 2=cez 3=mex 4=mez 5=ie
  #delta: effect size
  #m: cluster size
  #rho0
  #alpha
  #sigma0: sigma_0^2
  #b4: sigma_b_4^2
  #piz: pi_Z
  #pix: pi_X
  
  #Output:
  #n: the number of cluster 
  
  #type I error rate
  epsilon1 <- 0.05
  #power
  epsilon2 <- 0.2
  #sigma_Z^2
  sz <- piz*(1-piz)
  #calculate lambda
  lamb <- callamb(m, rho0, alpha, sigma0, b4, piz)
  l11 <- lamb[1]
  l12 <- lamb[2]
  l22 <- lamb[3]
  if (type==1){
    #variance part 1
    vp1 <- (sigma0*((m-1)*alpha+1)*((m-1)*alpha*(1-piz)+1-alpha))/(m*(1-pix)*(1-piz)*((m-2)*alpha+1))
    #variance part2
    vp2 <- l22/(pix*(l11*l22-l12^2))
  }else if (type==2){
    #variance part 1
    vp1 <- (sigma0*((m-1)*alpha+1)*(1-alpha))/(m*(1-pix)*sz*((m-2)*alpha+1))
    #variance part2
    vp2 <- 0
  }else if (type==3){
    #variance part 1
    vp1 <- (sigma0*((m-1)*alpha+1))/(m*(1-pix))
    #variance part2
    vp2 <- (piz^2*l11-2*piz*l12+l22)/(pix*(l11*l22-l12^2))
  } else if (type==4){
    #variance part 1
    vp1 <- (sigma0*(1-alpha)*((m-1)*alpha+1)*(1-pix))/(m*sz*((m-2)*alpha+1))
    #variance part2
    vp2 <- (pix*l11)/(l11*l22-l12^2)
  }else if (type==5){
    #variance part 1
    vp1 <- (sigma0*(1-alpha)*((m-1)*alpha+1))/(m*(1-pix)*sz*((m-2)*alpha+1))
    #variance part2
    vp2 <- l11/(pix*(l11*l22-l12^2))
  }
  n <- (qnorm(1-epsilon1/2)+qnorm(1-epsilon2))^2/delta^2*(vp1+vp2)
  return(n)
}


#Function calculating p0-p3 using naive approach 
approxp_na <- function(m, rho0, alpha, sigma0, b4, piz){
  #Input
  #m: cluster size
  #rho0
  #alpha
  #sigma0: sigma_0^@
  #b4: sigma_b_4^2
  #piz: pi_Z
  
  #Output:
  #app: a vector of (p0,p1,p2,p3)
  
  #sigma_epsilon^2
  ep <- sigma0*(1-alpha)
  #mean of m1i
  mb <- m*piz
  #sigma_Z^2
  sZ <- piz*(1-piz)
  
  #approximated p0
  app0 <- 1/(-b4*rho0*mb^2+f(rho0,m)*b4*mb+ep*f(rho0,m))
  #third moment for m1i
  tm <- m*(m-1)*(m-2)*piz^3+3*m*(m-1)*piz^2+m*piz
  #approximated p1
  app1 <- mb*app0
  #approximated p2
  app2 <- -1/(b4*rho0)+f(rho0,m)/rho0*app1+ep*f(rho0,m)/(b4*rho0)*app0
  #approximated p3
  app3 <- -mb/(b4*rho0)+f(rho0,m)/rho0*app2+ep*f(rho0,m)/(b4*rho0)*app1
  
  #output
  app <- c(app0,app1,app2,app3)
  return(app)
}

#Function calculating lambda11 lambda12 lambda22 via naive approach
callamb_na <- function(m, rho0, alpha, sigma0, b4, piz){
  #Input
  #m: cluster size
  #rho0
  #alpha
  #sigma0: sigma_0^@
  #b4: sigma_b_4^2
  #piz: pi_Z
  
  #output:
  #estlamb lambda11,lambda12,lambda22
  
  #sigma_epsilon^2
  ep <- sigma0*(1-alpha)
  esp <- approxp_na(m, rho0, alpha, sigma0, b4, piz)
  p0 <- esp[1]
  p1 <- esp[2]
  p2 <- esp[3]
  p3 <- esp[4]
  
  #calculate lambda11 lambda12 lambda22
  lamb11 <- (m+b4*(2*m*rho0-f(rho0,m))*p2-rho0*m^2*b4*p1-rho0*m^2*ep*p0)/(ep)
  lamb12 <- (m*piz+rho0*b4*p3-f(rho0,m)*b4*p2-m*rho0*ep*p1)/(ep)
  lamb22 <- (m*piz+rho0*b4*p3-(f(rho0,m)*b4+rho0*ep)*p2)/(ep)
  
  #output
  estlamb <- c(lamb11,lamb12,lamb22)
  return(estlamb)
}

#Function calculating the number of cluster via naive approach
nc_na <- function(type, delta, m, rho0, alpha, sigma0, b4, piz, pix){
  #Input
  #type: test type 1=cex 2=cez 3=mex 4=mez 5=ie
  #delta: effect size
  #m: cluster size
  #rho0
  #alpha
  #sigma0: sigma_0^2
  #b4: sigma_b_4^2
  #piz: pi_Z
  #pix: pi_X
  
  #Output:
  #n: the number of cluster 
  
  #type i error rate
  epsilon1 <- 0.05
  #power
  epsilon2 <- 0.2
  #sigma_Z^2
  sz <- piz*(1-piz)
  #calculate lambda
  lamb <- callamb_na(m, rho0, alpha, sigma0, b4, piz)
  l11 <- lamb[1]
  l12 <- lamb[2]
  l22 <- lamb[3]
  if (type==1){
    #variance part 1
    vp1 <- (sigma0*((m-1)*alpha+1)*((m-1)*alpha*(1-piz)+1-alpha))/(m*(1-pix)*(1-piz)*((m-2)*alpha+1))
    #variance part2
    vp2 <- l22/(pix*(l11*l22-l12^2))
  }else if (type==2){
    #variance part 1
    vp1 <- (sigma0*((m-1)*alpha+1)*(1-alpha))/(m*(1-pix)*sz*((m-2)*alpha+1))
    #variance part2
    vp2 <- 0
  }else if (type==3){
    #variance part 1
    vp1 <- (sigma0*((m-1)*alpha+1))/(m*(1-pix))
    #variance part2
    vp2 <- (piz^2*l11-2*piz*l12+l22)/(pix*(l11*l22-l12^2))
  } else if (type==4){
    #variance part 1
    vp1 <- (sigma0*(1-alpha)*((m-1)*alpha+1)*(1-pix))/(m*sz*((m-2)*alpha+1))
    #variance part2
    vp2 <- (pix*l11)/(l11*l22-l12^2)
  }else if (type==5){
    #variance part 1
    vp1 <- (sigma0*(1-alpha)*((m-1)*alpha+1))/(m*(1-pix)*sz*((m-2)*alpha+1))
    #variance part2
    vp2 <- l11/(pix*(l11*l22-l12^2))
  }
  n <- (qnorm(1-epsilon1/2)+qnorm(1-epsilon2))^2/delta^2*(vp1+vp2)
  return(n)
}
