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

#callamb(50,0.2,0.05,1,0.05,0.5)
#for (i in 1:5){
#  print(nc(i,0.25,50,0.2,0.1,1,0.05,0.5,0.5))
#}
#approxp(100,0.04,0.02,1,0.0107,0.5)

#Simulation Codes:

#Power calculation function 
power_cal <- function(type, nc, delta, m, rho00, rho0, sigma0, b4, piz, pix){
  
  #Input
  #type: test type: 1=cex 2=cez 3=mex 4=mez 5=ie
  #nc: number of clusters
  #delta: effect size
  #m: cluster size
  #rho00: rho_0,0
  #rho0: rho_0
  #sigma0: sigma_0^2
  #b4: sigma_b_4^2
  #piz: pi_Z
  #pix: pi_X
  
  
  #output:
  #power
  
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
  var <- vp1+vp2
  power <- pnorm(sqrt(nc*delta^2/var) - qnorm(1-alpha/2))
  return(power)
}

for (i in 1:5){
  print(power_cal(i,100,0.25,50,0.2,0.1,1,0.05,0.5,0.5))
}

#data-generation function
data_gene <- function(type, nc, delta, m, rho00, rho0, sigma0, b4){
  
  #Input
  #type: test type: 1=cex 2=cez 3=mex 4=mez 5=ie
  #nc: number of clusters
  #delta: effect size
  #m: cluster size
  #rho00: rho_0,0
  #rho0: rho_0
  #sigma0: sigma_0^2
  #b4: sigma_b_4^2
  
  
  #output:
  #sim_data: simulated data set
  
  #pi_X
  pix <- 0.5
  #pi_Z
  piz <- 0.5
  
  #effect size
  beta1 = 1
  if (type==1){
    beta2 <- delta
    beta3 <- 0.05
    beta4 <- 0.05
  } else if (type==2){
    beta2 <- 0.15
    beta3 <- delta
    beta4 <- 0.05
  } else if (type==3){
    beta2 <- 0.15
    beta3 <- 0.05
    beta4 <- 2*delta-0.3
  } else if (type==4){
    beta2 <- 0.15
    beta3 <- 0.05
    beta4 <- 2*delta-0.1
  } else if (type==5){
    beta2 <- 0.15
    beta3 <- 0.05
    beta4 <- delta
  }
  #cluster size
  cs <- rep(m, nc)
  #generate Xi
  ia <- rep(0, nc)
  rs <- sample(1:nc, size = nc*pix)
  ia[rs] = 1
  X <- rep(ia, cs)
  #generate Zij
  Z <- rbinom(sum(cs),1,piz)
  #generate outcome
  var_epsilon <- sigma0*(1-rho0)
  var_gamma <- sigma0-var_epsilon
  b2 <- sigma0*(rho00-rho0)/(1-rho00)   #sigma_b_2^2
  gamma <- rep(rnorm(nc,0,sqrt(var_gamma)) ,cs)
  epsilon <- rnorm(sum(cs),0,sqrt(var_epsilon))
  sb2 <- rep(rnorm(nc,0,sqrt(b2)) ,cs)
  sb4 <- rep(rnorm(nc,0,sqrt(b4)) ,cs)
  #alternative hypothesis
  Y_alt <- beta1 + (beta2+sb2)*X + beta3*Z + (beta4+sb4)*X*Z + gamma + epsilon
  #null hypothesis
  if (type==1){
    Y_null <-  beta1 + sb2*X + beta3*Z + (beta4+sb4)*X*Z + gamma + epsilon
  } else if (type==2){
    Y_null <- beta1 + (beta2+sb2)*X + (beta4+sb4)*X*Z + gamma + epsilon
  } else if (type==3){
    Y_null <- beta1 + (beta2+sb2)*X + beta3*Z + (-0.3+sb4)*X*Z + gamma + epsilon
  } else if (type==4){
    Y_null <- beta1 + (beta2+sb2)*X + beta3*Z + (-0.1+sb4)*X*Z + gamma + epsilon
  } else if (type==5){
    Y_null <- beta1 + (beta2+sb2)*X + beta3*Z + sb4*X*Z + gamma + epsilon
  }
  
  #whole index
  ind <- rep(1:nc, cs)
  
  sim_data <- data.frame(ind, X, Z, Y_alt, Y_null)
  
  return(sim_data)
}

#test1 <- data_gene(1,60,0.2,20,0.04,0.02,1,0.01074561)

#function implementing simulation
simula <- function(type, delta, m, rho00, rho0, sigma0, b4, B = 5000){
  
  #Input
  #type: test type: 1=cex 2=cez 3=mex 4=mez 5=ie
  #delta: effect size
  #m: cluster size
  #rho00: rho_0,0
  #rho0: rho_0
  #sigma0: sigma_0^2
  #b4: sigma_b_4^2
  #B: simulation times
  
  require(nlme)
  #type i error rate
  alpha <- 0.05
  #1-power
  mu <- 0.2
  #pi_X
  pix <- 0.5
  #pi_Z
  piz <- 0.5
  
  #calculate the sample size via our formula
  mb <- m*piz
  nc <- ceiling(nc(type, delta, m, rho00, rho0, sigma0, b4, piz, pix))
  #round up to the nearest even integer
  if (nc %% 2 ==1){nc <- nc + 1}
  #store p-values for null distribution
  #pvn <- array(NA,dim=B)
  pvn <- array(NA,dim=B)
  #store p-values for alternative distribution
  pva <- array(NA,dim=B)
  #position
  po2 <- 4
  if (type==1){
    po1 <- 2
    co <- 0
  } else if (type==2){
    po1 <- 3
    co <- 0
  } else if (type==3){
    po1 <- 2
    co <- piz
  } else if (type==4){
    po1 <- 3
    co <- pix
  } else if (type==5){
    po1 <- 4
    co <- 0
  }
  for (i in 1:B){
    
    #generate data
    #indicator for singular-fitting
    exit <- FALSE
    while (exit==FALSE){
      df <- data_gene(type, nc, delta, m, rho00, rho0, sigma0, b4)
      lmm1 <- try(lme(Y_null ~  X*Z, data = df, 
                      random = list(ind = pdDiag(~ X*Z-Z))))
      lmm2 <- try(lme(Y_alt ~  X*Z, data = df, 
                      random = list(ind = pdDiag(~ X*Z-Z))))
      if(class(lmm1)!="try-error"&class(lmm2)!="try-error"){
        exit <- TRUE
      }
    }
    #t-test
    pvn[i] <- coef(summary(lmm1))[po1,5]
    pva[i] <- coef(summary(lmm2))[po1,5]
  }
  #calculate true power
  tp <- power_cal(type, nc, delta, m, rho00, rho0, sigma0, b4, piz, pix)
  output <- as.data.frame(cbind(type, nc,round(mean(pvn[]<0.05, na.rm=T),3),round(tp,3),round(mean(pva[]<0.05, na.rm=T),3)))
  names(output) <- c("test.type","number.of.cluster","empirical.typeI.error","predicted.power","empirical.power")
  return(output)
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
library(nlme)
#simula(1, delta1[1], m[1], rho00[2], rho0[2], 1, b4[2], B = 1000)

#function generating simulation results
simre <- function(type, delta, m, rho00, rho0, sigma0, b4, B = 5000){
  result <- c()
  for (j in 1:length(delta)){
    for (k in 1:length(m)){
      for (l in 1:length(rho00)){
        result <- rbind(result,simula(type, delta[j], m[k], rho00[l], rho0[l], sigma0, b4[l], B = 5000))
      }
    }
  }
  return(result)
}

#function adding label
addla <- function(type,rho11,rho0,rho00,m,delta,resul){
  resul$type <- rep(type,48)
  resul$rho11 <- rep(rho11,6)
  resul$rho0 <- rep(rho0,6)
  resul$rho00 <- rep(rho00,6)
  resul$m <- rep(rep(m,c(8,8,8)),2)
  resul$delta <- rep(delta,c(24,24))
  return(resul)
}

rety1 <- simre(1,delta1,m,rho00,rho0,1,b4,B=5000)
finalty1 <- addla(1,rho11,rho0,rho00,m,delta1,rety1)
write.csv(finalty1,"/Users/deckard/desktop/Fan Li Project/Project 2/t-test simulation results for CEX.csv")

rety2 <- simre(2,delta2,m,rho00,rho0,1,b4,B=5000)
finalty2 <- addla(2,rho11,rho0,rho00,m,delta2,rety2)
write.csv(finalty2,"/Users/deckard/desktop/Fan Li Project/Project 2/t-test simulation results for CEZ.csv")

rety5 <- simre(5,delta5,m,rho00,rho0,1,b4,B=5000)
finalty5 <- addla(5,rho11,rho0,rho00,m,delta5,rety5)
write.csv(finalty5,"/Users/deckard/desktop/Fan Li Project/Project 2/t-test simulation results for IE.csv")