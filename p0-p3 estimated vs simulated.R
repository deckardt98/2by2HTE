#f(m)
f <- function(icc,m){
  return((icc*(m-1)+1))
}

#function simulate p0-p3
simp <- function(ep, icc, b4, m, piZ, B){
  
  #Input:
  #ep: sigma_ep^2
  #icc: rho_0,0
  #b4: sigma_b_4^2
  #m: cluster size
  #piZ: pi_Z
  #B: simulation time
  
#mean of m1i
mb <- m*piZ
#sigma_Z^2
sZ <- piZ*(1-piZ)

#generating m_1i
mi <- rbinom(B,m,piZ)
#simulate p0
p0 <- 1/(ep*f(icc,m)+mi*b4*f(icc,(m-mi)))
#simulate p1 p2 p3
p1 <- p0*mi
p2 <- p1*mi
p3 <- p2*mi
#simulated p0 p1 p2 p3
sp0 <- mean(p0)
sp1 <- mean(p1)
sp2 <- mean(p2)
sp3 <- mean(p3)

return(c(sp0,sp1,sp2,sp3))
}

appp <- function(ep, icc, b4, m, piZ){
  
  #Input:
  #ep: sigma_ep^2
  #icc: rho_0,0
  #b4: sigma_b_4^2
  #m: cluster size
  #piZ: pi_Z
  
#mean of m1i
mb <- m*piZ
#sigma_Z^2
sZ <- piZ*(1-piZ)
#determinant
D <- f(icc,m)^2*b4^2+4*b4*icc*ep*f(icc,m)
#roots
x1 <- (f(icc,m)*b4-sqrt(D))/(2*b4*icc)
x2 <- (f(icc,m)*b4+sqrt(D))/(2*b4*icc)
#factors
f1 <- ((x1+x2)*(x1+x2-3*mb)-x1*x2+3*mb^2)/((x1-mb)^3*(x2-mb)^3)
f2 <- (x1+x2-2*mb)/((x1-mb)^2*(x2-mb)^2)-2*mb*f1
f3 <- 3*mb^2*f1+(x1*x2-mb^2)/((x1-mb)^2*(x2-mb)^2)-
  2*mb*(x1*x2*(x1+x2-3*mb)+mb^3)/((x1-mb)^3*(x2-mb)^3)
#approximated p0
app0 <- -1/(b4*icc)*(f1*(m*sZ+mb^2)+f2*mb+f3)
#third moment for m1i
tm <- m*(m-1)*(m-2)*piZ^3+3*m*(m-1)*piZ^2+m*piZ
#approximated p1
app1 <- -1/(b4*icc)*(f1*tm+f2*(m*sZ+mb^2)+f3*mb)
#approximated p2
app2 <- -1/(b4*icc)+f(icc,m)/icc*app1+ep*f(icc,m)/(b4*icc)*app0
#approximated p3
app3 <- -mb/(b4*icc)+f(icc,m)/icc*app2+ep*f(icc,m)/(b4*icc)*app1
#app3 <- -mb/b4/icc-f(icc,m)/b4/icc^2+f(icc,m)/icc*(f(icc,m)/icc+ep/b4)*app1+ep*f(icc,m)^2/b4/icc^2*app0
return(c(app0,app1,app2,app3))
}
#test
simp(0.98,0.04,0.0107,100,0.5,100000)
appp(0.98,0.04,0.0107,100,0.5)

