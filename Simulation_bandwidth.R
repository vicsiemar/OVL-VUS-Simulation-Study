kernel.densidad=function(datos,puntos,h,kernel){
  ndatos=length(datos)
  npuntos=length(puntos)
  matk=kernel((puntos%*%t(rep(1,times=ndatos))-t(datos%*%t(rep(1,times=npuntos))))/h)
  as.vector((matk%*%rep(1,ndatos))/(ndatos*h))
}

kernel.distribucion=function(datos,puntos,h,kernel.df){
  ndatos=length(datos)
  npuntos=length(puntos)
  matk=kernel.df((puntos%*%t(rep(1,times=ndatos))-t(datos%*%t(rep(1,times=npuntos))))/h)
  as.vector((matk%*%rep(1,ndatos))/(ndatos)) 
}

kernel.df_gauss=function(u){ 
  pnorm(u)                   
}

kernel_gauss=function(u){ 
  dnorm(u)
}

#VUS calculado usando kernel gaussiano
vus=function(datos1,datos2,datos3,kernel,S) {
  if (kernel=="gauss") {
    kernel.df=kernel.df_gauss
    kernel=kernel_gauss
  } 
  else {
    stop("Definir nuevo kernel")
  }
  n=c(length(datos1),length(datos2),length(datos3))
  m=length(n)/sum(1/n)
  h1=S*sd(datos1)*m^(-1/5)
  h2=S*sd(datos2)*m^(-1/5)
  h3=S*sd(datos3)*m^(-1/5)
  puntos=seq(min(datos1,datos2,datos3),max(datos1,datos2,datos3),length=1000)
  return(sum(kernel.distribucion(datos1,puntos,h1,kernel.df)*(1-kernel.distribucion(datos3,puntos,h3,kernel.df))*kernel.densidad(datos2,puntos,h2,kernel))*(puntos[2]-puntos[1]))
}

#OVL calculado usando kernel gaussiano
ovl=function(datos1,datos2,datos3,kernel,S){
  if (kernel=="gauss") {
    kernel=kernel_gauss
  } 
  else {
    stop("Definir nuevo kernel")}
  n=c(length(datos1),length(datos2),length(datos3))
  m=length(n)/sum(1/n)
  h1=S*sd(datos1)*m^(-1/5)
  h2=S*sd(datos2)*m^(-1/5)
  h3=S*sd(datos3)*m^(-1/5)
  puntos=seq(min(datos1,datos2,datos3),max(datos1,datos2,datos3),length=1000)
  return(sum(pmin(kernel.densidad(datos1,puntos,h1,kernel),kernel.densidad(datos2,puntos,h2,kernel),kernel.densidad(datos3,puntos,h3,kernel)))*(puntos[2]-puntos[1]))
}

#Las dos funciones anteriores se pueden extender para utilizar otros kernels distintos.

nx=c(20,20,20,30,50,50,50,100)
ny=c(20,20,30,50,50,50,100,100)
nz=c(20,30,50,50,50,100,100,100)
S=c(1/2,1,3,6,9)

alpha=0.05

##1. Las tres clases siguen una distribución Normal
k=1
S_val=5
set.seed(1)

MC=1000
pvalor_VUS=numeric(MC)
pvalor_OVL=numeric(MC)
for(t in 1:MC){
  x=rnorm(nx[k],0,1)
  y=rnorm(ny[k],0,1.5)
  z=rnorm(nz[k],0,2)
  VUS_obs=vus(x,y,z,"gauss",S[S_val])
  OVL_obs=ovl(x,y,z,"gauss",S[S_val])

  combined=c(x,y,z)

  B=500
  VUS_B=numeric(B)
  OVL_B=numeric(B)
  for (b in 1:B){
    #Bajo H0 las tres muestras son intercambiables
    combinedB=sample(combined,length(combined),replace=TRUE)
    xnuevo=combinedB[1:nx[k]]
    ynuevo=combinedB[(nx[k]+1):(nx[k]+ny[k])]
    znuevo=combinedB[(nx[k]+ny[k]+1):(nx[k]+ny[k]+nz[k])]

    #Calculamos el VUS de estas muestras
    VUS_B[b]=vus(xnuevo,ynuevo,znuevo,"gauss",S[S_val])
    OVL_B[b]=ovl(xnuevo,ynuevo,znuevo,"gauss",S[S_val])
    }
  pvalor_VUS[t]=mean(VUS_B>VUS_obs)
  pvalor_OVL[t]=mean(OVL_B<OVL_obs)
}
mean(pvalor_OVL<alpha)
mean(pvalor_VUS<alpha)


## 2.Las tres clases siguen una distribución Log-Normal
k=1
S_val=1
set.seed(1)

MC=1000
pvalor_VUS=numeric(MC)
pvalor_OVL=numeric(MC)
for(t in 1:MC){
  x=rlnorm(nx[k],1,0.5)
  y=rlnorm(ny[k],1,1)
  z=rlnorm(nz[k],1,1.5)
  VUS_obs=vus(x,y,z,"gauss",S[S_val])
  OVL_obs=ovl(x,y,z,"gauss",S[S_val])
  
  combined=c(x,y,z)
  
  
  B=500
  VUS_B=numeric(B)
  OVL_B=numeric(B)
  for (b in 1:B){
    #Bajo H0 las tres muestras son intercambiables
    combinedB=sample(combined,length(combined),replace=TRUE)
    xnuevo=combinedB[1:nx[k]]
    ynuevo=combinedB[(nx[k]+1):(nx[k]+ny[k])]
    znuevo=combinedB[(nx[k]+ny[k]+1):(nx[k]+ny[k]+nz[k])]
    
    #Calculamos el VUS de estas muestras
    VUS_B[b]=vus(xnuevo,ynuevo,znuevo,"gauss",S[S_val])
    OVL_B[b]=ovl(xnuevo,ynuevo,znuevo,"gauss",S[S_val])
  }
  pvalor_VUS[t]=mean(VUS_B>VUS_obs)
  pvalor_OVL[t]=mean(OVL_B<OVL_obs)
}
mean(pvalor_OVL<alpha)
mean(pvalor_VUS<alpha)


## 3.Las tres clases siguen una distribución Gamma
k=1
S_val=1
set.seed(1)

MC=1000
pvalor_VUS=numeric(MC)
pvalor_OVL=numeric(MC)
for(t in 1:MC){
  x=rgamma(nx[k],shape=0.2,scale=0.6)
  y=rgamma(ny[k],shape=0.2,scale=0.7)
  z=rgamma(nz[k],shape=0.5,scale=0.5)
  VUS_obs=vus(x,y,z,"gauss",S[S_val])
  OVL_obs=ovl(x,y,z,"gauss",S[S_val])
  
  combined=c(x,y,z)
  
  B=500
  VUS_B=numeric(B)
  OVL_B=numeric(B)
  for (b in 1:B){
    #Bajo H0 las tres muestras son intercambiables
    combinedB=sample(combined,length(combined),replace=TRUE)
    xnuevo=combinedB[1:nx[k]]
    ynuevo=combinedB[(nx[k]+1):(nx[k]+ny[k])]
    znuevo=combinedB[(nx[k]+ny[k]+1):(nx[k]+ny[k]+nz[k])]
    
    #Calculamos el VUS de estas muestras
    VUS_B[b]=vus(xnuevo,ynuevo,znuevo,"gauss",S[S_val])
    OVL_B[b]=ovl(xnuevo,ynuevo,znuevo,"gauss",S[S_val])
  }
  pvalor_VUS[t]=mean(VUS_B>VUS_obs)
  pvalor_OVL[t]=mean(OVL_B<OVL_obs)
}
mean(pvalor_OVL<alpha)
mean(pvalor_VUS<alpha)

## 4.Las tres clases siguen distribuciones diferentes
k=7
S_val=1
set.seed(1)

MC=1000
pvalor_VUS=numeric(MC)
pvalor_OVL=numeric(MC)
for(t in 1:MC){
  x=rnorm(nx[k],0,1)
  y=rgamma(ny[k],2,1)
  z=rlnorm(nz[k],0,1)
  VUS_obs=vus(x,y,z,"gauss",S[S_val])
  OVL_obs=ovl(x,y,z,"gauss",S[S_val])
  
  combined=c(x,y,z)
  
  B=500
  VUS_B=numeric(B)
  OVL_B=numeric(B)
  for (b in 1:B){
    #Bajo H0 las tres muestras son intercambiables
    combinedB=sample(combined,length(combined),replace=TRUE)
    xnuevo=combinedB[1:nx[k]]
    ynuevo=combinedB[(nx[k]+1):(nx[k]+ny[k])]
    znuevo=combinedB[(nx[k]+ny[k]+1):(nx[k]+ny[k]+nz[k])]
    
    #Calculamos el VUS de estas muestras
    VUS_B[b]=vus(xnuevo,ynuevo,znuevo,"gauss",S[S_val])
    OVL_B[b]=ovl(xnuevo,ynuevo,znuevo,"gauss",S[S_val])
  }
  pvalor_VUS[t]=mean(VUS_B>VUS_obs)
  pvalor_OVL[t]=mean(OVL_B<OVL_obs)
}
mean(pvalor_OVL<alpha)
mean(pvalor_VUS<alpha)

## 5.1. Las tres clases siguen distribuciones mixtas (normal+normal)
library(distr)

k=6
S_val=1
set.seed(1)

MC=1000
pvalor_VUS=numeric(MC)
pvalor_OVL=numeric(MC)

f1=UnivarMixingDistribution(distr::Norm(0,1), Norm(3,1), mixCoeff=c(0.5,0.5))
f2=UnivarMixingDistribution(distr::Norm(0,1), Norm(3,1), mixCoeff=c(0.5,0.5))
f3=UnivarMixingDistribution(distr::Norm(0,1), Norm(3,1), mixCoeff=c(0.5,0.5))
rf1=slot(f1,"r")
rf2=slot(f2,"r")
rf3=slot(f3,"r")

for(t in 1:MC){
  x=rf1(nx[k])
  y=rf2(ny[k])
  z=rf3(nz[k])
  VUS_obs=vus(x,y,z,"gauss",S[S_val])
  OVL_obs=ovl(x,y,z,"gauss",S[S_val])
  
  combined=c(x,y,z)
  
  B=500
  VUS_B=numeric(B)
  OVL_B=numeric(B)
  for (b in 1:B){
    #Bajo H0 las tres muestras son intercambiables
    combinedB=sample(combined,length(combined),replace=TRUE)
    xnuevo=combinedB[1:nx[k]]
    ynuevo=combinedB[(nx[k]+1):(nx[k]+ny[k])]
    znuevo=combinedB[(nx[k]+ny[k]+1):(nx[k]+ny[k]+nz[k])]
    
    #Calculamos el VUS de estas muestras
    VUS_B[b]=vus(xnuevo,ynuevo,znuevo,"gauss",S[S_val])
    OVL_B[b]=ovl(xnuevo,ynuevo,znuevo,"gauss",S[S_val])
  }
  pvalor_VUS[t]=mean(VUS_B>VUS_obs)
  pvalor_OVL[t]=mean(OVL_B<OVL_obs)
}
mean(pvalor_OVL<alpha)
mean(pvalor_VUS<alpha)

## 5.2. Las tres clases siguen distribuciones mixtas (gamma+gamma)

library(distr)

k=6
S_val=1
set.seed(1)

MC=1000
pvalor_VUS=numeric(MC)
pvalor_OVL=numeric(MC)

f1=UnivarMixingDistribution(distr::Gammad(1,1), Gammad(4,1), mixCoeff=c(0.5,0.5))
f2=UnivarMixingDistribution(distr::Gammad(1,1), Gammad(4,1), mixCoeff=c(0.5,0.5))
f3=UnivarMixingDistribution(distr::Gammad(1,1), Gammad(4,1), mixCoeff=c(0.5,0.5))
rf1=slot(f1,"r")
rf2=slot(f2,"r")
rf3=slot(f3,"r")

for(t in 1:MC){
  x=rf1(nx[k])
  y=rf2(ny[k])
  z=rf3(nz[k])
  VUS_obs=vus(x,y,z,"gauss",S[S_val])
  OVL_obs=ovl(x,y,z,"gauss",S[S_val])
  
  combined=c(x,y,z)
  
  B=500
  VUS_B=numeric(B)
  OVL_B=numeric(B)
  for (b in 1:B){
    #Bajo H0 las tres muestras son intercambiables
    combinedB=sample(combined,length(combined),replace=TRUE)
    xnuevo=combinedB[1:nx[k]]
    ynuevo=combinedB[(nx[k]+1):(nx[k]+ny[k])]
    znuevo=combinedB[(nx[k]+ny[k]+1):(nx[k]+ny[k]+nz[k])]
    
    #Calculamos el VUS de estas muestras
    VUS_B[b]=vus(xnuevo,ynuevo,znuevo,"gauss",S[S_val])
    OVL_B[b]=ovl(xnuevo,ynuevo,znuevo,"gauss",S[S_val])
  }
  pvalor_VUS[t]=mean(VUS_B>VUS_obs)
  pvalor_OVL[t]=mean(OVL_B<OVL_obs)
}
mean(pvalor_OVL<alpha)
mean(pvalor_VUS<alpha)

## 5.3. Las tres clases siguen distribuciones mixtas (normal+gamma)
library(distr)

k=7
S_val=1
set.seed(1)

MC=1000
pvalor_VUS=numeric(MC)
pvalor_OVL=numeric(MC)

f1=UnivarMixingDistribution(distr::Norm(0,1), Gammad(4,1), mixCoeff=c(0.5,0.5))
f2=UnivarMixingDistribution(distr::Norm(0,1), Gammad(4,1), mixCoeff=c(0.5,0.5))
f3=UnivarMixingDistribution(distr::Norm(0,1), Gammad(4,1), mixCoeff=c(0.5,0.5))
rf1=slot(f1,"r")
rf2=slot(f2,"r")
rf3=slot(f3,"r")

for(t in 1:MC){
  x=rf1(nx[k])
  y=rf2(ny[k])
  z=rf3(nz[k])
  VUS_obs=vus(x,y,z,"gauss",S[S_val])
  OVL_obs=ovl(x,y,z,"gauss",S[S_val])
  
  combined=c(x,y,z)
  
  B=500
  VUS_B=numeric(B)
  OVL_B=numeric(B)
  for (b in 1:B){
    #Bajo H0 las tres muestras son intercambiables
    combinedB=sample(combined,length(combined),replace=TRUE)
    xnuevo=combinedB[1:nx[k]]
    ynuevo=combinedB[(nx[k]+1):(nx[k]+ny[k])]
    znuevo=combinedB[(nx[k]+ny[k]+1):(nx[k]+ny[k]+nz[k])]
    
    #Calculamos el VUS de estas muestras
    VUS_B[b]=vus(xnuevo,ynuevo,znuevo,"gauss",S[S_val])
    OVL_B[b]=ovl(xnuevo,ynuevo,znuevo,"gauss",S[S_val])
  }
  pvalor_VUS[t]=mean(VUS_B>VUS_obs)
  pvalor_OVL[t]=mean(OVL_B<OVL_obs)
}
mean(pvalor_OVL<alpha)
mean(pvalor_VUS<alpha)
