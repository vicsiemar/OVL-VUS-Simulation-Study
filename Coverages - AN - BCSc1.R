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

#OVL calculado usando kernel gaussiano
ovl_kernel=function(datos1,datos2,datos3,kernel){
  if (kernel=="gauss") {
    kernel=kernel_gauss
  } 
  else {
    stop("Definir nuevo kernel")}
  h1=(4/(3*length(datos1)))^(1/5)*min(sd(datos1),IQR(datos1)/1.349)
  h2=(4/(3*length(datos2)))^(1/5)*min(sd(datos2),IQR(datos2)/1.349)
  h3=(4/(3*length(datos3)))^(1/5)*min(sd(datos3),IQR(datos3)/1.349)
  puntos=seq(min(datos1,datos2,datos3),max(datos1,datos2,datos3),length=1000)
  return(sum(pmin(kernel.densidad(datos1,puntos,h1,kernel),kernel.densidad(datos2,puntos,h2,kernel),kernel.densidad(datos3,puntos,h3,kernel)))*(puntos[2]-puntos[1]))
}

#OVL paramétrico
OVL_parametrico <- function(mu1, sd1, mu2, sd2, mu3, sd3) {
  # Definición de las funciones de densidad f1(x), f2(x) y f3(x)
  f1 <- function(x) {
    return(dnorm(x, mu1, sd1))}
  
  f2 <- function(x) {
    return(dnorm(x, mu2, sd2))}
  
  f3 <- function(x) {
    return(dnorm(x, mu3, sd3))}
  
  # Calcula el mínimo de las tres funciones de densidad para un valor dado de x
  min_density <- Vectorize(function(x) {
    return(min(f1(x), f2(x), f3(x)))
  }, "x")
  
  # Calcula la integral desde menos infinito hasta más infinito del mínimo de las tres funciones
  res_OVL <- integrate(min_density, lower = -Inf, upper = Inf)
  return(res_OVL$value)  
}

#OVL teórico
OVL_teorico <- function(mu1, sd1, mu2, sd2, mu3, sd3) {
  # Definición de las funciones de densidad f1(x), f2(x) y f3(x)
  f1 <- function(x) {
    return(dnorm(x, mu1, sd1))}
  
  f2 <- function(x) {
    return(dnorm(x, mu2, sd2))}
  
  f3 <- function(x) {
    return(dnorm(x, mu3, sd3))}
  
  # Calcula el mínimo de las tres funciones de densidad para un valor dado de x
  min_density <- Vectorize(function(x) {
    return(min(f1(x), f2(x), f3(x)))
  }, "x")
  
  # Calcula la integral desde menos infinito hasta más infinito del mínimo de las tres funciones
  res_OVL <- integrate(min_density, lower = -Inf, upper = Inf)
  return(res_OVL$value)  
}

likbox=function(h,data,n){
  x <- data[1:n]
  y <- data[(n + 1):(n + n)]
  z <- data[(n + n + 1):(n + n + n)]
  if (abs(h)<1e-5){
    xh<-log(x)
    yh<-log(y)
    zh<-log(z)
  } else {
    xh<-((x^h)-1)/h
    yh<-((y^h)-1)/h
    zh<-((z^h)-1)/h
  }
  oout=-n/2*log(sum((xh-sum(xh)/n)^2)/n)-n/2*log(sum((yh-sum(yh)/n)^2)/n)-
    n/2*log(sum((zh-sum(zh)/n)^2)/n)+(h-1)*(sum(log(x))+sum(log(y))+sum(log(z)))
  return(-oout)
}

IC_AN <- function(x, y, z, B, alpha) {
xo = x
yo = y 
zo = z

  nx = length(xo)
  ny = length(yo)
  nz = length(zo)

h_ini=-0.6  #Parámetro a optimizar, valor inicial
  all_values<-c(x,y,z)
  if (any(all_values<=0)){
    x<-x+abs(min(all_values))+(max(all_values)-min(all_values))/2
    y<-y+abs(min(all_values))+(max(all_values)-min(all_values))/2
    z<-z+abs(min(all_values))+(max(all_values)-min(all_values))/2
  }
  
hhat_BFGS=try({optim(h_ini,likbox,data=c(x,y,z),n=length(x),method="BFGS")$par},silent=TRUE)
    if(class(hhat_BFGS)=="try-error"){
      hhat_LBFGSB= try({optim(h_ini,likbox,data=c(x,y,z),n=length(x),method="L-BFGS-B",lower=-2,upper=2)$par},silent=TRUE)
      if(class(hhat_LBFGSB)=="try-error"){
        hhat_NM= try({optim(h_ini,likbox,data=c(x,y,z),n=length(x),method="Nelder-Mead")$par},silent=TRUE)
        if(class(hhat_NM)=="try-error"){hhat=h_ini}else{hhat=hhat_NM}
      }else{hhat= hhat_LBFGSB}
    }else{hhat=hhat_BFGS}

if (abs(hhat)<1e-5){
    x=log(x)
    y=log(y)
    z=log(z)
  } else {
    x=((x^hhat)-2)/hhat
    y=((y^hhat)-2)/hhat
    z=((z^hhat)-2)/hhat
  }


mu1 <- mean(x); mu2 <- mean(y); mu3 <- mean(z)
sd1 <- sd(x); sd2 <- sd(y); sd3 <- sd(z)
#Estimadores en la muestra original transformada
OVL_B_parametrico_o <- OVL_parametrico(mu1, sd1, mu2, sd2, mu3, sd3)
OVL_B_kernel_o = ovl_kernel(x, y, z, "gauss")

  
  OVL_B_parametrico = numeric(B)
  OVL_B_kernel = numeric(B)
  
  for (b in 1:B) {
    # Bootstrap de las muestras
    xB = sample(xo, nx, replace = TRUE)
    yB = sample(yo, ny, replace = TRUE)
    zB = sample(zo, nz, replace = TRUE)
    
 h_ini=-0.6  #Parámetro a optimizar, valor inicial
  all_values<-c(xB,yB,zB)
  if (any(all_values<=0)){
    xB<-xB+abs(min(all_values))+(max(all_values)-min(all_values))/2
    yB<-yB+abs(min(all_values))+(max(all_values)-min(all_values))/2
    zB<-zB+abs(min(all_values))+(max(all_values)-min(all_values))/2
  }
  
hhat_BFGS=try({optim(h_ini,likbox,data=c(xB,yB,zB),n=length(xB),method="BFGS")$par},silent=TRUE)
    if(class(hhat_BFGS)=="try-error"){
      hhat_LBFGSB= try({optim(h_ini,likbox,data=c(xB,yB,zB),n=length(xB),method="L-BFGS-B",lower=-2,upper=2)$par},silent=TRUE)
      if(class(hhat_LBFGSB)=="try-error"){
        hhat_NM= try({optim(h_ini,likbox,data=c(xB,yB,zB),n=length(xB),method="Nelder-Mead")$par},silent=TRUE)
        if(class(hhat_NM)=="try-error"){hhat=h_ini}else{hhat=hhat_NM}
      }else{hhat= hhat_LBFGSB}
    }else{hhat=hhat_BFGS}

if (abs(hhat)<1e-5){
    x=log(xB)
    y=log(yB)
    z=log(zB)
  } else {
    x=((xB^hhat)-2)/hhat
    y=((yB^hhat)-2)/hhat
    z=((zB^hhat)-2)/hhat
  }

    # Paramétrico
    mu1 <- mean(x); mu2 <- mean(y); mu3 <- mean(z)
    sd1 <- sd(x); sd2 <- sd(y); sd3 <- sd(z)

    OVL_B_parametrico[b] <- OVL_parametrico(mu1, sd1, mu2, sd2, mu3, sd3)
    
    # Kernel
    OVL_B_kernel[b] = ovl_kernel(x, y, z, "gauss")
  }
  
  # Intervalos de confianza percentiles
  IC_parametrico = c(OVL_B_parametrico_o-qnorm(1 - alpha / 2)*sd(OVL_B_parametrico),OVL_B_parametrico_o+qnorm(1 - alpha / 2)*sd(OVL_B_parametrico))
  IC_kernel = c(OVL_B_kernel_o-qnorm(1 - alpha / 2)*sd(OVL_B_kernel),OVL_B_kernel_o+qnorm(1 - alpha / 2)*sd(OVL_B_kernel))
  
  return(list(IC_parametrico = IC_parametrico, IC_kernel = IC_kernel))
}

nx=c(20,50,100)
ny=c(20,50,100)
nz=c(20,50,100)

alpha=0.05


##1. Las tres clases siguen una distribución Normal
k=1
MC=1000
alpha=0.05
contador_param=0
contador_ker=0

set.seed(1)

OVL_teorico=OVL_teorico(0,1,1/2,1,1,1)
  
for (t in 1:MC){
  x=rnorm(nx[k],0,1)
  y=rnorm(ny[k],1/2,1)
  z=rnorm(nz[k],1,1)
  
  #Para calcular la cobertura
  resultados_IC=IC_AN(x,y,z,B=500,alpha=0.05)
  contador_param=contador_param + 
    (resultados_IC$IC_parametrico[1] <= OVL_teorico & OVL_teorico <= resultados_IC$IC_parametrico[2])
  contador_ker=contador_ker + 
    (resultados_IC$IC_kernel[1] <= OVL_teorico & OVL_teorico <= resultados_IC$IC_kernel[2])
  
  }

cat("Cobertura OVL paramétrico:", contador_param/MC, "\n" )
cat("Cobertura OVL kernel:", contador_ker/MC , "\n")
