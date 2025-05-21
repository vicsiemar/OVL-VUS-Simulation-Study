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

IC_AN <- function(x, y, z, B, alpha) {
mu1 <- mean(x); mu2 <- mean(y); mu3 <- mean(z)
sd1 <- sd(x); sd2 <- sd(y); sd3 <- sd(z)
#Estimadores en la muestra original
OVL_B_parametrico_o <- OVL_parametrico(mu1, sd1, mu2, sd2, mu3, sd3)
OVL_B_kernel_o = ovl_kernel(x, y, z, "gauss")


  nx = length(x)
  ny = length(y)
  nz = length(z)
  
  OVL_B_parametrico = numeric(B)
  OVL_B_kernel = numeric(B)
  
  for (b in 1:B) {
    # Bootstrap de las muestras
    xB = sample(x, nx, replace = TRUE)
    yB = sample(y, ny, replace = TRUE)
    zB = sample(z, nz, replace = TRUE)
    
    # Paramétrico
    mu1 <- mean(xB); mu2 <- mean(yB); mu3 <- mean(zB)
    sd1 <- sd(xB); sd2 <- sd(yB); sd3 <- sd(zB)
    OVL_B_parametrico[b] <- OVL_parametrico(mu1, sd1, mu2, sd2, mu3, sd3)
    
    # Kernel
    OVL_B_kernel[b] = ovl_kernel(xB, yB, zB, "gauss")
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
