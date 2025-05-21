### VUS y OVL
group <- AL$group
AL <- subset(AL, select = -group)
AL <- -AL

kfront <- as.numeric(AL$kfront)
kfront.list <- split(kfront, group)
clase=c("D-","D0","D+")

#VUS paramétrico
vus<-function(mu1, sd1, mu2, sd2, mu3, sd3){
  lim_int_n<-function(p1){1-pnorm(qnorm(p1,mu1,sd1),mu3,sd3)}
  roc_normal_est<-function(p1,p3){pnorm(qnorm(1-p3,mu3,sd3),mu2,sd2)-
      pnorm(qnorm(p1,mu1,sd1),mu2,sd2)}
  return(pracma::integral2(roc_normal_est,0,1,0,lim_int_n)$Q)
}

#VUS calculado usando estimación no paramétrica
F_emp<-function(fd,x){ecdf(fd)(x)}
q_emp<-function(fd,p){as.numeric(quantile(fd,p))}

VUS_NoParametrico<-function(x,y,z){
  lim_int<-function(p1){1-F_emp(z,F_emp(x,q_emp(x,p1)))}
  roc_emp<-function(p1,p3){F_emp(y,q_emp(z,1-p3))-F_emp(y,q_emp(x,p1))}
  return (pracma::integral2(roc_emp,0,1,0,lim_int)$Q)
}

#OVL paramétrico
ovl <- function(mu1, sd1, mu2, sd2, mu3, sd3) {
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

vectores_ordenados <- function(x, y, z) {
  # Calcula las medias y sus índices ordenados
  medias <- c(mean(x), mean(y), mean(z))
  orden <- order(medias)
  
  # Ordena los vectores en función de las medias
  vectores_ordenados <- list(x, y, z)[orden]
  
  return(vectores_ordenados)
}

set.seed(1)
alpha=0.05

x=na.exclude(as.numeric(kfront.list[[1]]))
y=na.exclude(as.numeric(kfront.list[[2]]))
z=na.exclude(as.numeric(kfront.list[[3]]))

orden <- vectores_ordenados(x, y, z)

x_ord <- orden[[1]]
y_ord <- orden[[2]]
z_ord <- orden[[3]]

mu1=mean(x_ord); mu2=mean(y_ord); mu3=mean(z_ord)
sd1=sd(x_ord); sd2=sd(y_ord); sd3=sd(z_ord)

VUS_obs_parametrico=vus(mu1, sd1, mu2, sd2, mu3, sd3)
VUS_obs_emp=VUS_NoParametrico(x_ord,y_ord,z_ord)
OVL_obs_parametrico=ovl(mu1, sd1, mu2, sd2, mu3, sd3)

B=10000
VUS_B_parametrico=numeric(B)
VUS_B_emp=numeric(B)
OVL_B_parametrico=numeric(B)
for (b in 1:B){
  
  xnuevo=sample(x,length(x),replace=TRUE)
  ynuevo=sample(y,length(y),replace=TRUE)
  znuevo=sample(z,length(z),replace=TRUE)
  
  orden <- vectores_ordenados(xnuevo, ynuevo, znuevo)
  
  xnuevo_ord <- orden[[1]]
  ynuevo_ord <- orden[[2]]
  znuevo_ord <- orden[[3]]
  
  mu1=mean(xnuevo_ord); mu2=mean(ynuevo_ord); mu3=mean(znuevo_ord)
  sd1=sd(xnuevo_ord); sd2=sd(ynuevo_ord); sd3=sd(znuevo_ord)
  
  #Calculamos el VUS de estas muestras
  VUS_B_parametrico[b]=vus(mu1, sd1, mu2, sd2, mu3, sd3)
  VUS_B_emp[b]=VUS_NoParametrico(xnuevo_ord,ynuevo_ord,znuevo_ord)
  OVL_B_parametrico[b]=ovl(mu1, sd1, mu2, sd2, mu3, sd3)
}

par(mfrow = c(1, 3))

#IC percentil y estadísticos para OVL_B (paramétrico)
hist(OVL_B_parametrico,freq=FALSE, main = "OVL (paramétrico)", col = "#E1DEFC")
IC_percentil_OVL_B_parametrico=quantile(OVL_B_parametrico,c(alpha/2,1-alpha/2))
points(IC_percentil_OVL_B_parametrico[1],0,col="red",pch=19)
points(IC_percentil_OVL_B_parametrico[2],0,col="red",pch=19)
abline(v=IC_percentil_OVL_B_parametrico[1],col="red")
abline(v=IC_percentil_OVL_B_parametrico[2],col="red")

mean_OVL_B_parametrico <- mean(OVL_B_parametrico)
sd_OVL_B_parametrico <- sd(OVL_B_parametrico)
min_OVL_B_parametrico <- min(OVL_B_parametrico)
max_OVL_B_parametrico <- max(OVL_B_parametrico)

#IC percentil y estadísticos para VUS_B (paramétrico)
hist(VUS_B_parametrico,freq=FALSE, main = "VUS (paramétrico)", col = "#E1DEFC")
IC_percentil_VUS_B_parametrico=quantile(VUS_B_parametrico,c(alpha/2,1-alpha/2))
points(IC_percentil_VUS_B_parametrico[1],0,col="red",pch=19)
points(IC_percentil_VUS_B_parametrico[2],0,col="red",pch=19)
abline(v=IC_percentil_VUS_B_parametrico[1],col="red")
abline(v=IC_percentil_VUS_B_parametrico[2],col="red")

mean_VUS_B_parametrico <- mean(VUS_B_parametrico)
sd_VUS_B_parametrico <- sd(VUS_B_parametrico)
min_VUS_B_parametrico <- min(VUS_B_parametrico)
max_VUS_B_parametrico <- max(VUS_B_parametrico)

#IC percentil y estadísticos para VUS_B (no paramétrico)
hist(VUS_B_emp,freq=FALSE, main = "VUS (no paramétrico)", col = "#E1DEFC")
IC_percentil_VUS_B_emp=quantile(VUS_B_emp,c(alpha/2,1-alpha/2))
points(IC_percentil_VUS_B_emp[1],0,col="red",pch=19)
points(IC_percentil_VUS_B_emp[2],0,col="red",pch=19)
abline(v=IC_percentil_VUS_B_emp[1],col="red")
abline(v=IC_percentil_VUS_B_emp[2],col="red")

mean_VUS_B_emp <- mean(VUS_B_emp)
sd_VUS_B_emp <- sd(VUS_B_emp)
min_VUS_B_emp <- min(VUS_B_emp)
max_VUS_B_emp <- max(VUS_B_emp)

# Resultados

cat("Estimación puntual del OVL (paramétrico):", OVL_obs_parametrico, "\n")
cat("Intervalo de confianza percentil para OVL_B (paramétrico):", IC_percentil_OVL_B_parametrico, "\n")
cat("OVL_B - Media (paramétrico):", mean_OVL_B_parametrico, "\n")
cat("OVL_B - Desviación Típica (paramétrico):", sd_OVL_B_parametrico, "\n")
cat("OVL_B - Mínimo (paramétrico):", min_OVL_B_parametrico, "\n")
cat("OVL_B - Máximo (paramétrico):", max_OVL_B_parametrico, "\n")

cat("Estimación puntual del VUS (paramétrico):", VUS_obs_parametrico, "\n")
cat("Intervalo de confianza percentil para VUS_B (paramétrico):", IC_percentil_VUS_B_parametrico, "\n")
cat("VUS_B - Media (paramétrico):", mean_VUS_B_parametrico, "\n")
cat("VUS_B - Desviación Típica (paramétrico):", sd_VUS_B_parametrico, "\n")
cat("VUS_B - Mínimo (paramétrico):", min_VUS_B_parametrico, "\n")
cat("VUS_B - Máximo (paramétrico):", max_VUS_B_parametrico, "\n")

cat("Estimación puntual del VUS (no paramétrico):", VUS_obs_emp, "\n")
cat("Intervalo de confianza percentil para VUS_B (no paramétrico):", IC_percentil_VUS_B_emp, "\n")
cat("VUS_B - Media (no paramétrico):", mean_VUS_B_emp, "\n")
cat("VUS_B - Desviación Típica (no paramétrico):", sd_VUS_B_emp, "\n")
cat("VUS_B - Mínimo (no paramétrico):", min_VUS_B_emp, "\n")
cat("VUS_B - Máximo (no paramétrico):", max_VUS_B_emp, "\n")

