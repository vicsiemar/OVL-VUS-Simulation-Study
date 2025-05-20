library(plotly)
library(distr)

group <- AL$group
AL <- subset(AL, select = -group)
AL <- -AL

kfront <- as.numeric(AL$kfront)
kfront.list <- split(kfront, group)

kernel.distribucion=function(datos,puntos,h,kernel.df){
  ndatos=length(datos)
  npuntos=length(puntos)
  matk=kernel.df((puntos%*%t(rep(1,times=ndatos))-t(datos%*%t(rep(1,times=npuntos))))/h)
  as.vector((matk%*%rep(1,ndatos))/(ndatos)) 
}

kernel.df_gauss <- function(u) {
  pnorm(u)
}

x=na.exclude(as.numeric(kfront.list[[1]]))
y=na.exclude(as.numeric(kfront.list[[2]]))
z=na.exclude(as.numeric(kfront.list[[3]]))

h1=(4/(3*length(x)))^(1/5)*min(sd(x),IQR(x)/1.349)
h2=(4/(3*length(y)))^(1/5)*min(sd(y),IQR(y)/1.349)
h3=(4/(3*length(z)))^(1/5)*min(sd(z),IQR(z)/1.349)
puntos=seq(min(x,y,z),max(x,y,z),length=1000)

q_emp<-function(x,p){as.numeric(quantile(x,p))}

ROC_NoParametrico <- function(p1, p3) {
  result <- ifelse(q_emp(x,p1) <= q_emp(z,1-p3),
                   kernel.distribucion(y, q_emp(z,1-p3), h2, kernel.df_gauss) -
                     kernel.distribucion(y, q_emp(x,p1), h2, kernel.df_gauss),
                   0)
  return(result)
}

p1 <- seq(0, 1, 0.01)
p3 <- p1 

z1 <- outer(1-p3, p1, ROC_NoParametrico)

# Crear la figura
p <- plot_ly(x = p1, y = 1-p3, z = z1) %>% add_surface(colorscale = "Rainbow") %>% 
  layout(
    title = "",
    scene = list(
      xaxis = list(title = list(text = "p1", font = list(size = 16, color = "black", family = "Arial")), range = c(0, 1), tickfont = list(size = 12)),
      yaxis = list(title = list(text = "1-p3", font = list(size = 16, color = "black", family = "Arial")), range = c(0, 1), tickfont = list(size = 12)),
      zaxis = list(title = list(text = "ROC", font = list(size = 16, color = "black", family = "Arial")), range = c(0, 1), tickfont = list(size = 12))
    ),
    margin = list(l = 0, r = 0, b = 0, t = 0)  # Eliminar márgenes
  )

# Mostrar la superficie ROC
hide_colorbar(p)
