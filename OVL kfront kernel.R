# --- FUNCIONES KERNEL ---

kernel.densidad = function(datos, puntos, h, kernel){
  ndatos = length(datos)
  npuntos = length(puntos)
  matk = kernel((puntos %*% t(rep(1, times = ndatos)) - t(datos %*% t(rep(1, times = npuntos)))) / h)
  as.vector((matk %*% rep(1, ndatos)) / (ndatos * h))
}

kernel.distribucion = function(datos, puntos, h, kernel.df){
  ndatos = length(datos)
  npuntos = length(puntos)
  matk = kernel.df((puntos %*% t(rep(1, times = ndatos)) - t(datos %*% t(rep(1, times = npuntos)))) / h)
  as.vector((matk %*% rep(1, ndatos)) / ndatos) 
}

kernel.df_gauss = function(u){ 
  pnorm(u)                   
}

kernel_gauss = function(u){ 
  dnorm(u)
}

# --- FUNCIONES VUS y OVL ---

vus = function(datos1, datos2, datos3, kernel) {
  if (kernel == "gauss") {
    kernel.df = kernel.df_gauss
    kernel = kernel_gauss
  } else {
    stop("Definir nuevo kernel")
  }
  datos1 = na.exclude(as.numeric(datos1))
  datos2 = na.exclude(as.numeric(datos2))
  datos3 = na.exclude(as.numeric(datos3))
  
  h1 = (4 / (3 * length(datos1)))^(1/5) * min(sd(datos1), IQR(datos1, na.rm = TRUE) / 1.349)
  h2 = (4 / (3 * length(datos2)))^(1/5) * min(sd(datos2), IQR(datos2, na.rm = TRUE) / 1.349)
  h3 = (4 / (3 * length(datos3)))^(1/5) * min(sd(datos3), IQR(datos3, na.rm = TRUE) / 1.349)
  
  puntos = seq(min(datos1, datos2, datos3), max(datos1, datos2, datos3), length = 1000)
  return(sum(kernel.distribucion(datos1, puntos, h1, kernel.df) *
               (1 - kernel.distribucion(datos3, puntos, h3, kernel.df)) *
               kernel.densidad(datos2, puntos, h2, kernel)) * (puntos[2] - puntos[1]))
}

ovl = function(datos1, datos2, datos3, kernel){
  if (kernel == "gauss") {
    kernel = kernel_gauss
  } else {
    stop("Definir nuevo kernel")
  }
  datos1 = na.exclude(as.numeric(datos1))
  datos2 = na.exclude(as.numeric(datos2))
  datos3 = na.exclude(as.numeric(datos3))
  
  h1 = (4 / (3 * length(datos1)))^(1/5) * min(sd(datos1), IQR(datos1, na.rm = TRUE) / 1.349)
  h2 = (4 / (3 * length(datos2)))^(1/5) * min(sd(datos2), IQR(datos2, na.rm = TRUE) / 1.349)
  h3 = (4 / (3 * length(datos3)))^(1/5) * min(sd(datos3), IQR(datos3, na.rm = TRUE) / 1.349)
  
  puntos = seq(min(datos1, datos2, datos3), max(datos1, datos2, datos3), length = 1000)
  return(sum(pmin(kernel.densidad(datos1, puntos, h1, kernel),
                  kernel.densidad(datos2, puntos, h2, kernel),
                  kernel.densidad(datos3, puntos, h3, kernel)))*(puntos[2] - puntos[1]))
}

# --- PROCESAMIENTO Y VISUALIZACIÓN ---

# Extraer grupo y variable
group <- AL$group
kfront <- -as.numeric(AL$kfront)

# Dividir por grupo
kfront.list <- split(kfront, group)
clase <- c("D-", "D0", "D+")
colores <- c("green", "blue", "red")

# Datos limpios por grupo
x <- na.exclude(as.numeric(kfront.list[[1]]))
y <- na.exclude(as.numeric(kfront.list[[2]]))
z <- na.exclude(as.numeric(kfront.list[[3]]))

# Bandwidths
h1 <- (4 / (3 * length(x)))^(1/5) * min(sd(x), IQR(x, na.rm = TRUE) / 1.349)
h2 <- (4 / (3 * length(y)))^(1/5) * min(sd(y), IQR(y, na.rm = TRUE) / 1.349)
h3 <- (4 / (3 * length(z)))^(1/5) * min(sd(z), IQR(z, na.rm = TRUE) / 1.349)

# Puntos comunes
puntos <- seq(-10, 10, length = 1000)

# Densidades por kernel
densidades <- list(
  D1 = list(x = puntos, y = kernel.densidad(x, puntos, h1, kernel_gauss)),
  D2 = list(x = puntos, y = kernel.densidad(y, puntos, h2, kernel_gauss)),
  D3 = list(x = puntos, y = kernel.densidad(z, puntos, h3, kernel_gauss))
)

# Rango común
xlim_range <- c(-10, 10)
ylim_range <- range(sapply(densidades, function(d) range(d$y)))

# Función para sombrear OVL
sombrear_ovl <- function(d1, d2, col, alpha = 0.3) {
  x_common <- d1$x
  y1 <- d1$y
  y2 <- d2$y
  y_min <- pmin(y1, y2)
  polygon(c(x_common, rev(x_common)),
          c(y_min, rep(0, length(y_min))),
          col = adjustcolor(col, alpha.f = alpha), border = NA)
}

# Iniciar gráfico vacío para sombrear primero
plot(NULL, xlim = xlim_range, ylim = ylim_range,
     main = "Kernel-based OVL",
     xlab = "kfront", ylab = "Density")

# Sombrear áreas de solapamiento (OVL)
sombrear_ovl(densidades[[1]], densidades[[2]], col = colores[1])
sombrear_ovl(densidades[[2]], densidades[[3]], col = colores[2])
sombrear_ovl(densidades[[1]], densidades[[3]], col = colores[3])

# Dibujar densidades por encima
lines(densidades[[1]]$x, densidades[[1]]$y, lty = "solid", lwd = 2, col = colores[1])
lines(densidades[[2]]$x, densidades[[2]]$y, lty = "dashed", lwd = 2, col = colores[2])
lines(densidades[[3]]$x, densidades[[3]]$y, lty = "dotdash", lwd = 2, col = colores[3])

# Leyenda
legend("topright", legend = clase,
       col = colores,
       lty = c("solid", "dashed", "dotdash"), lwd = 2)
