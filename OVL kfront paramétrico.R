# Extraer grupo y variable
group <- AL$group
kfront <- -as.numeric(AL$kfront)

# Obtener valores de los grupos sin NA
x <- na.exclude(as.numeric(kfront.list[[1]]))
y <- na.exclude(as.numeric(kfront.list[[2]]))
z <- na.exclude(as.numeric(kfront.list[[3]]))

# Calcular medias y desviaciones estándar
mu1 <- mean(x); sd1 <- sd(x)
mu2 <- mean(y); sd2 <- sd(y)
mu3 <- mean(z); sd3 <- sd(z)

# Crear un rango común para el eje x
x_seq <- seq(-10, 10, length.out = 1000)

# Calcular las densidades normales teóricas
dens1 <- dnorm(x_seq, mean = mu1, sd = sd1)
dens2 <- dnorm(x_seq, mean = mu2, sd = sd2)
dens3 <- dnorm(x_seq, mean = mu3, sd = sd3)

# Rango para el eje y
ylim_range <- range(0, dens1, dens2, dens3)

# Función para sombrear OVL
sombrear_ovl <- function(d1, d2, col, alpha = 0.3) {
  y_min <- pmin(d1, d2)
  polygon(c(x_seq, rev(x_seq)),
          c(y_min, rep(0, length(y_min))),
          col = adjustcolor(col, alpha.f = alpha), border = NA)
}

# Iniciar gráfico vacío
plot(NULL, xlim = c(-10,10), ylim = ylim_range,
     main = "Parametric OVL",
     xlab = "kfront", ylab = "Density")

# Sombrear áreas de solapamiento
sombrear_ovl(dens1, dens2, col = colores[1])
sombrear_ovl(dens2, dens3, col = colores[2])
sombrear_ovl(dens1, dens3, col = colores[3])

# Dibujar curvas normales teóricas
lines(x_seq, dens1, col = colores[1], lwd = 2, lty = "solid")
lines(x_seq, dens2, col = colores[2], lwd = 2, lty = "dashed")
lines(x_seq, dens3, col = colores[3], lwd = 2, lty = "dotdash")

# Leyenda
legend("topright", legend = clase,
       col = colores,
       lty = c("solid", "dashed", "dotdash"), lwd = 2)

