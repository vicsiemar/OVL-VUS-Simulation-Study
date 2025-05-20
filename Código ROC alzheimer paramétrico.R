# Representacion de la superfice ROC

library(plotly)
library(distr)

group <- AL$group
AL <- subset(AL, select = -group)
AL <- -AL

kfront <- as.numeric(AL$kfront)
kfront.list <- split(kfront, group)

roc_normal <- function(p1, p3) {
  result <- ifelse(qnorm(p1, mu1, sd1) < qnorm(1 - p3, mu3, sd3),
                   pnorm(qnorm(1 - p3, mu3, sd3), mu2, sd2) - pnorm(qnorm(p1, mu1, sd1), mu2, sd2),
                   0)
  return(result)
}

x=na.exclude(as.numeric(kfront.list[[1]]))
y=na.exclude(as.numeric(kfront.list[[2]]))
z=na.exclude(as.numeric(kfront.list[[3]]))

mu1=mean(x); mu2=mean(y);mu3=mean(z)
sd1=sd(x); sd2=sd(y); sd3=sd(z)

p1 <- seq(0, 1, 0.01)
p3 <- p1 

z1 <- outer(1-p3, p1, roc_normal)

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