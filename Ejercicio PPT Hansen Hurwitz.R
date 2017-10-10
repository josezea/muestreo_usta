library(TeachingSampling)

?S.PPS

x <- c(6, 2, 10, 1, 5)
set.seed(1234)
indicadora <- S.PPS(m = 3,x)[,1]
indicadora

library(pps)
?ppswr
ppswr(sizes = x, n = 3)

library(sampling)
data("MU284")
p_k <- MU284$ME84 / sum(MU284$ME84)
2135 /  sum(MU284$ME84)




datos_quiz <- read.delim("clipboard", dec = ",")
datos_quiz <- read.delim2("clipboard")
library(TeachingSampling)
?E.PPS
options(scipen = 99999)
E.PPS(y = datos_quiz$RMT85, pk = datos_quiz$p_k)
Z <- Domains(datos_quiz$Zona)
YZ <- datos_quiz$RMT85 * Z
datos <- cbind(Z, YZ, datos_quiz$RMT85)
colnames(datos)<- c("Z1", "Z2", "YZ1", "YZ2", "y")
head(datos)
E.PPS(y = datos, pk = datos_quiz$p_k)
