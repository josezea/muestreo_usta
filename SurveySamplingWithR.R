# Todos los derechos reservados: José Fernando Zea ©
rm(list = ls())
library(sampling)
library(survey)
data(api)

#-------------  Simple random sampling without replacement -------------

# Interest variable: api00

# Sample size (Ideal conditions)
mu <- mean(apipop$api00); sd <- sd(apipop$api00); conf <- 0.95;  relme = 0.03
N <- nrow(apipop)
me <- relme * mu
z <- qnorm(Z <- 1 - ((1 - conf)/2))
n <- (z^2 * sd^2) / me^2
n_srswor  <-  ceiling(n / (1 + (n/N))) # Finite Correction

# Sampling selection
set.seed(12345)
selunits_srswor  <- sample(1:nrow(apipop), n_srswor, replace = F)
sample_srswor  <- apipop[selunits_srswor,] # Sample from a Simple Random Sampling Design

# Estimation using survey package
sample_srswor$N <- nrow(apipop)
sample_srswor$pw <- nrow(apipop) / nrow(sample_srswor) # sampling weight
sample_srswor$pi_k <- nrow(sample_srswor) / nrow(apipop) # first order probability national

design_srswor  <- svydesign(ids = ~1,  fpc = ~N, data = sample_srswor) 
estima_srswor <- svytotal(~api00, design_srswor, deff = T)
estima_srswor[1] / N
100 * cv(estima_srswor) # Coeficient of variation

# Other alternative
design_srswor  <- svydesign(ids = ~1,  fpc = ~pi_k, data = sample_srswor) 
(estima_srswor <- svytotal(~api00, design_srswor, deff = T))
estima_srswor[1] / N
100 * cv(estima_srswor) # Coeficient of variation

# fpc argument cannot be omited for without replacement survey design.

# Another ...
design_srswor  <- svydesign(ids = ~1,  weights = ~pw, fpc = ~N, data = sample_srswor) 
(estima_srswor <- svytotal(~api00, design_srswor, deff = T))
estima_srswor[1] / N
100 * cv(estima_srswor) # Coeficient of variation


# Another (2)...
design_srswor  <- svydesign(ids = ~1,  probs = ~pi_k, fpc = ~N, data = sample_srswor) 
(estima_srswor <- svytotal(~api00, design_srswor, deff = T))
estima_srswor[1] / N
100 * cv(estima_srswor) # Coeficient of variation

# -------------------------------------------------------------------- #


#-------------  Simple random sampling with replacement -------------
data(api)
set.seed(1234)

# Selection
selunits_srswr <- sample(1:nrow(apipop), n_srswor, replace = T)
sample_srswr <- apipop[selunits_srswr,] # Sample from a Simple Random Sampling Design

# Estimation
sample_srswr$N <- nrow(apipop)
sample_srswr$pw <- nrow(apipop) / nrow(sample_srswr) # sampling weight (N/m)
sample_srswr$pi_k <- nrow(sample_srswr) / nrow(apipop) # first order probability national


design_srswr <- svydesign(ids = ~1,  probs =  ~pi_k, data = sample_srswr) 
(estima_srswr <- svytotal(~api00, design_srswr, deff = T))
estima_srswr[1] / N
100 * cv(estima_srswr) # Coeficient of variation

#E.WR(N, m = n_srs, y = sample_srswr$api00)


# Another (Use horvitz - thompson estimator)
design_srswr <- svydesign(ids = ~1,  weights = ~pw,  data = sample_srswr) 
(estima_srswr <- svytotal(~api00, design_srswr, deff = T))
estima_srswr[1] / N
100 * cv(estima_srswr) # Coeficient of variation


# -------------------------------------------------------------------- #


#-------------  Bernoulli Design -------------

data(api)

# Selection

pi_k <- n_srswor / N
set.seed(12345)
indica_bern <- UPpoisson(rep(pi_k, N))
sample_bern <- apipop[indica_bern,]
sample_bern$pi_k <- pi_k

# Estimation
M_pi_kl <-  as.matrix(sample_bern$pi_k ) %*% t(as.matrix(sample_bern$pi_k ))
diag(M_pi_kl) <- sqrt(diag(M_pi_kl))
design_bern <- svydesign(id = ~1, fpc= ~pi_k, data = sample_bern, pps = ppsmat(M_pi_kl))
(estima_bern <- svytotal(~api00, design_bern))
estima_bern[1] / N
100 * cv(estima_bern) # Coeficient of variation

# -------------------------------------------------------------------- #


#-------------  Systematic Design -------------


# Selection

(a <- ceiling(N / n_srswor)) # Number of systematics groups
# 34 * n_srswor + 6 * (n_srswor - 1) 

pi_k <- 1 / a

set.seed(12345)
indica_syst<- as.logical(UPsystematic(pik = rep(pi_k, N)))

length(which(indica_syst))
which(indica_syst)

sample_syst <- apipop[indica_syst,]

# Estimation
sample_syst$pi_k <- pi_k
design_syst <- svydesign(id = ~1, fpc= ~pi_k, data = sample_syst)

(estima_syst <- svytotal(~api00, design_syst))
estima_syst[1] / N
100 * cv(estima_syst) # Coeficient of variation

# -------------------------------------------------------------------- #


#-------------  Poisson Design -------------

# Auxiliar variable: api99
# Interest variable: api00
data(api)
cor(apipop$api99, apipop$api00)
lm(apipop$api00~apipop$api99)


# Selection
apipop$pi_k <- inclusionprobabilities(a = apipop$api99, n = n_srswor)
set.seed(12345)
indica_poisson <- UPpoisson(apipop$pi_k )
sample_poisson <- apipop[indica_poisson,]

# Estimation
M_pi_kl <-  as.matrix(sample_poisson$pi_k ) %*% t(as.matrix(sample_poisson$pi_k ))
diag(M_pi_kl) <- sqrt(diag(M_pi_kl))

library(Matrix)
design_poisson <- svydesign(id = ~1, fpc= ~pi_k, data = sample_poisson, pps = ppsmat(M_pi_kl))
(estima_poisson <- svytotal(~api00, design_poisson))
estima_poisson[1] / N
100 * cv(estima_poisson) # Coeficient of variation

# -------------------------------------------------------------------- #


#--------------------------  piPS sampling design (Brewer Method) --------------------------
data(api)

# Selection
apipop$pi_k <- inclusionprobabilities(a = apipop$api99, n = n_srswor)
set.seed(12345)
indica_piPS <- as.logical(UPbrewer(apipop$pi_k ))
sample_piPS <- apipop[indica_piPS,]

# Estimation
M_pi_kl <-  as.matrix(sample_piPS$pi_k ) %*% t(as.matrix(sample_piPS$pi_k ))
diag(M_pi_kl) <- sqrt(diag(M_pi_kl))

library(Matrix)
design_piPS <- svydesign(id = ~1, fpc= ~pi_k, data = sample_piPS, pps = "brewer")
(estima_piPS <- svytotal(~api00, design_piPS))
estima_piPS[1] / N
100 * cv(estima_piPS) # Coeficient of variation

# -------------------------------------------------------------------- #


#--------------------------  pPS sampling design (with replacement) --------------------------
data(api)

# Selection

apipop$p_k <- apipop$api99 / sum(apipop$api99)
set.seed(1234)
selunits_pPS <- sample(1:nrow(apipop), n_srswor, replace = T, prob = apipop$p_k)
sample_pPS <- apipop[selunits_pPS,] # Sample from a Simple Random Sampling Design
sample_pPS$pw <- 1/(sample_pPS$p_k * n_srswor)


# Estimation
# Hansen - Horvitz estimator
design_pPS <- svydesign(ids = ~1,  weight =  ~pw, data = sample_pPS) 
(estima_pPS <- svytotal(~api00, design_pPS, deff = T))
estima_pPS[1] / N
100 * cv(estima_pPS) # Coeficient of variation
# ------------------------------------------------------------------------ #


#- Stratified Sampling (sample random sampling without replacement in each strata) -#
data(api)

# Selection
# Order dataset by strata
apipop <- apipop[order(apipop$stype),]

# Proportional allocation
n_h <- ceiling(n_srswor * (table(apipop$stype) / sum(table(apipop$stype))))

selunits_strata <- strata(data = apipop , stratanames = "stype",
       size = n_h, method = "srswor", description = T)

sample_strata <- getdata(apipop, selunits_strata)


# Estimation
design_strata <- svydesign(ids = ~1, strata = ~stype,  fpc =  ~Prob, data = sample_strata)
(estima_strata<- svytotal(~api00, design_strata, deff = T))
estima_strata[1] / N
100 * cv(estima_strata) # Coeficient of variation
# ------------------------------------------------------------------------ #


# ----- Cluster Sampling (Clusters are selected by sample random sampling ) ----- #


# We will select 25 cluster

# Selection
set.seed(12345)
selunits_cluster <- cluster(data = apipop, clustername = "dnum",
                           size = 25, method = "srswor", description = T)

sample_cluster <- getdata(apipop, selunits_cluster)
sample_cluster$N_I <- length(table(apipop$dnum))

# Estimation

options(survey.lonely.psu="adjust") # Adjust estimator when there is 1 unit in a cluster
design_cluster <- svydesign(ids = ~dnum, fpc =  ~N_I, data = sample_cluster)

(estima_cluster <- svytotal(~api00, design_cluster, deff = T))
estima_cluster[1] / N
100 * cv(estima_strata) # Coeficient of variation

