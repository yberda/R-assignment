library(OncoSimulR)


## NO CHEMO

set.seed(2)
RNGkind("L'Ecuyer-CMRG")


# Healthy     Sensitive   Resistant
a <- 1;       b <- 0.5;   c <- 0.5    # Healthy
d <- 1;       e <- 1.2;   f <- 0.7    # Sensitive
g <- 0.9;     h <- -0.5;  i <- 0.8    # Resistant


wt_fitness <- paste0(a, "*f_+", b, "*f_S+", c, "*f_S_R")
sens_fitness <- paste0(d,"*f_+",e,"*f_S+",f,"*f_S_R")
res_fitness <- paste0(g,"*f_+",h,"*f_S+",i,"*f_S_R")


## Fitness definition
players <- data.frame(Genotype = c("WT","S","R","S,R"),
                      Fitness = c(wt_fitness, #WT
                                  sens_finess,      #S
                                  "0",              #R
                                  res_fitness),     #S,R
                      stringsAsFactors = FALSE)


game <- allFitnessEffects(genotFitness = players,
                          frequencyDependentFitness = TRUE,
                          frequencyType = "rel")


## Plot: first scenario
eag <- evalAllGenotypes(game, spPopSizes = c(10,1,0,10))[c(1, 3, 4),]

plot(eag)

## Simulation
gamesimul <- oncoSimulIndiv(game,
                            model = "McFL",
                            onlyCancer = FALSE,
                            finalTime = 50,
                            mu = 0.01,
                            initSize = 5000,
                            keepPhylog = FALSE,
                            seed = NULL)
## Plot 2
plot(gamesimul, show = "genotypes", type = "line",
     col = c("black", "green", "red"), ylim = c(20, 50000))
plot(gamesimul, show = "genotypes")

##################################################################
####### CHEMO: FIXED DOSE

# Effect of drug on fitness sensible tumor cells

drug_eff <- 0.01  ## drug is only affecting sensible cells

wt_fitness <- paste0(a, "*f_+", b, "*f_S+", c, "*f_S_R")
sens_fitness <- paste0(d, "*f_+", e, "*f_S+", f, "*f_S_R")
res_fitness <- paste0(g, "*f_+", h, "*f_S+", i, "*f_S_R")

players_1 <- data.frame(Genotype = c("WT", "S", "R", "S, R"),
                        Fitness = c(wt_fitness,                                               #WT
                                    paste0("if (T>50) ", drug_eff, "*(",sens_fitness, ")",";
                                           else ", sens_fitness, ";"),                        #S
                                    "0",                                                      #R
                                    res_fitness),                                             #S,R
                        stringsAsFactors = FALSE)

period_1 <- allFitnessEffects(genotFitness = players_1,
                              frequencyDependentFitness = TRUE,
                              frequencyType = "rel")

set.seed(2)

final_time <- 170 ## for speed
simul_period_1 <- oncoSimulIndiv(period_1,
                                 model = "McFL",
                                 onlyCancer = FALSE,
                                 finalTime = final_time,
                                 mu = 0.01,
                                 initSize = 5000,
                                 keepPhylog = FALSE,
                                 seed = NULL)

# ylim has been adapted to number of cells
plot(simul_period_1, show = "genotypes", type = "line",
     col = c("black", "green", "red"), ylim = c(20, 300000),
     thinData = TRUE)
plot(simul_period_1, show = "genotypes", ylim = c(20, 12000))


############################################################################
### FIXED DOSE EVERY 2 DAYS - ADAPTATIVE THERAPY

set.seed(2)
RNGkind("L'Ecuyer-CMRG")

# Healthy     Sensitive   Resistant
a <- 1;       b <- 0.5;   c <- 0.5    # Healthy
d <- 1;       e <- 1.2;   f <- 0.7    # Sensitive
g <- 0.9;     h <- -0.5;  i <- 0.8    # Resistant

wt_fitness   <- paste0(a, "*f_+", b, "*f_S+", c, "*f_S_R")
sens_fitness <- paste0(d, "*f_+", e, "*f_S+", f, "*f_S_R")
res_fitness  <- paste0(g, "*f_+", h, "*f_S+", i, "*f_S_R")

fitness_df <-data.frame(Genotype = c("WT", "S", "R", "S, R"),
                        Fitness = c(wt_fitness,                                              #WT
                                    paste0("if (T>52 & T<60) 0.01 * (", sens_fitness,")",
                                          "; else if (T>62 & T<70) 0.01*(", sens_fitness,")",
                                          "; else if (T>72 & T<80) 0.01*(", sens_fitness,")",
                                          "; else if (T>82 & T<90) 0.01*(", sens_fitness,")",
                                          "; else if (T>92 & T<100) 0.01*(", sens_fitness,")",
                                          "; else if (T>102 & T<110) 0.01*(", sens_fitness,")",
                                          "; else if (T>112 & T<120) 0.01*(", sens_fitness,")",
                                          "; else if (T>122 & T<130) 0.01*(", sens_fitness,")",
                                          "; else if (T>132 & T<140) 0.01*(", sens_fitness,")",
                                          "; else if (T>142 & T<150) 0.01*(", sens_fitness,")",
                                          "; else if (T>152 & T<160) 0.01*(", sens_fitness,")",
                                          "; else if (T>162 & T<170) 0.01*(", sens_fitness,")",
                                          "; else if (T>172 & T<180) 0.01*(", sens_fitness,")",
                                          "; else if (T>182 & T<190) 0.01*(", sens_fitness,")",
                                          "; else if (T>192 & T<200) 0.01*(", sens_fitness,")",
                                          "; else if (T>202 & T<210) 0.01*(", sens_fitness,")",
                                          "; else ", sens_fitness, ";"),                    #S
                                    "0",                                                     #R
                                    res_fitness),                                            #S,R
                        stringsAsFactors = FALSE)

afe <- allFitnessEffects(genotFitness = fitness_df,
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel")

switching_sim  <- oncoSimulIndiv(afe,
                                 model = "McFL",
                                 onlyCancer = FALSE,
                                 finalTime = 220,
                                 mu = 0.01,
                                 initSize = 5000,
                                 keepPhylog = FALSE,
                                 seed = NULL)

plot(switching_sim, show = "genotypes", type = "line",
     col = c("black", "green", "red"), ylim = c(20, 200000))
plot(switching_sim, show = "genotypes", ylim = c(20, 1000))

############################################################################
### METRONOMIC THERAPY

set.seed(2)
RNGkind("L'Ecuyer-CMRG")

# Healthy     Sensitive   Resistant
a <- 1;       b <- 0.5;   c <- 0.5    # Healthy
d <- 1;       e <- 0.9;   f <- 0.6    # Sensitive
g <- 0.9;     h <- -0.5;  i <- 0.8    # Resistant

wt_fitness   <- paste0(a, "*f_+", b, "*f_S+", c, "*f_S_R")
sens_fitness <- paste0(d, "*f_+", e, "*f_S+", f, "*f_S_R")
res_fitness  <- paste0(g, "*f_+", h, "*f_S+", i, "*f_S_R")

fitness_df <-data.frame(Genotype = c("WT", "S", "R", "S, R"),
                        Fitness = c(wt_fitness,                                              #WT
                                    paste0("if (T>50 & T<80) 0.01 * (", sens_fitness,")",
                                           "; else if (T>125 & T<133) 0.01*(", sens_fitness,")",
                                           "; else if (T>170 & T<178) 0.01*(", sens_fitness,")",
                                           "; else if (T>215 & T<223) 0.01*(", sens_fitness,")",
                                           "; else ", sens_fitness, ";"),                    #S
                                    "0",                                                     #R
                                    res_fitness),                                            #S,R
                        stringsAsFactors = FALSE)

afe <- allFitnessEffects(genotFitness = fitness_df,
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel")

switching_sim  <- oncoSimulIndiv(afe,
                                 model = "McFL",
                                 onlyCancer = FALSE,
                                 finalTime = 230,
                                 mu = 0.01,
                                 initSize = 5000,
                                 keepPhylog = FALSE,
                                 seed = NULL)

plot(switching_sim, show = "genotypes", type = "line",
     col = c("black", "green", "red"), ylim = c(20, 200000))
plot(switching_sim, show = "genotypes", ylim = c(20, 1000))
