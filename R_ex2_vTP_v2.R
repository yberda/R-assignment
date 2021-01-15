library(OncoSimulR)

# A --> androgen dependent cells
# B --> androgen producing cells
# C --> androgen independent cells


## NO TREATMENT

andr_eff <- 1.5  #androgen level

# Lotka-Volterra equation

LV_fx1 <- function(r1, r2, r3, K2, K3,
                    a_12, a_13, a_21, a_23, a_31, a_32, awt = 1e-4,
                    gt = c("WT", "S1", "S2", "S3")) {
  data.frame(Genotype = gt,
             Fitness = c(
              paste0("max(0.1, 1 - ", awt, " * (n_1 + n_2 + n_3))"),
              paste0("max(1, log(1 + n_1/(", andr_eff, "*n_2))) + ", r1,
                     "* ( 1 - (n_1 + ", a_12, " * n_2 + ", a_13, "* n_3)/(", 
                     andr_eff, "*n_2)",
                      ")"),
              paste0("max(1, log(1 + n_2/", K2, ")) + ", r2,
                     "* ( 1 - (n_2 + ", a_21, " * n_1 + ", a_23, "* n_3)/", K2,
                      ")"),
              paste0("max(1, log(1 + n_3/", K3, ")) + ", r3,
                     "* ( 1 - (n_3 + ", a_31, " * n_1 + ", a_32, "* n_2)/", K3,
                      ")")
             ))
}


fe1 <-
  allFitnessEffects(
    genotFitness =
      LV_fx1(r1 = 0.3, r2 = 0.4, r3 = 0.7, 
              K2 = 10000, K3 = 10000,
              a_12 = -0.5, a_13 = 0.1, 
              a_21 = 0, a_23 = 0.1, 
              a_31 = 0.6, a_32 = 0.5, 
              gt = c("WT","A", "B", "C")),
    frequencyDependentFitness = TRUE,
    frequencyType = "abs")

fe1

# Simulation

set.seed(2)

Simul1 <- oncoSimulIndiv(fe1,
                       model = "McFL",
                       onlyCancer = FALSE, 
                       finalTime = 170,
                       mu = 1e-4,
                       initSize = 10000, 
                       keepPhylog = TRUE,
                       seed = NULL, 
                       errorHitMaxTries = FALSE, 
                       errorHitWallTime = FALSE)

# Plot

plot(Simul1, show = "genotypes", type = "line",
     col = c("black", "green", "red", "blue"),
     thinData = TRUE,
     ylim = c(100, 100000))


## SCENARIO WITH CONTINUOUS TREATMENT: INHIBITION OF ADROGEN PRODUCTION AND
## ELIMINATION OF ANDROGEN-DEPENDENT AND PRODUCING CELLS

andr_eff <- 1.5  # androgen level in absence of treatment

andr_eff2 <- 0.5  # androgen level in presence of treatment 
                  # (inhibition of production)

drug_eff <- 0.01  # drug effect on androgen dependent and producing cells

# Lotka-Volterra equation

LV_fx2 <- function(r1, r2, r3, K2, K3,
                    a_12, a_13, a_21, a_23, a_31, a_32, awt = 1e-4,
                    gt = c("WT", "S1", "S2", "S3")) {
  data.frame(Genotype = gt,
             Fitness = c(
               paste0("max(0.1, 1 - ", awt, " * (n_1 + n_2 + n_3))"),
               
               paste0("if (T>100) max(1, 
                      log(1 + N/(", andr_eff2*drug_eff, "*n_2))) + ", r1,
                      "* ( 1 - (n_1 + ", a_12, " * n_2 + ", a_13, "* n_3)/(", 
                      andr_eff2*drug_eff, "*n_2));
                      
                      else max(1, log(1 + N/(", andr_eff, "*n_2))) + ", r1,
                      "* ( 1 - (n_1 + ", a_12, " * n_2 + ", a_13, "* n_3)/(", 
                      andr_eff, "*n_2))"),
               
               paste0("if (T>100) max(1, log(1 + N/", drug_eff*K2, ")) + ", r2,
                      "* ( 1 - (n_2 + ", a_21, " * n_1 + ", a_23, "* n_3))/(",
                      drug_eff*K2, ");
                      
                      else max(1, log(1 + N/", K2, ")) + ", r2,
                      "* ( 1 - (n_2 + ", a_21, " * n_1 + ", a_23, "* n_3)/", 
                      K2,")"),
               
               paste0("max(1, log(1 + N/", K3, ")) + ", r3,
                      "* (1 - (n_3 + ", a_31, " * n_1 + ", a_32, "* n_2)/", K3,
                      ")")
             ))
}

fe2 <-
  allFitnessEffects(
    genotFitness =
      LV_fx2(r1 = 0.3, r2 = 0.4, r3 = 0.7, 
              K2 = 10000, K3 = 10000,
              a_12 = -0.5, a_13 = 0.1, 
              a_21 = 0, a_23 = 0.1, 
              a_31 = 0.6, a_32 = 0.5, 
              gt = c("WT","A", "B", "C")),
    frequencyDependentFitness = TRUE,
    frequencyType = "abs")

fe2

# Simulation

set.seed(2)

Simul2 <- oncoSimulIndiv(fe2,
                       model = "McFL",
                       onlyCancer = FALSE, 
                       finalTime = 170,
                       mu = 1e-4,
                       initSize = 10000, 
                       keepPhylog = TRUE,
                       seed = NULL, 
                       errorHitMaxTries = FALSE, 
                       errorHitWallTime = FALSE)

# Plot

plot(Simul2, show = "genotypes", type = "line",
     col = c("black", "green", "red", "blue"),
     thinData = TRUE,
     ylim = c(100, 100000))


## SCENARIO WITH ADAPTIVE THERAPY: SWITCHING DRUG DOSES

andr_eff <- 1.5  # androgen level in absence of treatment

# (1 + sin(T + 21)/2) --> androgen level ranges from 0.5 to 1.5, depending 
# on whether the drug is applied or withdrawn 

# (0.41 + 2 * sin(T + 21)/5) --> drug effect on androgen dependent and 
# producing cells (ranges from approximately 0.01 to 0.8, depending on
# whether the drug is applied or withdrawn)


# Lotka-Volterra equation

LV_fx3 <- function(r1, r2, r3, K2, K3,
                    a_12, a_13, a_21, a_23, a_31, a_32, awt = 1e-4,
                    gt = c("WT", "S1", "S2", "S3")) {
  data.frame(Genotype = gt,
          Fitness = c(
            paste0("max(0.1, 1 - ", awt, " * (n_1 + n_2 + n_3))"),
            
            paste0("if (T > 100 & n_3 < 1500) max(1, log(1 + 
                   N/((0.41 + 2*sin(T+21)/5) *(1 + sin(T + 21)/2) *n_2)))+",
                   r1, "* ( 1 - (n_1 + ", a_12, " * n_2 + ", a_13, "* n_3)/
                   (( 0.41 + 2 *sin(T + 21)/5) *(1 + sin(T + 21)/2) *n_2));
                   
                   else max(1, log(1 + N/(", andr_eff, "* n_2))) + ", r1,
                   "* ( 1 - (n_1 + ", a_12, " * n_2 + ", a_13, "* n_3)/(", 
                   andr_eff, "* n_2))"),
            
            paste0("if (T > 100 & n_3 < 1500) max(1, log(1 +
                   N/((0.41 + 2 * sin(T + 21)/5)*", K2, "))) + ", r2,
                   "* ( 1 - (n_2 + ", a_21, " * n_1 + ", a_23, "* n_3)/",
                   "(0.41 + 2 * sin(T + 21)/5)*", K2, ");
                   
                   else max(1, log(1 + N/", K2, ")) + ", r2,
                   "*( 1 - (n_2 + ", a_21, " *n_1 + ", a_23, "*n_3)/", K2,")"),
            
            paste0("max(1, log(1 + N/", K3, ")) + ", r3,
                   "* ( 1 - (n_3 + ", a_31, " *n_1 + ", a_32, "*n_2)/", K3,")")
             ))
}

fe3 <-
  allFitnessEffects(
    genotFitness =
      LV_fx3(r1 = 0.3, r2 = 0.4, r3 = 0.7, 
              K2 = 10000, K3 = 10000,
              a_12 = -0.5, a_13 = 0.1, 
              a_21 = 0, a_23 = 0.1, 
              a_31 = 0.6, a_32 = 0.5, 
              gt = c("WT","A", "B", "C")),
    frequencyDependentFitness = TRUE,
    frequencyType = "abs")

fe3

# Simulation

set.seed(2)

Simul3 <- oncoSimulIndiv(fe3,
                       model = "McFL",
                       onlyCancer = FALSE, 
                       finalTime = 250,
                       mu = 1e-4,
                       initSize = 10000, 
                       keepPhylog = TRUE,
                       seed = NULL, 
                       errorHitMaxTries = FALSE, 
                       errorHitWallTime = FALSE)

# Plot

plot(Simul3, show = "genotypes", type = "line",
     col = c("black", "green", "red", "blue"),
     thinData = TRUE,
     ylim = c(100, 100000))
