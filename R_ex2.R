library(OncoSimulR)
packageVersion("OncoSimulR")

r <- data.frame(Genotype = c("WT", "A", "B", "A, B"),
                Fitness = c("1",
                            "1 + 2*(f_A + f_B + 2*f_A_B)", #creo que eq a  "if (T>50) 1.5; else 0;"
                            "1 + 2*(f_A + f_B + 2*f_A_B)",
                            "2 + 4*(f_A + f_B + 2*f_A_B)"),
                stringsAsFactors = FALSE)

afe <- allFitnessEffects(genotFitness = r,
                        frequencyDependentFitness = TRUE,
                        frequencyType = "rel")
evalAllGenotypes(fe, spPopSizes = c(WT=100, A=100, B=100))
evalAllGenotypes(fe, spPopSizes = c(WT=100, A=100, B=100), currentTime = 50)
evalAllGenotypes(fe, spPopSizes = c(WT=100, A=100, B=100), currentTime = 51)


## Simulation
set.seed(1)
osi <- oncoSimulIndiv(afe,
                      model = "McFL",
                      onlyCancer = FALSE,
                      finalTime = 1000,
                      #verbosity = 0,
                      mu = c(A=1e-5, B=1e-6),
                      initSize = 1000,
                      keepPhylog = TRUE,
                      seed = NULL,
                      #detectionProb = NA,
                      #detectionSize = NA,
                      errorHitMaxTries = FALSE,
                      errorHitWallTime = FALSE)
osi
plot(osi, show = "genotypes", type = "line")
plotClonePhylog(osi, N = 0)

## target --> A

drug_eff <- 0.01

r2 <- data.frame(Genotype = c("WT", "A", "B", "A, B"),
                Fitness = c("1",
                            "0.1*(1 + 2*(f_A + f_B + 2*f_A_B))",
                            "1 + 2*(f_A + f_B + 2*f_A_B)",
                            "2 + 4*(f_A + f_B + 2*f_A_B)"),
                stringsAsFactors = FALSE)

afe2 <- allFitnessEffects(genotFitness = r2,
                          frequencyDependentFitness = TRUE,
                          frequencyType = "rel")

## Simulation
set.seed(1)
osi2 <- oncoSimulIndiv(afe2,
                      model = "McFL",
                      onlyCancer = FALSE,
                      finalTime = 5000,
                      #verbosity = 0,
                      mu = c(A=1e-5, B=1e-6),
                      initSize = 1000,
                      keepPhylog = TRUE,
                      seed = NULL,
                      #detectionProb = NA,
                      #detectionSize = NA,
                      errorHitMaxTries = FALSE,
                      errorHitWallTime = FALSE)
osi2
plot(osi2, show = "genotypes", type = "line")
plotClonePhylog(osi2, N = 0)
