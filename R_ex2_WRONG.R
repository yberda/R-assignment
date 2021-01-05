library(OncoSimulR)
packageVersion("OncoSimulR")

r <- data.frame(Genotype = c("WT", "A", "B", "A, B"),
                Fitness = c("1",
                            "1 + 2*(f_A + f_B + 2*f_A_B)", 
                            "1 + 2*(f_A + f_B + 2*f_A_B)",
                            "2 + 4*(f_A + f_B + 2*f_A_B)"),
                stringsAsFactors = FALSE)

afe <- allFitnessEffects(genotFitness = r,
                        frequencyDependentFitness = TRUE,
                        frequencyType = "rel")
afe
evalAllGenotypes(afe, spPopSizes = c(100, 100, 100, 100))
evalAllGenotypes(afe, spPopSizes = c(100, 100, 100, 100), currentTime = 1000)

## Simulation
set.seed(1)
osi <- oncoSimulIndiv(afe,
                      model = "McFL",
                      onlyCancer = FALSE,
                      finalTime = 20,
                      #verbosity = 0,
                      mu = c(A=0.01, B=0.001),
                      initSize = 5000,
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

r2 <- data.frame(Genotype = c("WT", "A", "B", "A, B"),
                Fitness = c("1",
                            "0.01*(1 + 2*(f_A + f_B + 2*f_A_B))",
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
                      finalTime = 50,
                      #verbosity = 0,
                      mu = c(A=0.01, B=0.001),
                      initSize = 5000,
                      keepPhylog = TRUE,
                      seed = NULL,
                      #detectionProb = NA,
                      #detectionSize = NA,
                      errorHitMaxTries = FALSE,
                      errorHitWallTime = FALSE)
osi2
plot(osi2, show = "genotypes", type = "line")
plotClonePhylog(osi2, N = 0)

## target --> B

r3 <- data.frame(Genotype = c("WT", "A", "B", "A, B"),
                 Fitness = c("1",
                             "1 + 2*(f_A + f_B + 2*f_A_B)",
                             "0.01*(1 + 2*(f_A + f_B + 2*f_A_B))",
                             "2 + 4*(f_A + f_B + 2*f_A_B)"),
                 stringsAsFactors = FALSE)

afe3 <- allFitnessEffects(genotFitness = r3,
                          frequencyDependentFitness = TRUE,
                          frequencyType = "rel")

## Simulation
set.seed(1)
osi3 <- oncoSimulIndiv(afe3,
                       model = "McFL",
                       onlyCancer = FALSE,
                       finalTime = 50,
                       #verbosity = 0,
                       mu = c(A=0.01, B=0.001),
                       initSize = 5000,
                       keepPhylog = TRUE,
                       seed = NULL,
                       #detectionProb = NA,
                       #detectionSize = NA,
                       errorHitMaxTries = FALSE,
                       errorHitWallTime = FALSE)
osi3
plot(osi3, show = "genotypes", type = "line")
plotClonePhylog(osi3, N = 0)

## target --> A + B

r4 <- data.frame(Genotype = c("WT", "A", "B", "A, B"),
                 Fitness = c("1",
                             "0.01*(1 + 2*(f_A + f_B + 2*f_A_B))",
                             "0.01*(1 + 2*(f_A + f_B + 2*f_A_B))",
                             "2 + 4*(f_A + f_B + 2*f_A_B)"),
                 stringsAsFactors = FALSE)

afe4 <- allFitnessEffects(genotFitness = r4,
                          frequencyDependentFitness = TRUE,
                          frequencyType = "rel")

## Simulation
set.seed(1)
osi4 <- oncoSimulIndiv(afe4,
                       model = "McFL",
                       onlyCancer = FALSE,
                       finalTime = 500,
                       #verbosity = 0,
                       mu = c(A=0.01, B=0.001),
                       initSize = 5000,
                       keepPhylog = TRUE,
                       seed = NULL,
                       #detectionProb = NA,
                       #detectionSize = NA,
                       errorHitMaxTries = FALSE,
                       errorHitWallTime = FALSE)
osi4
plot(osi4, show = "genotypes", type = "line")
plotClonePhylog(osi4, N = 0)
