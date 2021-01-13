library(OncoSimulR)


##  COEFFICIENTS

# Healthy     Sensitive   Resistant
a <- 1;       b <- 0.5;   c <- 0.5    # Healthy
d <- 1;       e <- 1.2;   f <- 0.7    # Sensitive
g <- 0.9;     h <- -0.5;  i <- 0.8    # Resistant


wt_fitness <- paste0(a, "*f_+", b, "*f_S+", c, "*f_S_R")
sens_fitness <- paste0(d,"*f_+",e,"*f_S+",f,"*f_S_R")
res_fitness <- paste0(g,"*f_+",h,"*f_S+",i,"*f_S_R")


## SCENARIO WITHOUT CHEMOTHERAPY

## Fitness definition
fit_cells <- data.frame(Genotype = c("WT","S","R","S,R"),
                       Fitness = c(wt_fitness,       #WT
                                   sens_fitness,     #S
                                   "0",              #R
                                   res_fitness),     #S,R
                       stringsAsFactors = FALSE)


all_fe <- allFitnessEffects(genotFitness = fit_cells,
                            frequencyDependentFitness = TRUE,
                            frequencyType = "rel")


## Plot: first scenario
sc1 <- evalAllGenotypes(all_fe, spPopSizes = c(10,1,0,10))[c(1, 3, 4),]

plot(sc1)

## Simulation

set.seed(2)

simul <- oncoSimulIndiv(all_fe,
                            model = "McFL",
                            onlyCancer = FALSE,
                            finalTime = 50,
                            mu = 0.01,
                            initSize = 5000,
                            keepPhylog = FALSE,
                            seed = NULL)

## Plots: number of cells vs time
plot(simul, show = "genotypes", type = "line",
     col = c("black", "red", "blue"), ylim = c(20, 50000))
plot(simul, show = "genotypes", col = c("grey", "Red3", "Sky blue"))


## SCENARIO WITH CONTINUOUS CHEMOTHERAPY (FIXED DOSE)

## Fitness definition

drug_eff <- 0.01  # effect of drug on fitness sensible tumor cells

fit_cells2 <- data.frame(Genotype = c("WT", "S", "R", "S, R"),
                         Fitness = c(wt_fitness,                       #WT
                                      paste0("if (T>50) ", drug_eff, 
                                      "*(",sens_fitness, ")",";
                                      else ", sens_fitness, ";"),      #S
                                    "0",                               #R
                                    res_fitness),                      #S,R
                        stringsAsFactors = FALSE)

all_fe2 <- allFitnessEffects(genotFitness = fit_cells2,
                              frequencyDependentFitness = TRUE,
                              frequencyType = "rel")

## Simulation

set.seed(2)

simul2 <- oncoSimulIndiv(all_fe2,
                         model = "McFL",
                         onlyCancer = FALSE,
                         finalTime = 170,
                         mu = 0.01,
                         initSize = 5000,
                         keepPhylog = FALSE,
                         seed = NULL)


## Plots: number of cells vs time
# ylim has been adapted to number of cells
plot(simul2, show = "genotypes", type = "line",
     col = c("black", "red", "blue"), ylim = c(20, 300000),
     thinData = TRUE)
plot(simul2, show = "genotypes", ylim = c(20, 10000), 
     col = c("grey", "Red3", "Sky blue"))


## SCENARIO 1 WITH ADAPTIVE THERAPY: INTERVENTION THAT
## DEPENDS ON TIME

## Fitness definition


fit_cells3 <-data.frame(Genotype = c("WT", "S", "R", "S, R"),
                        Fitness = c(wt_fitness,                         #WT
                                paste0("if (T>50) ((0.41 + 2*sin(T+21)/5))
                                                   *(", sens_fitness, ");
                                       else ", sens_fitness, ";"),      #S
                                   "0",                                 #R
                                   res_fitness),                        #S,R
                        stringsAsFactors = FALSE)


all_fe3 <- allFitnessEffects(genotFitness = fit_cells3,
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel")

##Simulation

set.seed(2)

simul3  <- oncoSimulIndiv(all_fe3,
                          model = "McFL",
                          onlyCancer = FALSE,
                          finalTime = 170,
                          mu = 0.01,
                          initSize = 5000,
                          keepPhylog = FALSE,
                          seed = NULL)

## Plots: number of cells vs time
# ylim has been adapted to number of cells
plot(simul3, show = "genotypes", type = "line",
     col = c("black", "red", "blue"), ylim = c(20, 200000))
plot(simul3, show = "genotypes", ylim = c(20, 9000),
     col = c("grey", "Red3", "Sky blue"))


## SCENARIO 2 WITH ADAPTIVE THERAPY: INTERVENTION THAT
## DEPENDS ON TOTAL NUMBER OF CELLS

## Fitness definition

drug_eff2 <- 0.015

fit_cells4 <-data.frame(Genotype = c("WT", "S", "R", "S, R"),
                        Fitness = c(wt_fitness,                              #WT
                                    paste0("if (T>50 & N>500)", 
                                           drug_eff2, "* (", sens_fitness,")",
                                           "; else ", sens_fitness, ";"),    #S
                                    "0",                                     #R
                                    res_fitness),                            #S,R
                        stringsAsFactors = FALSE)

all_fe4 <- allFitnessEffects(genotFitness = fit_cells4,
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel")

## Simulation

set.seed(2)

simul4  <- oncoSimulIndiv(all_fe4,
                                 model = "McFL",
                                 onlyCancer = FALSE,
                                 finalTime = 170,
                                 mu = 0.01,
                                 initSize = 5000,
                                 keepPhylog = FALSE,
                                 seed = NULL)

## Plots: number of cells vs time
# ylim has been adapted to number of cells
plot(simul4, show = "genotypes", type = "line",
     col = c("black", "red", "blue"), ylim = c(20, 200000))
plot(simul4, show = "genotypes", ylim = c(20, 10000),
     col = c("grey", "Red3", "Sky blue"))

