if (!require("BiocManager"))
        install.packages("BiocManager")

# BiocManager::install("OncoSimulR", version = "3.11")

if (!require("devtools"))
        install.packages("devtools")

# library(devtools)

# install_github("rdiaz02/OncoSimul/OncoSimulR", ref = "freq-dep-fitness")

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
     col = c("black", "red", "blue"), ylim = c(20, 50000),
     main = "No chemotherapy")
plot(simul, show = "genotypes", col = c("grey", "Red3", "Sky blue"),
     main = "No chemotherapy")


## SCENARIO WITH CONTINUOUS CHEMOTHERAPY (FIXED DOSE)

## Fitness definition

drug_eff <- 0.01  # effect of drug on sensitive cells

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
     thinData = TRUE, main = "Fixed dose")
plot(simul2, show = "genotypes", ylim = c(20, 10000), 
     col = c("grey", "Red3", "Sky blue"), main = "Fixed dose")


## FIRST APPROACH WITH ADAPTIVE THERAPY: FIXED DOSES 10 DAYS EVERY TWO DAYS

## Fitness definition

drug_eff2 <- 0.015 

fit_cells3 <-data.frame(Genotype = c("WT", "S", "R", "S, R"),
                        Fitness = c(wt_fitness,                         #WT
                                    paste0("if (T>48 & T<58)", 
                                           drug_eff2, "*(", sens_fitness,")",
                                           "; else if (T>60 & T<70) ", 
                                           drug_eff2, "*(", sens_fitness,")",
                                           "; else if (T>72 & T<82) ", 
                                           drug_eff2, "*(", sens_fitness,")",
                                           "; else if (T>84 & T<94) ", 
                                           drug_eff2, "*(", sens_fitness,")",
                                           "; else if (T>96 & T<106) ", 
                                           drug_eff2, "*(", sens_fitness,")",
                                           "; else if (T>108 & T<118) ", 
                                           drug_eff2, "*(", sens_fitness,")",
                                           "; else if (T>120 & T<130) ",
                                           drug_eff2, "*(", sens_fitness,")",
                                           "; else if (T>132 & T<142) ", 
                                           drug_eff2, "*(", sens_fitness,")",
                                           "; else if (T>144 & T<154) ", 
                                           drug_eff2, "*(", sens_fitness,")",
                                           "; else if (T>156 & T<166) ", 
                                           drug_eff2, "*(", sens_fitness,")",
                                           "; else if (T>168 & T<178) ", 
                                           drug_eff2, "*(", sens_fitness,")",
                                           "; else if (T>180 & T<190) ", 
                                           drug_eff2, "*(", sens_fitness,")",
                                           "; else if (T>192 & T<202) ", 
                                           drug_eff2, "*(", sens_fitness,")",
                                           "; else if (T>204 & T<214) ", 
                                           drug_eff2, "*(", sens_fitness,")",
                                           "; else ", sens_fitness, ";"),   #S
                                    "0",                                    #R
                                    res_fitness),                           #S,R
                        stringsAsFactors = FALSE)


all_fe3 <- allFitnessEffects(genotFitness = fit_cells3,
                             frequencyDependentFitness = TRUE,
                             frequencyType = "rel")

##Simulation

set.seed(2)

simul3  <- oncoSimulIndiv(all_fe3,
                          model = "McFL",
                          onlyCancer = FALSE,
                          finalTime = 220,
                          mu = 0.01,
                          initSize = 5000,
                          keepPhylog = FALSE,
                          seed = NULL)

## Plots: number of cells vs time
# ylim has been adapted to number of cells
plot(simul3, show = "genotypes", type = "line",
     col = c("black", "red", "blue"), ylim = c(20, 200000),
     main = "Fixed dose (10 d) every 2 days")
plot(simul3, show = "genotypes", ylim = c(20, 9000),
     col = c("grey", "Red3", "Sky blue"),
     main = "Fixed dose (10 d) every 2 days")


## SECOND APPROACH WITH ADAPTIVE THERAPY: TREATMENT ON/OFF DEPENDING
## ON TOTAL NUMBER OF CELLS

drug_eff2 <- 0.015

fit_cells4 <-data.frame(Genotype = c("WT", "S", "R", "S, R"),
                       Fitness = c(wt_fitness,                              #WT
                                   # the next expression uses a tag to turn
                                   # the drug effect on and off (it does not 
                                   # work because variables are not saved after
                                   # each iteration and drug is always assigned
                                   # 'true')
                                   paste0("var drug := true;
                                           if (T>50 & N>500 & drug = true)",
                                           drug_eff2, "* (", sens_fitness,");
                                           else if (T>50 & N=500) {
                                           drug := false;",
                                           sens_fitness, ";
                                           }
                                           else if (T>50 & N=2500) {
                                           drug := true;",
                                           drug_eff2, "* (", sens_fitness,");
                                           }
                                           else ", sens_fitness, ";"),      #S
                                   "0",                                     #R   
                                   res_fitness),                           #S,R
                       stringsAsFactors = FALSE)

fit_cells4

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
     col = c("black", "red", "blue"), ylim = c(20, 200000),
     main = "Attempt to turn treatment on/off")
plot(simul4, show = "genotypes", ylim = c(20, 10000),
     col = c("grey", "Red3", "Sky blue"),
     main = "Attempt to turn treatment on/off")


## THIRD APPROACH WITH ADAPTIVE THERAPY: INTERVENTION THAT DEPENDS ON TOTAL
## NUMBER OF CELLS

## Fitness definition

drug_eff2 <- 0.015

fit_cells5 <-data.frame(Genotype = c("WT", "S", "R", "S, R"),
                        Fitness = c(wt_fitness,                              #WT
                                    paste0("if (T>50 & N>500)", 
                                           drug_eff2, "* (", sens_fitness,")",
                                           "; else ", sens_fitness, ";"),    #S
                                    "0",                                     #R
                                    res_fitness),                            #S,R
                        stringsAsFactors = FALSE)

all_fe5 <- allFitnessEffects(genotFitness = fit_cells5,
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel")

## Simulation

set.seed(2)

simul5  <- oncoSimulIndiv(all_fe5,
                          model = "McFL",
                          onlyCancer = FALSE,
                          finalTime = 170,
                          mu = 0.01,
                          initSize = 5000,
                          keepPhylog = FALSE,
                          seed = NULL)

## Plots: number of cells vs time
# ylim has been adapted to number of cells
plot(simul5, show = "genotypes", type = "line",
     col = c("black", "red", "blue"), ylim = c(20, 200000),
     main = "Intervention depending on N")
plot(simul5, show = "genotypes", ylim = c(20, 10000),
     col = c("grey", "Red3", "Sky blue"),
     main = "Intervention depending on N")


## FOURTH APPROACH WITH ADAPTIVE THERAPY: CHANGING DRUG EFFECT AS A FUNCTION 
## OF TIME

## Fitness definition

fit_cells6 <-data.frame(Genotype = c("WT", "S", "R", "S, R"),
                        Fitness = c(wt_fitness,                         #WT
                                    paste0("if (T>50) ((0.41 + 2*sin(T+21)/5))
                                                   *(", sens_fitness, ");
                                       else ", sens_fitness, ";"),      #S
                                    "0",                                 #R
                                    res_fitness),                        #S,R
                        stringsAsFactors = FALSE)


all_fe6 <- allFitnessEffects(genotFitness = fit_cells6,
                             frequencyDependentFitness = TRUE,
                             frequencyType = "rel")

##Simulation

set.seed(2)

simul6  <- oncoSimulIndiv(all_fe6,
                          model = "McFL",
                          onlyCancer = FALSE,
                          finalTime = 170,
                          mu = 0.01,
                          initSize = 5000,
                          keepPhylog = FALSE,
                          seed = NULL)

## Plots: number of cells vs time
# ylim has been adapted to number of cells
plot(simul6, show = "genotypes", type = "line",
     col = c("black", "red", "blue"), ylim = c(20, 200000),
     main = "Intervention depending on T")
plot(simul6, show = "genotypes", ylim = c(20, 9000),
     col = c("grey", "Red3", "Sky blue"),
     main = "Intervention depending on T")


## FIFTH APPROACH WITH ADAPTIVE THERAPY: FIXED DOSE APPLICATION AS A FUNCTION
## OF TIME --> ~3 DAYS TREATMENT - ~3 DAYS NO TREATMENT

## Fitness definition

drug_eff2 <- 0.015

fit_cells7 <-data.frame(Genotype = c("WT", "S", "R", "S, R"),
                        Fitness = c(wt_fitness,                   #WT
                                    paste0("if (T>50 & sin(T)>0)", 
                                    drug_eff2, "*(", sens_fitness, ");
                                    else ", sens_fitness, ";"),    #S
                                    "0",                           #R
                                    res_fitness),                  #S,R
                        stringsAsFactors = FALSE)


all_fe7 <- allFitnessEffects(genotFitness = fit_cells7,
                             frequencyDependentFitness = TRUE,
                             frequencyType = "rel")

##Simulation

set.seed(2)

simul7  <- oncoSimulIndiv(all_fe7,
                          model = "McFL",
                          onlyCancer = FALSE,
                          finalTime = 170,
                          mu = 0.01,
                          initSize = 5000,
                          keepPhylog = FALSE,
                          seed = NULL)

## Plots: number of cells vs time
# ylim has been adapted to number of cells
plot(simul7, show = "genotypes", type = "line",
     col = c("black", "red", "blue"), ylim = c(20, 200000),
     main = "~3d treatment - ~3d no treatment")
plot(simul7, show = "genotypes", ylim = c(20, 9000),
     col = c("grey", "Red3", "Sky blue"),
     main = "~3d treatment - ~3d no treatment")


## SIXTH APPROACH WITH ADAPTIVE THERAPY: FIXED DOSE APPLICATION AS A FUNCTION
## OF TIME --> ~4.5 DAYS TREATMENT - ~2 DAYS NO TREATMENT

## Fitness definition

drug_eff2 <- 0.015

fit_cells8 <-data.frame(Genotype = c("WT", "S", "R", "S, R"),
                        Fitness = c(wt_fitness,                   #WT
                                    paste0("if (T>50 & (sin(T) + 0.6)>0)", 
                                           drug_eff2, "*(", sens_fitness, ");
                                    else ", sens_fitness, ";"),    #S
                                    "0",                           #R
                                    res_fitness),                  #S,R
                        stringsAsFactors = FALSE)


all_fe8 <- allFitnessEffects(genotFitness = fit_cells8,
                             frequencyDependentFitness = TRUE,
                             frequencyType = "rel")

##Simulation

set.seed(2)

simul8  <- oncoSimulIndiv(all_fe8,
                          model = "McFL",
                          onlyCancer = FALSE,
                          finalTime = 170,
                          mu = 0.01,
                          initSize = 5000,
                          keepPhylog = FALSE,
                          seed = NULL)

## Plots: number of cells vs time
# ylim has been adapted to number of cells
plot(simul8, show = "genotypes", type = "line",
     col = c("black", "red", "blue"), ylim = c(20, 200000),
     main = "~4.5d treatment - ~2d no treatment")
plot(simul8, show = "genotypes", ylim = c(20, 9000),
     col = c("grey", "Red3", "Sky blue"),
     main = "~4.5d treatment - ~2d no treatment")