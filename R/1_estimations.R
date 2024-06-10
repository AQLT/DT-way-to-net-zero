source("R/0_functions.R")
if (!dir.exists("est_mod")) 
  dir.create("est_mod")
library(future)
library(dplyr)
library(stringr)
plan(multisession, workers= min(length(future::availableWorkers()), 10))

# Fixation des paramètres
list_parameters <- read.csv("R/parameters.csv", encoding = "utf-8")
ntir <- 50 

for (parameters_nb in 1:nrow(list_parameters)) {
  set.seed(151071)
  parameters <- list_parameters[list_parameters$param_nb == parameters_nb,]
  for (name in colnames(list_parameters)){
    assign(name, parameters[[name]])
  }
  ### Paramètres de base supplémentaires
  Tdeb <- as.numeric(str_split_fixed(parameters$seqt, ":",2)[1])
  Tfin <- as.numeric(str_split_fixed(parameters$seqt, ":",2)[2])
  seqt <- Tdeb:Tfin
  nmax <- length(seqt) # utile lorsque l'on prend un pas de temps qui augmente
  nbpoints <- 2*parameters$Tmax+nmax # plus d'investissement ni de dépreciation du brun après Tmax et pas de dépréciation du vert
  delta <- c(delta_1 = parameters$delta_1, delta_2 = parameters$delta_2) # taux de dépreciation du capital
  ## deprecated : e_t_l_max <- parameters$E * (1-0.38) # cible : baisse de 55 % des émissions en 2030 p/r à 1990
  e_t_l_max <- parameters$E_net_1990 * (1-0.55) # cible : baisse de 55 % des émissions nettes en 2030 p/r à 1990.
  e_t_l_max2 <- parameters$E_net_1990 * (1-0.9)  # cible : baisse de 90% des émissions nettes en 2040 p/r à 1990
  
#  sigma <- parameters$sigma
  
  # choix d'une valeur de sigma compatible avec un cout d'ajustement marginal
#  CA_fin <- 0.776
  
  f1 = sprintf("est_mod/simul_rev_%s_%s_S1.RDS",parameters_nb,list_parameters[list_parameters$param_nb == parameters_nb, "param_name"])
  f2 = sprintf("est_mod/simul_rev_%s_%s_S2.RDS",parameters_nb,list_parameters[list_parameters$param_nb == parameters_nb, "param_name"])
  f2b = sprintf("est_mod/simul_rev_%s_%s_S2b.RDS",parameters_nb,list_parameters[list_parameters$param_nb == parameters_nb, "param_name"])
  f3 = sprintf("est_mod/simul_rev_%s_%s_S3.RDS",parameters_nb,list_parameters[list_parameters$param_nb == parameters_nb, "param_name"])
  f4 = sprintf("est_mod/simul_rev_%s_%s_S4.RDS",parameters_nb,list_parameters[list_parameters$param_nb == parameters_nb, "param_name"])
  if (file.exists(f1)) 
    next

  initpar <- init_param(rho = parameters$rho,
                        delta = delta,
                        sigma = parameters$sigma,
                        alpha0 = 1,
                        beta0 = 1,
                        pib = parameters$pib,
                        K_t = parameters$K_t,
                        part_K_b = parameters$part_K_b,
                        epuits = parameters$epuits,
                        Emax = parameters$Emax,
                        E = parameters$E,
                        e_t_l_max = e_t_l_max)

  k0 <- initpar$k0
  i0 <- initpar$i0
  E0 <- initpar$E0
  c0 <- initpar$c0
  em <- initpar$em
  a <- initpar$a

# le coefficient du coût d'ajustement est rajouté à la fin de a, en position 6  
  a <- c(initpar$a,parameters$Cout_Ajust)
  
  is <- initpar$is
  kminb <- a[3]
  alpha0 <- 1

  init <- list()
  for (j in 1:ntir) {
    i2_i <- runif(nbpoints,0,1)
    i2_i[2*Tmax+(1:nmax)] <- is[2]
    init[[j]] <- i2_i
  }

  ########################
  # Zen 2050, Scénario 1 #
  ########################
  
  # Zero contrainte, émissions nettes limitées au puits de carbone en 2050
  Emax <- 20
  lambda0 <- c(0,0)
  Elimit <- rep(Emax, nmax+1)  #On contraint toutes les dates à Emax, comme borne supérieure
  Elimit[(parameters$Tmax+2):(nmax+1)] <- 0

  if (a[6]>0) {
    # On prend la baseline avec coût d'ajustement nul comme jeu de paramètres initiaux possible
    f0 = sprintf("est_mod/simul_rev_%s_%s_S1.RDS",1,list_parameters[list_parameters$param_nb == 1, "param_name"])    
    res0 <- readRDS(f0)
    init[[ntir+1]] <- res0$par
  }
  
  res1 <- estim_robust(init,1,lambda0,Elimit,seqt,alpha0,k0,E0,em,a,parameters$Tmax,delta,is,rho)  
  res1$Elimit <- Elimit
  res1$lambda <- lambda0
  initpar$Emax <- Emax
  res1$initpar <- initpar
  res1$parameters <- parameters

  saveRDS(res1, f1)
  voir1 <- format_result(res1)
  Emismin <- max(voir1$Em_stock)
  
  ################################
  # Fit for 55 cappé, Scénario 2 #
  ################################
  lambda0 <- c(0,0)
  
  ElimitF55 <- Elimit
  ElimitF55[(tl+1):(Tmax+1)] <- e_t_l_max
  
  init[[ntir+1]] <- res1$par
  
  resF55 <- estim_robust(init,1,lambda0,ElimitF55,seqt,alpha0,k0,E0,em,a,Tmax,delta,is,rho)  
  
  resF55$Elimit <- ElimitF55
  resF55$lambda <- lambda0
  initpar$Emax <- Emax
  resF55$initpar <- initpar
  resF55$parameters <- parameters
  
  #Lissage du résultat pour éviter des phénomènes d'oscillation
  
  param <- resF55$par
  param[tl:(tl+9)] <- mean(resF55$par[tl:(tl+9)])

  res2 <- optim(param,crit,numdeb=1,lambda=lambda0, Elimit=ElimitF55, seqt = seqt,
                    alpha0 = alpha0, k0 = k0, E0 = E0,
                    em = em, a = a, Tmax = Tmax, delta = delta, is = is, rho = rho,
                    lower=rep(0,4*nmax),upper=rep(1,4*nmax),method="L-BFGS-B", ) #gr=dcrit)
  res2$Elimit <- ElimitF55
  res2$lambda <- lambda0
  initpar$Emax <- Emax
  res2$initpar <- initpar
  res2$parameters <- parameters
  
  voir2 <- format_result(res2)
  max(voir2$Em_stock)
  
  saveRDS(res2, f2)
  

  #####################################
  # Fit for 55 non cappé, Scénario 2b #
  #####################################
  lambda0 <- c(0,0)
  
  ElimitF55b <- Elimit
  ElimitF55b[(tl+1)] <- e_t_l_max
  ElimitF55b[(tl+2):(Tmax+1)] <- Emax
  
  init[[ntir+1]] <- res1$par
  init[[ntir+2]] <- res2$par
  
  res2b <- estim_robust(init,1,lambda0,ElimitF55b,seqt,alpha0,k0,E0,em,a,Tmax,delta,is,rho)  
  
  res2b$Elimit <- ElimitF55b
  res2b$lambda <- lambda0
  res2b$initpar <- initpar
  res2b$parameters <- parameters
  
  voir2b <- format_result(res2b)
  max(voir2b$Em_stock)
  saveRDS(res2b, f2b)
  
  #############################################
  # Fit for 55 + Fit for 90 cappé, Scénario 3 #
  #############################################
  
  lambda0 <- c(0,0)
  
  # verifier les annees : tl:tl2-1 plutôt ?
  ElimitF5590 <- Elimit
  ElimitF5590[(tl+1):(tl2)] <- e_t_l_max
  ElimitF5590[(tl2+1):(Tmax+1)] <- e_t_l_max2
  
  init[[ntir+1]] <- res2$par
  
  resF5590 <- estim_robust(init,1,lambda0,ElimitF5590,seqt,alpha0,k0,E0,em,a,Tmax,delta,is,rho)  
  
  resF5590$Elimit <- ElimitF5590
  resF5590$lambda <- lambda0
  initpar$Emax <- Emax
  resF5590$initpar <- initpar
  resF5590$parameters <- parameters
  
  
  res3 <- resF5590
  
  voir3 <- format_result(res3)
  max(voir3$Em_stock)
  
  saveRDS(res3, f3)
  
  ##############################
  # Budget Carbone, Scénario 4 #
  ##############################
  
  Emax <- parameters$Emax
  Tmax <- parameters$Tmax
  init[[ntir+2]] <- res1$par
  
  if (Emismin>Emax) {
    lambdamin <- 0
    lambdamax <- 0.5
    nstop <- 0
    #recherche d'une valeur de lambda pour laquelle les émissions sont plus basses que Emax
    repeat{
      
      res4 <- estim_robust(init,1,c(lambdamax,0),Elimit,seqt,alpha0,k0,E0,em,a,Tmax,delta,is,rho)  
      
      res4$Elimit <- Elimit
      res4$lambda <- c(lambdamax,0)
      initpar$Emax <- 20
      res4$initpar <- initpar
      res4$parameters <- parameters
      
      voir4 <- format_result(res4)
      Emismax <- max(voir4$Em_stock)
      
      if (Emismax<Emax) {break}
      nstop <- nstop+1
      if (nstop>10) {break}
      lambdamin <- lambdamax
      Emismin <- Emismax
      lambdamax <- 2*lambdamax
    }
    # Recherche de la valeur de lambda qui sature le budget carbone
    nstop <- 0
    repeat{
      lambdae <- lambdamin+(lambdamax-lambdamin)*abs((Emax-Emismin)/(Emismax-Emismin))
      
      init[[ntir+2]] <- res4$par
      res4 <- estim_robust(init,1,c(lambdae,0),Elimit,seqt,alpha0,k0,E0,em,a,Tmax,delta,is,rho)  
      
      res4$Elimit <- Elimit
      res4$lambda <- c(lambdae,0)
      initpar$Emax <- 20
      res4$initpar <- initpar
      res4$parameters <- parameters
      
      voir4 <- format_result(res4)
      Emis <- max(voir4$Em_stock)
      if (Emis<Emax) {
        lambdamax <- lambdae
        Emismax <- Emis
      }
      if (Emis>=Emax) {
        lambdamin <- lambdae
        Emismin <- Emis
      }
      print(c(list_parameters[list_parameters$param_nb == parameters_nb, "param_name"],
              nstop,lambdae,Emis))
      if ((abs(lambdamax-lambdamin)<1e-4) & (Emax > Emis)) 
        break
      nstop <- nstop+1
      if (nstop>=30) {break}
    }
    saveRDS(res4, f4)
  } else {
    res4 <- res1
    saveRDS(res4, f4)    
  }
}
  




