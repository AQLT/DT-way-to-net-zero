source("R/0_functions.R")
library(future)
library(dplyr)
library(stringr)
list_parameters <- read.csv("R/parameters.csv", encoding = "utf-8")

###############################################################
####### Budget carbone différé (2023, 2028, 2033, 2038) #######
###############################################################

#Choix de la spécification à partir de laquelle on construit les scénarios
# Pour budget 6.23, specif 82

param_nb <- which(list_parameters$Emax==6.23)
f0 <- paste("est_mod/simul_rev_",param_nb,"_var_Emax_S1.RDS",sep="")
f4 <- paste("est_mod/simul_rev_",param_nb,"_var_Emax_S4.RDS",sep="")

f0 <- "est_mod/simul_rev_1_baseline_S1.RDS"
f4 <- "est_mod/simul_rev_1_baseline_S4.RDS"
res4 <- readRDS(f4)

# Trajectoire ZEN de base, à partir de laquelle on part pour faire le BC

res0 <- readRDS(f0)
tab0 <- format_result(res0)

parameters <- res0$parameters
for (name in colnames(parameters)){
  assign(name, parameters[[name]])
}
Emax <- res0$parameters$Emax
saveRDS(res4, sprintf("est_mod/ZEN_BC%s_Emax_%s_S4.RDS",2023,Emax))


alpha0 <- 1
beta0 <-1
nmax <- length(seqt)

Tdeb <- as.numeric(str_split_fixed(parameters$seqt, ":",2)[1])
Tfin <- as.numeric(str_split_fixed(parameters$seqt, ":",2)[2])
seqt <- Tdeb:Tfin
nmax <- length(seqt) # utile lorsque l'on prend un pas de temps qui augmente
nbpoints <- 2*parameters$Tmax+nmax # plus d'investissement ni de dépreciation du brun après Tmax et pas de dépréciation du vert
delta <- c(delta_1 = parameters$delta_1, delta_2 = parameters$delta_2) # taux de dépreciation du capital
e_t_l_max <- parameters$E * (1-0.38) # cible : baisse de 55 % des émissions en 2030 p/r à 1990

initpar <- res0$initpar
k0 <- initpar$k0
i0 <- initpar$i0
E0 <- initpar$E0
c0 <- initpar$c0
em <- initpar$em
a <- initpar$a
is <- initpar$is

a[6] <- parameters$Cout_Ajust

#Sentier: passage à un budget carbone en 2028,33, ... cf. table 5 du papier
set.seed(151071)

Elimit <- rep(Emax, nmax+1)  #On contraint toutes les dates à Emax, comme borne supérieure
Elimit[(Tmax+2):(nmax+1)] <- 0

ldate_BC <- c(2023,2028,2033,2038)
n_max <- list()
for(date_BC_i in seq_along(ldate_BC)) {
  set.seed(151071)
  date_BC <- ldate_BC[date_BC_i]
  print(date_BC)
  fl = sprintf("est_mod/ZEN_BC%s_Emax_%s_S4.RDS",date_BC,Emax)
  if (file.exists(fl)) 
    next
  
  numdeb <- date_BC-2023+1
  res_ref <- res0
  tab_ref <- format_result(res_ref)
  
  kn <- unlist(tab_ref[numdeb+1,paste("capital",1:2,sep="")])
  En <- tab_ref[numdeb+1,"Em_stock"]
  alphan <- alpha0/(1+rho)^(numdeb-1)
  lambda <- c(0,0)
  Emismin <- max(tab_ref$Em_stock)
  
  lambdae <- 1
  nstop <- 1
  ntir <- 10
  repeat{
    
    init <- list()
    for (j in 1:ntir) {
      i2_i <- runif(nbpoints,0,1)
      i2_i[1:numdeb] <- res_ref$par[1:numdeb]
      i2_i[Tmax+1:numdeb] <- res_ref$par[Tmax+1:numdeb]
      i2_i[Tmax+nmax+1:numdeb] <- res_ref$par[Tmax+nmax+1:numdeb]
      init[[j]] <- i2_i
    }
    
    init[[length(init)+1]] <- res_ref$par
    
    res2 <- estim_robust(init,numdeb,c(lambdae,0),Elimit,seqt,alphan,kn,En,em,a,Tmax,delta,is,rho)
    
    res2$Elimit <- Elimit
    res2$lambda <- c(lambdae,0)
    res2$initpar <- res_ref$initpar
    res2$parameters <- res_ref$parameters
    
    tab2 <- format_result(res2)
    Emismax <- max(tab2$Em_stock)
    print(c(nstop,lambdae,Emismax))
    
    if (Emismax<Emax | nstop>=20) {break}
    nstop <- nstop+1
    lambdae <- 2*lambdae
  }
  
  if (Emismax > Emax) {
    warning("Le fichier suivant n'a pas pu être calculé :\n", fl)
    next 
  }
  
  nstop <- 0
  lambdamin <- 0
  lambdamax <- lambdae
  
  resl <- res2
  
  repeat{
    lambdae <- lambdamin+(lambdamax-lambdamin)*abs((Emax-Emismin)/(Emismax-Emismin))
    lambda <- c(lambdae,0)
    
    init <- list()
    for (j in 1:ntir) {
      i2_i <- runif(nbpoints,0,1)
      i2_i[1:numdeb] <- res_ref$par[1:numdeb]
      i2_i[Tmax+1:numdeb] <- res_ref$par[Tmax+1:numdeb]
      i2_i[Tmax+(numdeb+10):nmax] <- is[2]
      i2_i[Tmax+nmax+1:numdeb] <- res_ref$par[Tmax+nmax+1:numdeb]
      init[[j]] <- i2_i
    }
    init[[ntir+1]] <- res2$par
    init[[ntir+2]] <- resl$par
    
    resl <- estim_robust(init,numdeb,lambda,Elimit,seqt,alphan,kn,En,em,a,Tmax,delta,is,rho)
    resl$Elimit <- Elimit
    resl$lambda <- c(lambdae,0)
    resl$initpar <- res_ref$initpar
    resl$parameters <- res_ref$parameters
    
    tabl <- format_result(resl)
    Emis <- max(tabl$Em_stock)
    if (Emis<Emax) {
      lambdamax <- lambdae
      Emismax <- Emis
    }
    if (Emis>=Emax) {
      lambdamin <- lambdae
      Emismin <- Emis
    }
    print(c(nstop,lambdae,Emis))
    if ((abs(lambdamax-lambdamin)<1e-4) & (Emax > Emis)) 
      break
    #if (abs(Emax-Emis)<1e-4) {break}
    nstop <- nstop+1
    n_max[[date_BC_i]] <- nstop+1
    if (nstop>=20) 
      break
    
  }
  saveRDS(resl,fl)
}

# On ne garde que 2023 ici
for (Emax in c(4, 4.5, seq(5, 6, by = 0.1))) {
  param_nb <- which(list_parameters$Emax==Emax)
  saveRDS(readRDS(paste("est_mod/simul_rev_",param_nb,"_var_Emax_S4.RDS",sep="")), 
          sprintf("est_mod/ZEN_BC%s_Emax_%s_S4.RDS",2023,Emax))
}

########################################################
####### Cibles imposées tous les 2, 5, ou 10 ans #######
########################################################


for (Emax in c(3.93, 6.23, 5.5, 4.5, seq(5, 6, by = 0.1))) {
  print(sprintf("Emax = %s", Emax))
  if (Emax == 3.93) {
    f0 <- "est_mod/simul_rev_1_baseline_S1.RDS"
    f1 <- "est_mod/simul_rev_1_baseline_S4.RDS"
  } else {
    param_nb <- which(list_parameters$Emax==Emax)
    f0 <- "est_mod/simul_rev_1_baseline_S1.RDS"
    f1 <- paste("est_mod/simul_rev_",param_nb,"_var_Emax_S4.RDS",sep="")
  }
  res0 <- readRDS(f0)
  res1 <- readRDS(f1)
  tab1 <- format_result(res1)
  
  parameters <- res1$parameters
  for (name in colnames(parameters)){
    assign(name, parameters[[name]])
  }
  
  alpha0 <- 1
  beta0 <-1
  nmax <- length(seqt)
  
  Tdeb <- as.numeric(str_split_fixed(parameters$seqt, ":",2)[1])
  Tfin <- as.numeric(str_split_fixed(parameters$seqt, ":",2)[2])
  seqt <- Tdeb:Tfin
  nmax <- length(seqt) # utile lorsque l'on prend un pas de temps qui augmente
  nbpoints <- 2*parameters$Tmax+nmax # plus d'investissement ni de dépreciation du brun après Tmax et pas de dépréciation du vert
  delta <- c(delta_1 = parameters$delta_1, delta_2 = parameters$delta_2) # taux de dépreciation du capital
  e_t_l_max <- parameters$E * (1-0.38) # cible : baisse de 55 % des émissions en 2030 p/r à 1990
  
  initpar <- res1$initpar
  k0 <- initpar$k0
  i0 <- initpar$i0
  E0 <- initpar$E0
  c0 <- initpar$c0
  em <- initpar$em
  a <- initpar$a
  is <- initpar$is
  
  a[6] <- parameters$Cout_Ajust
  
  # Construction des cibles d'émissions intermédiaires fondées sur la trajectoire budget carbone
  Emit_BC <- tab1$Em_flux1b-epuits
  numdeb <- 1
  
  ntir <- 10
  lnint <- c(2,5,10)
  numint <- 1
  for (numint in (1:length(lnint))) {
    set.seed(151071)
    nint <- lnint[numint]
    ndeb <- 1
    struct <- c(1:ndeb)
    count <- ndeb
    for (j in (1:(round(nmax/nint)+1))) {
      struct <- c(struct,rep(count+1,nint))
      count <- count+nint
    }
    struct <- struct[1:(nmax+1)]
    Elimit_nint <- Emit_BC[struct]  # Plafonds d'émission
    Elimit_nint[(Tmax+2):(nmax+1)] <- 0 # Contrainte ZEN
    lambda <- c(0,0)
    
    init <- list()
    for (j in 1:ntir) {
      i2_i <- runif(nbpoints,0,1)
      i2_i[1:(numdeb-1)] <- res0$par[1:(numdeb-1)]
      i2_i[Tmax+1:(numdeb-1)] <- res0$par[Tmax+1:(numdeb-1)]
      i2_i[Tmax+((numdeb+10):nmax)] <- is[2]
      i2_i[Tmax+nmax+1:(numdeb-1)] <- res0$par[Tmax+nmax+1:(numdeb-1)]
      init[[j]] <- i2_i
    }
    init[[ntir+1]] <- res1$par
    
    res2 <- estim_robust(init,1,lambda,Elimit_nint,seqt,alpha0,k0,E0,em,a,Tmax,delta,is,rho)
    res2$Elimit <- Elimit_nint
    res2$lambda <- lambda
    res2$initpar <- res1$initpar
    res2$parameters <- res1$parameters
    
    f2=sprintf("est_mod/ZEN_Interm%s_Emax_%s_S4.RDS",nint,Emax)
    saveRDS(res2,f2)
  }
}
