#functional assumptions : u(c)=ln(c),
#CES F(k)=((a[1]*(kminb+k[1]))^((sigma-1)/sigma)+(a[2]*k[2])^((sigma-1)/sigma)^(alpha*sigma/(sigma-1))
fprod <- function(k, a, delta) {
  kminb <- a[3]
  alpha <- a[4]
  sigma <- a[5]
  as.numeric(
    (a[1] * (kminb + k[1])) ^ ((sigma - 1) / sigma) + (a[2] * k[2]) ^
      ((sigma - 1) / sigma)
  ) ^ (alpha * sigma / (sigma - 1)) - delta[1] * kminb
}


dfprod1 <- function(k, a) {
  kminb <- a[3]
  alpha <- a[4]
  sigma <- a[5]
  unlist(
    alpha * (a[1]) ^ ((sigma - 1) / sigma) * (kminb + k[1]) ^ (-1 / sigma) *
      ((a[1] * (kminb + k[1])) ^ ((sigma - 1) / sigma) + (a[2] * k[2]) ^ 
         ((sigma - 1) / sigma)) ^ (alpha * sigma / (sigma - 1) - 1)
  )
}

dfprod2 <- function(k, a) {
  kminb <- a[3]
  alpha <- a[4]
  sigma <- a[5]
  unlist(
    alpha * (a[2]) ^ ((sigma - 1) / sigma) * (k[2]) ^ (-1 / sigma) *
      ((a[1] * (kminb + k[1])) ^ ((sigma - 1) / sigma) + (a[2] * k[2]) ^ 
         ((sigma - 1) / sigma)) ^ (alpha * sigma / (sigma - 1) - 1)
  )
}


util <- function(c)
{
  log(c)
}

dutil <- function(c)
{
  1 / c
}

d2util <- function(c)
{
  -1 / c ^ 2
}

cout_aj <- function(c, a) {
  if (length(a) < 6) {
    int <- 0
  }
  if (length(a) >= 6) {
    int <- a[6]
  }
  as.numeric(int * c ^ 2)
}

dcout_aj <- function(c, a) {
  if (length(a) < 6) {
    int <- 0
  }
  if (length(a) >= 6) {
    int <- a[6]
  }
  as.numeric(2 * int * c)
}



#définition récursive de la valeur intertemporelle
#rajout de la contrainte sur les émissions, au-delà de la période max, seul k2 est actif

objectif <- function(i,lambda,Elimit, seqt,alphai,kt,Et,numseq, em, a, Tmax, delta, is, rho){
  nmax <- length(seqt)
  Phit <- i[numseq,c(3,4)]
  kminb <- a[3]
  alpha <- a[4]
  kte <- kt*(1-Phit) #actifs échoués
  #Contrainte à chaque date sur les émissions nettes. Donc Elimit=0 au-delà d'une certaine date
  if (Elimit[1]>0) {
    kte[1] <- min(kte[1],Elimit[seqt[numseq]+1]/em[1])
  }
  kte[kte<=0] <- 0
  
  et <- (kte[1])*em[1] #émissions nettes, le puits de carbone est déjà comptabilisé
  fpr <- fprod(kte, a, delta)
  if (numseq==nmax) {
    #on utilise la solution obtenue sous hypothèse stationnaire pour la dernière valeur du capital obtenue
    ct <- fpr-delta[2]*kte[2]-delta[1]*kte[1]
    alphaf <- alphai*(1+rho*(seqt[numseq]-seqt[numseq-1]))
    if (ct<=10^(-10)) {ct <- 10^(-10)}
    v <- alphaf*1/rho*util(ct)
  }
  if (numseq<nmax) {
    Deltan <- seqt[numseq+1]-seqt[numseq]
    It <- cbind(i[numseq,1],i[numseq,2]*(1-i[numseq,1]))*(fpr)
    ct <- (1-i[numseq,2])*(1-i[numseq,1])*fpr
    ktp <- (kte+Deltan*(-delta*kte+It))
    Etp <- Et+et
    alphap <- alphai/(1+rho*Deltan) # alpha utilisé à la date suivante
    if (ct <= 10^(-10)) { ct <- 10^(-10) }
    v <- (alphai*Deltan*util(ct)-cout_aj(kt[1]-kte[1],a)- 
            #(lambda[1]+lambda[2]*(seqt[numseq]==tl))*et +
            lambda[1]*et +
            objectif(i,lambda,Elimit, seqt,alphap,ktp,Etp,numseq+1,
                     em = em, a = a, Tmax = Tmax, delta = delta, is = is, rho = rho))
  }
  v
}


#i2: séquence des variables de contrôle
#numdeb: date de début d'optimisation, numdeb=1, on part de 2020

crit <- function(i2,numdeb,lambda,Elimit,seqt, alpha0, k0, E0, em, a, Tmax, delta, is, rho) {
  # simplification du problème
  nmax <- length(seqt)
  i <- matrix(rep(0,4*nmax),nrow=nmax,ncol=4)
  i[1:Tmax,1] <- i2[1:Tmax]
  i[1:nmax,2] <- i2[Tmax+1:nmax]
  i[1:Tmax,3] <- i2[Tmax+nmax+1:Tmax]
  -objectif(i,lambda,Elimit,seqt,alpha0,k0,E0,numdeb, 
            em = em, a = a, Tmax = Tmax, delta = delta, is = is, rho = rho)
}

#définition récursive du gradient

#on définit dkt qui est la dérivée du capital kt1 et kt2 par rapport à toutes les variables de contrôles
#dkt est de taille ncol=2 nrow=4*nmax

dobjectif <- function(i,lambda,Elimit, seqt,alphai,kt,dkt,Et,numseq, em, a, Tmax, delta, is, rho){
  nmax <- length(seqt)
  Phit <- t(as.matrix(i[numseq,c(3,4)]))
  kminb <- a[3]
  alpha <- a[4]
  kte <- kt*(1-Phit) #capital efficace
  
  dkte <- dkt*((1-Phit) %x% as.matrix(rep(1,4*nmax)))
  dkte[2*nmax+numseq,1] <- unlist(dkte[2*nmax+numseq,1] - kt[1])
  dkte[3*nmax+numseq,2] <- unlist(dkte[3*nmax+numseq,2] - kt[2])
  #Contrainte à chaque date sur les émissions nettes. Donc Elimit=0 au-delà d'une certaine date
  if (Elimit[1]>0) {
    if (kte[1]>Elimit[seqt[numseq]+1]/em[1]) {
      kte[1] <- Elimit[seqt[numseq]+1]/em[1]
      dkte[,1] <- rep(0,nrow(dkte)) #car le niveau de capital est alors constant, limité par la borne supérieure des émissions
    }
  }
  kte[kte<=0] <- 0
  #On passe de t-1 à t
  et <- (kte[1])*em[1] #émissions nettes, le puits de carbone est déjà comptabilisé
  det <- as.matrix(dkte[,1]*em[1])
  #fpr <- fprod(kt[1]+kt[2]*a[2]/a[1],a[1]^alpha)
  fpr <- fprod(kte, a, delta)
  dfpr <- dkte %*% rbind(dfprod1(kte,a),dfprod2(kte,a))
  if (numseq==nmax) {
    #on utilise la solution obtenue sous hypothèse stationnaire pour la dernière valeur du capital obtenue
    ct <- fpr- as.matrix(kte) %*% as.matrix(delta)
    dct <- (dfpr - dkte %*% as.matrix(delta))
    alphaf <- alphai*(1+rho*(seqt[numseq]-seqt[numseq-1])) 
    if (ct<=10^(-10)) {ct <- 10^(-10)}
    v <- alphaf*1/rho*util(ct)
    dv <- drop(alphaf/rho*dutil(ct)) * dct 
  }
  if (numseq<nmax) {
    Deltan <- seqt[numseq+1]-seqt[numseq]
    It <- cbind(i[numseq,1],i[numseq,2]*(1-i[numseq,1]))*(fpr)
    dIt<- (cbind(i[numseq,1],i[numseq,2]*(1-i[numseq,1])) %x% rep(1,4*nmax))* (cbind(1,1) %x% dfpr)
    dIt[numseq,1] <- dIt[numseq,1]+fpr
    dIt[numseq,2] <- dIt[numseq,2]-i[numseq,2]*fpr
    dIt[nmax+numseq,2] <-dIt[nmax+numseq,2]+(1-i[numseq,1])*fpr
    ct <- (1-i[numseq,2])*(1-i[numseq,1])*fpr
    dct <- (1-i[numseq,2])*(1-i[numseq,1])*dfpr
    dct[numseq] <- dct[numseq]-(1-i[numseq,2])*fpr
    dct[nmax+numseq] <- dct[nmax+numseq]-(1-i[numseq,1])*fpr
    ktp <- (1-Deltan*delta)*kte+Deltan*It
    dktp <- (1-Deltan*delta)*dkte+Deltan*dIt
    Etp <- Et+et
    alphap <- alphai/(1+rho*Deltan) # alpha utilisé à la date suivante
    if (ct <= 10^(-10)) { ct <- 10^(-10) }
    v <- (alphai*Deltan*util(ct) -cout_aj(kt[1]-kte[1],a)
          #          - (lambda[1]+lambda[2]*(seqt[numseq]==tl))*et+
          - lambda[1]*et+
            objectif(i,lambda,Elimit, seqt,alphap,ktp,Etp,numseq+1, 
                     em = em, a = a, Tmax = Tmax, delta = delta, is = is, rho = rho))
    dvp <- dobjectif(i,lambda, Elimit, seqt,alphap,ktp,dktp, Etp, numseq+1, em, a, Tmax, delta, is, rho)
    dv <- alphai*Deltan*dutil(ct)*dct-dcout_aj(kt[1]-kte[1],a)*(dkt[,1]-dkte[,1])-
      lambda[1]*det+dvp
    #      (lambda[1]+lambda[2]*(seqt[numseq]==tl))*det+dvp
    
  }                                                 
  #if (numseq==numdeb) {print(v)}
  # print(c(v,ct,fpr,kte))
  # print(dv)
  dv
}

dcrit <- function(i2,numdeb,lambda,Elimit,seqt, alpha0, k0, E0, em, a, Tmax, delta, is, rho) {
  # simplification du problème
  nmax <- length(seqt)
  i <- matrix(rep(0,4*nmax),nrow=nmax,ncol=4)
  i[1:Tmax,1] <- i2[1:Tmax]
  i[1:nmax,2] <- i2[Tmax+1:nmax]
  i[1:Tmax,3] <- i2[Tmax+nmax+1:Tmax]
  dk0 <- matrix(rep(0,2*4*nmax),nrow=4*nmax, ncol=2)
  dv <- -dobjectif(i,lambda,Elimit,seqt,alpha0,k0,dk0,E0,numdeb, 
                   em = em, a = a, Tmax = Tmax, delta = delta, is = is, rho = rho)
  dvfin <- dv[c((1:Tmax),(nmax+(1:nmax)),(2*nmax+(1:Tmax)))]
  dvfin
}

estim_robust <- function(init,numdeb,lambda,Elimit,seqt,alpha0,k0,E0,em,a,Tmax,delta,is,rho, print = TRUE){
  ninit <- length(init)
  res <- list()
  for (j in 1:ninit) {
    res[[j]] <- future({  
      
      res2 <- optim(init[[j]],crit,numdeb=numdeb,lambda=lambda, Elimit=Elimit, seqt = seqt,
                    alpha0 = alpha0, k0 = k0, E0 = E0,
                    em = em, a = a, Tmax = Tmax, delta = delta, is = is, rho = rho,
                    lower=rep(0,4*nmax),upper=rep(1,4*nmax),method="L-BFGS-B", gr=dcrit)
      print(c(j,res2$value))
      res2
    })
  }
  res <- lapply(res, future::value)
  valres <- rep(0,ninit)  
  for (j in 1:(ninit)) {
    valres[j] <- res[[j]]$value
  }
  if (print)
    print(c(which.min(valres),min(valres)))
  res2 <- res[[which.min(valres)]]
  res2$Elimit <- Elimit
  res2$lambda <- lambda
  res2
}


#' Tracer des graphiques avec dygraphs
#' 
#' @param data un objet `ts`
#' @param titre,sous_titre titre et sous_titre
#' @param legende si l'on veut changer le nom des variables sur le graphique
#' @param xlab,ylab label des axes.
#' @param outDec séparateur décimal.
#' @param digits nombre de décimale.
#' @param color couleurs des séries.
#' @param date_window un éventuel vecteur avec date de début et de fin de la fenêtre temporelle sélectionnée dans l'affichage au départ.
#' @param ylimits,color_ylimits vecteur pour ajouter des un trait horizontal et la couleur associée. Si le vecteur est nommé, les noms sont utilisés pour la légende.
#' 

dy_plot <- function(data, titre = NULL, sous_titre = NULL,
                    legende = NULL,
                    xlab = NULL, ylab = NULL,
                    outDec = ".",
                    digits = 2,
                    color = NULL,
                    date_window = NULL,
                    ylimits = NULL,
                    color_ylimits = "red",
                    yaxis_limits = NULL,
                    ...){
  if (!require(dygraphs)){
    install.packages("dygraphs")
    require(dygraphs)
  }
  if (!is.ts(data))
    stop("Il faut que la table en entrée soit de type ts !")
  if (!is.null(legende) & is.mts(data)){
    colnames(data) <- legende
  }
  
  if (!is.null(titre)) {
    if (!is.null(sous_titre)) {
      titre <- sprintf("%s <br> <small><small>%s</small></small>", titre, sous_titre)
    }
  }
  
  
  dy <- dygraph(data,
                main = titre,
                xlab = xlab, ylab = ylab) %>% 
    dyOptions(colors = color) %>%
    dyRangeSelector(dateWindow = date_window)
  
  if( !is.null(ylimits)) {
    data_col_limit <- data.frame(ylimits, 
                                 names(ylimits),
                                 color_ylimits,
                                 stringsAsFactors = FALSE)
    for (i in seq_len(nrow(data_col_limit))) {
      dy <- dy %>% 
        dyLimit(as.numeric(data_col_limit[i, 1]), 
                label = data_col_limit[i, 2],
                color = data_col_limit[i, 3])
    }
    if(max(ylimits) > max(data)) {
      dy <- dy %>% 
        dyAxis("y", valueRange = c(0, max(data, ylimits)+ 0.5))
      
    }
  }
  
  
  if(!is.null(digits)){
    valueFormatter = sprintf("function formatValue(v) {  
  precision = Math.pow(10, %i)
  return (Math.ceil(v * precision) / precision).toFixed(%i).toString().replace('.','%s') 
  }",digits + 1,
                             digits,
                             outDec)
    dy <- dy %>%
      dyOptions(digitsAfterDecimal = digits) %>% 
      dyAxis("y",axisLabelFormatter = htmlwidgets::JS(valueFormatter),
             valueFormatter = htmlwidgets::JS(valueFormatter))
  }
  if (!is.null(yaxis_limits)) {
    dy <- dy %>% 
      dyAxis("y", valueRange = yaxis_limits)
  }
  dy
}

#' Tracer des graphiques avec ggplot
#' 
#' @param data donnée en sortie de [format_result()]
#' @param titre,sous_titre titre et sous_titre
#' @param legende si l'on veut changer le nom des variables sur le graphique
#' @param xlab,ylab label des axes.
#' @param outDec séparateur décimal.
#' @param digits nombre de décimale.
#' @param color couleurs des séries.
#' @param ylimits,color_ylimits vecteur pour ajouter des un trait horizontal et la couleur associée. Si le vecteur est nommé, les noms sont utilisés pour la légende.
#' @param size épaisseur des traits sur les graphiques
#' @param add_points booléen pour ajouter des points à chaque année sur les graphiques
#' 
gg_plot <- function(data, titre = NULL, sous_titre = NULL,
                    legende = NULL,
                    xlab = NULL, ylab = NULL,
                    outDec = ".",
                    digits = 2,
                    ylimits = NULL,
                    color_ylimits = "red",
                    yaxis_limits = NULL,
                    size=1,
                    add_points = FALSE,
                    ...){
  if (!require(ggplot2)){
    install.packages("ggplot2")
    require(ggplot2)
  }
  time <- time(data)
  dataGraph <- data.frame(cbind(time, data))
  if (is.null(legende)){
    if(is.mts(data)){
      legende <- colnames(data)
    }else{
      legende <- ""
    }
  }
  colnames(dataGraph) <- c("date", legende)
  
  dataGraph <- reshape2::melt(dataGraph, id="date")  # convert to long format
  
  p <- ggplot(data = dataGraph, aes(x = date, y = value, group = variable, colour = variable,
                                    linetype = variable)) +
    geom_line(linewidth = size) +
    labs(title = titre, subtitle = sous_titre,
         x = xlab, y = ylab) +
    scale_y_continuous(labels = function(x) format(x, decimal.mark = outDec))+
    theme_bw()
  if (add_points) {
    p <- p +
      geom_point(aes(shape = variable))
  }
  if( !is.null(ylimits)) {
    data_y_lim <- data.frame(yintercept = ylimits,
                             x = start(data)[1],
                             labels = names(ylimits))
    p + geom_hline(mapping = aes(yintercept = yintercept),
                   data = data_y_lim,
                   color = color_ylimits,
                   linetype = "dashed") +
      geom_text(aes(x = x, y = yintercept, label = labels),
                data = data_y_lim,
                vjust = -1,
                col = "red",
                inherit.aes = FALSE)
  }
  p 
}

result2ts <- function(x, start = andeb) {
  if (colnames(x)[[1]] == "date") {
    ts(x[,-1], start = x[1,1], frequency = 1)
  } else {
    res <- do.call(rbind, lapply(seq(x[1,"t"], x[nrow(x),"t"], by = 1), 
                                 function(i) {
                                   indices <- i- x[,"t"]
                                   indices[indices < 0] <- max(indices)
                                   
                                   x[which.min(indices),]
                                 }))
    ts(res[,-1], start = start, frequency = 1)
  }
  
}

# Combinaison des fonctions précédentes pour tracer des graphiques
graph_results <- function(x,
                          var =  c("prod", "conso", "capital2", "capital1", "Em_stock", "Em_flux1", "ech2", "ech1", "inv2", "inv1", "CAcorr"),
                          date_window,
                          end = NULL,
                          sous_titre = NULL,
                          ...,
                          digits = 2,
                          html = getOption("html.output"),
                          size = 1.2) {
  data <- result2ts(x)
  data <- window(data, end = end)
  data <- data[, var, drop = FALSE]
  colnames(data) <- gsub("1", "_brun", colnames(data))
  colnames(data) <- gsub("2", "_vert", colnames(data))
  if (!is.null(names(var))) {
    colnames(data)[names(var) != ""] <- names(var)[names(var) != ""]
  }
  if (missing(date_window))
    date_window <- c(sprintf("%s-01-01", start(data)[1]),
                     "2053-01-01")
  if (is.null(html) || html) {
    dy_plot(data,
            digits = digits,
            sous_titre = sous_titre,
            date_window = date_window,
            ...)
    
  } else {
    gg_plot(data,
            digits = digits,
            sous_titre = sous_titre,
            date_window = date_window,
            size = size,
            ...)
  }
}


# Fonction interne à format_result()
verif_res <- function(i,lambda,Elimit,seqt,alphai,betai,kt,Et,numseq, em, a, Tmax, Emax, delta, is, rho) {
  nmax <- length(seqt)
  Phit <- i[numseq,c(3,4)]
  kminb <- a[3]
  alpha <- a[4]
  kte <- kt*(1-Phit) #actifs échoués
  #Contrainte à chaque date
  if (Elimit[1]>0) {
    kte[1] <- min(kte[1],Elimit[seqt[numseq]+1]/em[1])
  }
  #Contrainte sur le cumul total d'émissions
  kte[1] <- min(kte[1],(Emax-Et)/em[1])
  kte[kte<=0] <- 0
  et <- kte[1]*em[1]
  fpr <- unlist(fprod(kte, a, delta))
  if (numseq==nmax) {
    #on utilise la solution obtenue sous hypothèse stationnaire pour la dernière valeur du capital obtenue
    ct <- fpr-delta[2]*kte[2]-delta[1]*kte[1]
    alphaf <- alphai*(1+rho*(seqt[numseq]-seqt[numseq-1]))
    betaf <- betai*(1+(rho+delta[1])/(1-delta[1])*(seqt[numseq]-seqt[numseq-1]))
    # res <- c(seqt[numseq],fpr,ct,kt,Et,is,Phit)
    CA_inst <- delta[1]/em[1]*(dfprod1(kte,a)-dfprod2(kte,a))
    CA <- betaf*(1+rho)/(rho+delta[1])*CA_inst
    
    if (ct<=10^(-10)) {ct <- 10^(-10)}
    v <- alphaf*1/rho*util(ct) 
    
    res_int <- unlist(c(seqt[numseq],v,fpr,ct,kt,kte[1],Et,is,Phit,CA,CA_inst/delta[1]))
    #colnames(res_int) <- NULL
    res <- t(as.matrix(rep(0,length(res_int))))
    res[1,1:length(res_int)] <- res_int
  }
  if (numseq<nmax) {
    Deltan <- seqt[numseq+1]-seqt[numseq]
    
    It <- cbind(i[numseq,1],i[numseq,2]*(1-i[numseq,1]))*(fpr)
    ct <- (1-i[numseq,2])*(1-i[numseq,1])*(fpr)
    ktp <- (kte+Deltan*(-delta*kte+It))
    
    Etp <- Et+et
    alphap <- alphai/(1+rho*Deltan)
    betap <- betai/(1+(rho+delta[1])/(1-delta[1])*Deltan)
    resp <- verif_res(i,lambda,Elimit,seqt,alphap,betap,ktp,Etp,numseq+1,
                      em, a, Tmax, Emax, delta, is, rho)
    CA_inst <- delta[1]/em[1]*(dfprod1(kte,a)-dfprod2(kte,a))
    CA <- betai*CA_inst+resp[1,ncol(resp)-1] #on prend la valeur correspondant à CA
    
    v <- alphai*Deltan*util(ct)+resp[1,2]
    
    res_int <- unlist(c(seqt[numseq],v,fpr,ct,kt,kte[1],Et,i[numseq,1:2],Phit,CA,CA_inst/delta[1]))
    colnames(res_int) <- NULL
    res <- t(as.matrix(rep(0,length(res_int))))
    res[1,1:length(res_int)] <- res_int
    res_int <- res
    res <- rbind(res_int,resp)
  }
  #colnames(res) <- names_r
  res
}

format_result <- function(res) {
  if (is.list(res))
    lambda <- res$lambda
  Elimit <- res$Elimit
  
  parameters <- res$parameters
  pib <- parameters$pib
  
  Tdeb <- as.numeric(str_split_fixed(parameters$seqt, ":",2)[1])
  Tfin <- as.numeric(str_split_fixed(parameters$seqt, ":",2)[2])
  seqt <- Tdeb:Tfin
  
  Tmax <- parameters$Tmax
  delta <- c(delta_1 = parameters$delta_1, delta_2 = parameters$delta_2)
  rho <- parameters$rho
  
  alpha0 <- 1
  beta0 <-1
  nmax <- length(seqt)
  
  initpar <- res$initpar
  k0 <- initpar$k0
  i0 <- initpar$i0
  E0 <- initpar$E0
  c0 <- initpar$c0
  em <- initpar$em
  a <- initpar$a
  Emax <-  initpar$Emax
  is <- initpar$is
  
  x <- res$par
  
  kminb <- a[3]
  alpha <- a[4]
  
  respar <- matrix(rep(0,4*nmax),nrow=nmax,ncol=4)
  respar[1:Tmax,1] <- x[1:Tmax]
  respar[1:nmax,2] <- x[Tmax+1:nmax]
  respar[1:Tmax,3] <- x[Tmax+nmax+1:Tmax]
  nprod <- 2
  result <- verif_res(respar,lambda, Elimit, seqt,alpha0,beta0,k0,Et = 0, numseq=1, em, a, Tmax, Emax, delta, is, rho)
  
  #valeurs à la période initiale, en état stationnaire
  res_int <- unlist(c(0, util(c0)/rho, pib-delta[1]*kminb, c0, k0, k0[1], 0, i0, 0, 0, 0, 0))
  res <- t(as.matrix(rep(0,length(res_int))))
  res[1,1:length(res_int)] <- res_int
  res_int <- res
  
  colnames(res_int)  <- NULL
  result <- rbind(res_int, result)
  
  result <- data.frame(result)
  colnames(result) <- c("t","crit","prod","conso","capital1","capital2","capital1_e","Em_stock","inv_shr1","inv_shr2","ech_shr1","ech_shr2","CA","CA_inst")
  
  result$prodc <- result$prod+delta[1]*kminb
  result$CAcorr <- result$CA*((1+rho)/(1-delta[1]))^(result$t-1) #Pour comparer à t=1
  
  result$prixcarb <- lambda[1]/dutil(result$conso)
  result$prixcarb_form <-(dfprod1(result[,c("capital1","capital2")],a)-dfprod2(result[,c("capital1","capital2")],a))/em[1]
  
  result$critcorr <- result$crit*(1+rho)^(result$t) #à comparer à t=0
  result[, sprintf("ech%i", 1:nprod)] <- as.matrix(result[, sprintf("ech_shr%i", 1:nprod)] * result[, sprintf("capital%i", 1:nprod)]) %*% diag(1-delta)
  
  # Calcul du capital échoué par solde  
  
  #Capital au début de chaque période
  kp <- result[2:nrow(result),sprintf("capital%i", 1:nprod)]
  echv <- result[1:(nrow(result)-1),sprintf("capital%i", 1:nprod)]-kp
  result[, sprintf("echv%i", 1:nprod)] <- 0
  result[1:(nrow(result)-1), sprintf("echv%i", 1:nprod)] <- echv*(echv>=0)
  
  result[,"inv1"] <- result[,"inv_shr1"] * result[, "prod"]
  result[,"inv2"] <- result[,"inv_shr2"] * result[, "prod"] * (1 - result[,"inv_shr1"])
  
  #Calcul du capital effectivement échoué
  
  result$echv1 <- result$capital1-result$capital1_e
  result$part_echv1 <- result$echv1/result$capital1
  result[is.na(result$part_echv1),"part_echv1"] <- 0
  
  result$capital1b <- result$capital1_e+kminb
  result$inv1b <- result$inv1+delta[1]*kminb
  #Le capital émissif est limité par la borne supérieure des émissions
  
  result[, sprintf("Em_flux%s",c("1b",2))] <- as.matrix(result[, sprintf("capital%s", c("1b",2))]) %*% 
    diag(em)
  
  
  Emaxf <- result$Em_stock[nmax+1]
  Tmaxf <- min(result[result$Em_stock==Emaxf,"t"])
  duc <- dutil(result$conso[nmax+1])
  
  lambda_0 <- lambda[1]
  
  dut <- duc*dfprod1(c(result$capital1_e[nmax+1],result$capital2[nmax+1]),a)
  mut <- dut/em[1]-lambda_0*(1+rho)^(Tmaxf-1)
  LagBC <- lambda_0*(1+rho)^(Tmaxf-1)
  
  lambda_t <- LagBC*rho/(1+rho)
  
  varphit <- 0
  LagKb <- 0
  LagKv <- (1+rho)*duc
  #LagE <- -(1+rho)/rho*mut
  BC <- Emaxf - result$Em_stock[nmax+1]
  contr_mu <- (em[1]*result$capital1_e[nmax+1]>=Elimit[nmax+1])
  contr_lambd <- (BC<=0)
  
  for (l in (nmax:2)){
    keffb <- result$capital1_e[l]
    duc <- dutil(result$conso[l])
    LagKv <- rbind(duc*dfprod2(c(keffb,result$capital2[l]),a)+
                     (1-delta[1])/(1+rho)*LagKv[1],LagKv)
    emt <- em[1]*keffb
    dut_int <- duc*dfprod1(c(keffb,result$capital2[l]),a)
    dut  <- rbind(dut_int,dut)
    
    BC <- Emaxf - result$Em_stock[l]
    contr_mu <- (em[1]*result$capital1_e[l]>=Elimit[l])
    contr_lambd <- (BC<=0)
    
    Interm2 <- (1-delta[1])/(1+rho)*LagKb[1]
    Interm3 <- em[1]/(1+rho)*LagBC[1]
    Interm <- dut[1]+Interm2-Interm3
    
    mutp <- 0
    lambda_tp <- contr_lambd*lambda_t
    if (contr_mu==TRUE) {
      mutp <- Interm/em[1]-lambda_tp
    }
    mut <- rbind(mutp,mut)
    Interm <- Interm-(lambda_tp+mut[1])*em[1]
    
    LagKb <- rbind((1-result$part_echv1[l])*(Interm>0)*Interm,LagKb)
    LagBC <- rbind(lambda_tp+LagBC[1]/(1+rho),LagBC)
  }
  #On rajoute les situations stationnaires de départ
  keffb <- result$capital1_e[1]
  duc <- dutil(result$conso[1])
  result$mut <- rbind(0,mut)
  result$LagKb <- rbind((1+rho)*duc,LagKb)
  result$LagKv <- rbind((1+rho)*duc,LagKv)
  
  result$LagBC <- rbind(0,LagBC)
  result$BC <- Emaxf-result$Em_stock
  
  result$prix_Kb <- result$LagKb/dutil(result$conso)
  result$prix_Kv <- result$LagKv/dutil(result$conso)
  result$prix_BC <- result$LagBC/dutil(result$conso)
  
  result$LagKv_th <- (1+rho)*dutil(result$conso)
  
  result$Wealth_Kb <- result$prix_Kb*result$capital1_e
  result$Wealth_Kv <- result$prix_Kv*result$capital2
  result$Wealth_BC <- result$prix_BC*result$BC
  
  ENA_Kb <- result$prix_Kb[1:nmax]*(result$capital1_e[2:(nmax+1)]-result$capital1_e[1:nmax])
  ENA_Kv <- result$prix_Kv[1:nmax]*(result$capital2[2:(nmax+1)]-result$capital2[1:nmax])
  ENA_BC <- result$prix_BC[1:nmax]*(result$BC[2:(nmax+1)]-result$BC[1:nmax])
  
  result$ENA_Kb <- c(0,ENA_Kb)
  result$ENA_Kv <- c(0,ENA_Kv)
  result$ENA_BC <- c(0,ENA_BC)
  ENA <- (ENA_Kb+ENA_Kv+ENA_BC)
  
  result$ENA <- c(0,ENA)
  result$Pib_adj <- result$conso+result$ENA
  result$Wealth <- (result$Wealth_Kb+result$Wealth_Kv+result$Wealth_BC)
  a_df <- as.data.frame(t(a))
  colnames(a_df) <- paste0("a", 1:length(a))
  result <- merge(result, a_df)
  
  result
}



#Résolution de l'équation de stationnarité
fonc_Kv <-function(kv,a,rho,delta) {
  (log(dfprod2(c(0,kv),a))-log(rho+delta[1]))^2
}

init_param <- function(rho,delta,sigma,alpha0,beta0,pib,K_t,part_K_b,epuits,Emax,E,e_t_l_max){
  K_b <- part_K_b*K_t #Part du brun dans le capital total
  K_v <- K_t-K_b
  
  #alpha calé pour que 2019 soit solution stationnaire exp(log(K_b)-log(pib)+log(rho+delta[1]))
  # on en déduit aussi a_b et a_v, à sigma fixé
  
  alpha <- (rho+delta[1])*K_t/pib
  
  a_b <- (pib)^(1/alpha)*(K_b)^(1/(sigma-1))/(K_t)^(sigma/(sigma-1))
  a_v <- (pib)^(1/alpha)*(K_v)^(1/(sigma-1))/(K_t)^(sigma/(sigma-1))
  
  ##### contraintes environnementales 
  
  
  #Emissions
  E0 <- E-epuits # Emissions nettes à la première date
  # D'après SDES : 522 Mt en 1990, 405 en 2019
  # Contrainte d'émissions nettes
  
  #puits de carbone: 0.035Gt
  #lpuits <- 0.035
  
  e <- E/K_b # Emissivité: on rapporte les émissions totales au capital brun
  em <- c(e,0)
  
  #valeur du capital brun minimal à conserver compensant le puits de carbone
  
  kminb <- epuits/e 
  
  a <- c(a_b,a_v,kminb,alpha,sigma)
  
  K_bi <- K_b-kminb
  
  #ca <- 0.3 # 90 €/tonne = 0,09  , 300€/tonne = 0,3 
  #A l'état stationnaire
  
  #pib <- fprod(c(K_bi,K_v),a, delta)+delta[1]*kminb #pour se rassurer
  Ib0 <- delta[1]*(K_b-kminb)  #Investissement endogène en capital brun
  Iv0 <- delta[2]*K_v
  c0 <- pib-delta[1]*kminb-Ib0-Iv0  #Il faut prévoir l'investissement pour le capital intouchable
  i0 <- c(Ib0 / (pib-delta[1]*kminb), Iv0 / (pib-delta[1]*kminb-Ib0))
  
  # initialisation
  k0 <- cbind(K_b-kminb, K_v)
  
  # Solution stationnaire à l'arrivée
  
  
  #recherche de la valeur de K_v annulant la fonction
  kvmax <- 2*K_v
  repeat{
    if (log(dfprod2(c(0,kvmax),a))-log(rho+delta[1])<0) {break}
    kvmax <- 2*kvmax
  }
  
  result <- optimize(fonc_Kv,a=a,rho=rho,delta=delta, interval = c(K_v, kvmax), tol = 0.0001)
  
  #La solution stationnaire est aussi intéressante à considérer à l'arrivée
  ks2 <- result$minimum
  
  pibs <- fprod(c(0,ks2),a, delta)+delta[1]*kminb
  cs <- pibs-delta[1]*kminb-delta[2]*ks2
  
  is <- c(0,delta[2]*ks2/(pibs-delta[1]*kminb))
  
  list(k0=k0, E0=E0,  i0=i0, c0=c0, em=em, a=a, ks=ks2, cs=cs, is=is)  
}

#Résolution de l'équation de stationnarité
fonc_sigma0 <- function(sigma,alpha,pib,K_b,K_t,kminb,e,CA_fin,rho,delta) {
  K_v <- K_t-K_b
  a_b <- (pib)^(1/alpha)*(K_b)^(1/(sigma-1))/(K_t)^(sigma/(sigma-1))
  a_v <- (pib)^(1/alpha)*(K_v)^(1/(sigma-1))/(K_t)^(sigma/(sigma-1))
  kv <- kminb*(a_v/a_b)^(sigma-1)*(1+e*CA_fin/(rho+delta[1]))^sigma
  a <- c(a_b,a_v,kminb,alpha,sigma)
  list(a=a,kv=kv)
}


table_res <- function(basen,a){
  if (missing(a) || is.null(a)) {
    a <- basen[1, grep("^a\\d+$", colnames(basen))]
  }
  anfin <- andeb+nrow(basen)
  pib2050 <- basen[2050-andeb,"prodc"]
  pib2051 <- basen[2051-andeb,"prodc"]
  pibfin <- basen[nrow(basen),"prodc"]
  conso2050 <- basen[2050-andeb,"conso"]
  conso2051 <- basen[2051-andeb,"conso"]
  consofin <- basen[nrow(basen),"conso"]
  pibdeb <- basen[2,"prodc"]
  consodeb <- basen[2,"conso"]
  conso_init <- basen[1, "conso"]
  diffconso2050_conso_init <- conso2050 / conso_init - 1
  diffconso2020_conso_init <- consodeb / conso_init - 1
  
  diffprod2050 <- dfprod1(basen[2050-andeb,c("capital1","capital2")],a)-
    dfprod2(basen[2050-andeb,c("capital1","capital2")],a)
  diffprod2051 <- dfprod1(basen[2051-andeb,c("capital1","capital2")],a)-
    dfprod2(basen[2050-andeb,c("capital1","capital2")],a)
  diffprodfin <- dfprod1(basen[2100-andeb,c("capital1","capital2")],a)-
    dfprod2(basen[2100-andeb,c("capital1","capital2")],a)
  cumEmission <- sum(basen[1:(2050-andeb),"Em_flux1b"])
  cumEmission_n <- sum(basen[1:(2050-andeb),"Em_flux1b"]-epuits)
  
  vec <- (2:(anfin-andeb))
  cumconso <- sum(basen[vec,"conso"]/((1+rho)^(vec-2)))+
    basen[nrow(basen),"conso"]/rho/(1+rho)^(nrow(basen)-2)
  lcumconso <- sum(log(basen[vec,"conso"])/((1+rho)^(vec-2)))+
    log(basen[nrow(basen),"conso"])/rho/(1+rho)^(nrow(basen)-2)
  ave_conso_gap <- (mean(basen[2:(2050-andeb), "conso"] / conso_init) - 1)*100
  pechv1 <- sum(basen[vec,"echv1"])/basen[1,"capital1b"]
  pechv1_pib <- sum(basen[vec,"echv1"]) / basen[1, "prodc"] # stranded assets, in % of initial GDP
  cuminv1_pib <- sum(basen[2:(2050-andeb), "inv1"]) / basen[1, "prodc"]
  cuminv2_pib <- sum(basen[2:(2050-andeb), "inv2"]) / basen[1, "prodc"]
  last_inv1 <- max(basen[basen$inv1 > 0 & basen$t <= (2050-andeb), "t"]) + andeb
  
  # carbon price and abatement costs
  CAcorr_2020 <- basen[2, "CAcorr"]
  CAcorr_2030 <- basen[12, "CAcorr"]
  prixcarb_2019 <- basen[1, "prixcarb"]
  res <- cbind(pibdeb,pib2050,pib2051,pibfin,consodeb,conso2050,conso2051,consofin, conso_init, diffconso2050_conso_init, diffconso2020_conso_init, ave_conso_gap, cuminv1_pib, last_inv1, cuminv2_pib,
               diffprod2050,diffprod2051,diffprodfin,cumEmission,cumEmission_n,cumconso,lcumconso,pechv1, pechv1_pib, prixcarb_2019, CAcorr_2020, CAcorr_2030)
  
  res
}

############################################################################
### Fonctions pour les graphiques articles (fichier 3_figures_article.R) ### 
############################################################################

var_fiche_pib1 <- c("Production" = "prodc", 
                    "Consommation" = "conso")
var_fiche_pib2 <- c("Investissement vert" = "inv2",
                    "Investissement brun" = "inv1b", 
                    "Capital brun échoué" = "echv1")
var_fiche_pib <- c(var_fiche_pib1,
                   var_fiche_pib2)
titre_pib <- "Fiche de PIB"
var_capital <- c("Capital vert" = "capital2", 
                 "Capital brun" = "capital1b", 
                 "Capital brun\néchoué" = "echv1")
titre_capital <- "Stock de Capital"
var_emissions <- c("Émissions\n(stock)" = "Em_stock", 
                   "Émissions\n(flux)" = "Em_flux1b")
titre_emissions <- "Stock et flux d’émission"

fiche_pib_anc <- function(...){
  p_pib <- graph_results(var = var_fiche_pib,
                         html = FALSE,
                         ...)+
    scale_color_grey()
  p_data <- layer_data(p_pib,1)
  p_pib <- p_pib + 
    geom_point(data = subset(p_data, group == 5),
               aes(x = x, y = y),
               colour =  subset(p_data, group == 5)[1,"colour"],
               inherit.aes = FALSE)
  p_pib
}

fiche_pib <- function(data, titre = NULL, ...){
  p_pib <- graph_results(var = var_fiche_pib,
                         data,
                         titre = NULL,
                         end = 2053,
                         size = 1,
                         add_points = points_in_plots,
                         html = FALSE)+
    scale_color_grey()
  
  p_data <- layer_data(p_pib,1)
  data_col <- unique(p_data[,c("group", "colour", "linetype")])
  p_pib <- p_pib + 
    geom_point(data = subset(p_data, group == 5),
               aes(x = x, y = y),
               colour =  subset(p_data, group == 5)[1,"colour"],
               inherit.aes = FALSE)
  p_pib %+% subset(p_data, group %in% 1:2)
  p_pib1 <- graph_results(var = var_fiche_pib1,
                          data,
                          titre = NULL,
                          end = 2053,
                          size = 1,
                          add_points = points_in_plots,
                          html = FALSE)+
    scale_color_manual(values = data_col$colour[1:2]) +
    scale_linetype_manual(values = data_col$linetype[1:2])
  p_pib2 <- graph_results(var = var_fiche_pib2,
                          data,
                          titre = NULL,
                          end = 2053,
                          size = 1,
                          add_points = points_in_plots,
                          html = FALSE)+
    scale_color_manual(values = data_col$colour[-(1:2)]) +
    scale_linetype_manual(values = data_col$linetype[-(1:2)])
  p_data <- layer_data(p_pib2,1)
  p_pib2 <- p_pib2 + 
    geom_point(data = subset(p_data, group == 3),
               aes(x = x, y = y),
               colour =  subset(p_data, group == 3)[1,"colour"],
               inherit.aes = FALSE)
  p_pib1 + p_pib2 + plot_layout(guides = 'collect') + 
    plot_annotation(title = titre) & (
      theme(legend.title = element_blank(),
            legend.text = element_text(size=8),
            legend.spacing = unit(0, "pt"),
            legend.spacing.x = unit(0, "pt"),
            legend.background = element_blank(), 
            legend.justification = "top",
            legend.box.spacing =  unit(0, "pt"),
            plot.title = element_text(size = 12,
                                      face="bold") 
      )
    )
}

fiche_pib_p_serie <- function(data){
  dates <- data[[1]][["date"]]
  res <- lapply(names(var_fiche_pib), function(yvar){
    data_graph <- ts(sapply(data,`[[`, var_fiche_pib[yvar]),
                     start = dates[1], frequency = 1)
    data_graph <- window(data_graph, end = c(2053))
    
    gg_plot(data_graph,
            digits = 2,
            titre = yvar) +
      scale_color_grey() +
      scale_x_continuous(breaks = c(2022, seq(2030, 2050, by = 10)))
  })
  c(res, list(guide_area()))
}
pres_scenario <- function(data) {
  
  (fiche_pib(data,
             titre = titre_pib) /
     (graph_results(data,
                    var = var_capital,
                    titre = titre_capital,
                    end = 2053,
                    size = 1,
                    add_points = points_in_plots,
                    html = FALSE)+
        scale_color_grey() +
        graph_results(data,
                      var = var_emissions,
                      titre = titre_emissions,
                      end = 2053,
                      size = 1,
                      add_points = points_in_plots,
                      html = FALSE) +
        scale_color_grey()) +  
     plot_annotation(title = titre_pib)
  ) & (theme(legend.title = element_blank(),
             legend.text = element_text(size=8),
             legend.spacing = unit(0, "pt"),
             legend.spacing.x = unit(0, "pt"),
             legend.background = element_blank(), 
             legend.justification = "top",
             legend.box.spacing =  unit(0, "pt"),
             plot.title = element_text(size = 12,
                                       face="bold") 
  )
  ) 
}

inrange <- function(x, lower, upper, incbounds=TRUE){
  if (incbounds) {
    x >= lower & x <= upper
  } else {
    x > lower & x < upper
  }
  
}
