library(dplyr)
library(patchwork)
library(gt)
library(ggplot2)
library(stringr)
source("R/0_functions.R")
points_in_plots <- FALSE

# parameters
list_parameters <- read.csv("R/parameters.csv", encoding = "utf-8")
# On enlève les variantes pour sigma < 1
list_parameters <- list_parameters[list_parameters$sigma > 1,]
variantes <- read.csv("R/variantes.csv", encoding = "utf8")

#########################################################
##### Chargement fichiers

liste_fichiers <- list.files("est_mod", "^ZEN_")
ZEN <- lapply(liste_fichiers, function(f){
  voir <- format_result(readRDS(file.path("est_mod", f)))
  voir <- result2ts(voir, start =list_parameters$andeb[1])
  data.frame(date = time(voir), voir)
})
names(ZEN) <- gsub(".RDS", "", liste_fichiers)


S1 <- S2 <- S2b <- S3 <-  S4 <- sens <- list()
# construct dataset for figures
for (parameters_name in unique(list_parameters[!list_parameters$param_name %in% c("baseline"), "param_name"])) {
  # parameters_name = "var_Emax"
  print(parameters_name)
  loop <- 1
  for (parameters_nb in unique(list_parameters[list_parameters$param_name %in% c("baseline", parameters_name), "param_nb"])) {
    #parameters_nb = 2
    parameters <- list_parameters[list_parameters$param_nb == parameters_nb,]
    
    if (parameters$param_name == "baseline") { # to add the baseline scenario in the figures
      f1 = sprintf("est_mod/simul_rev_%s_%s_S1.RDS", parameters_nb, "baseline")
      f2 = sprintf("est_mod/simul_rev_%s_%s_S2.RDS", parameters_nb, "baseline")
      f2b = sprintf("est_mod/simul_rev_%s_%s_S2b.RDS", parameters_nb, "baseline")
      f3 = sprintf("est_mod/simul_rev_%s_%s_S3.RDS", parameters_nb, "baseline")
      f4 = sprintf("est_mod/simul_rev_%s_%s_S4.RDS", parameters_nb, "baseline")  
    } else{
      f1 = sprintf("est_mod/simul_rev_%s_%s_S1.RDS", parameters_nb, parameters_name)
      f2 = sprintf("est_mod/simul_rev_%s_%s_S2.RDS", parameters_nb, parameters_name)
      f2b = sprintf("est_mod/simul_rev_%s_%s_S2b.RDS", parameters_nb, parameters_name)
      f3 = sprintf("est_mod/simul_rev_%s_%s_S3.RDS", parameters_nb, parameters_name)
      f4 = sprintf("est_mod/simul_rev_%s_%s_S4.RDS", parameters_nb, parameters_name)
    }
    if (!file.exists(f1)) 
      next
    
    for (name in colnames(list_parameters)){
      assign(name, parameters[[name]])
    }
    
    variable <- variantes[variantes$param_name == parameters_name, "variable"] # variable selon laquelle on étudie une robustesse
    nom_simul <- paste0(param_name, "_nb", param_nb)
    S1[[nom_simul]] <- format_result(readRDS(f1))
    S2[[nom_simul]] <- format_result(readRDS(f2))
    S2b[[nom_simul]] <- format_result(readRDS(f2b))
    S3[[nom_simul]] <- format_result(readRDS(f3))
    S4[[nom_simul]] <- format_result(readRDS(f4)) 
    
    if (loop == 1) {
      tab <- cbind(scenario = 1, name_var = variable, var = get(variable), parameters_nb = parameters_nb, param_name = parameters_name, table_res(S1[[nom_simul]]))
    } else {
      tab <- rbind(tab, cbind(scenario = 1, name_var = variable, var = get(variable), parameters_nb = parameters_nb, param_name = parameters_name, table_res(S1[[nom_simul]])))
    }
    loop <- 2
    
    tab <- rbind(
      tab, 
      cbind(scenario = 2, name_var = variable, var = get(variable), parameters_nb = parameters_nb, param_name = parameters_name, table_res(S2[[nom_simul]])), 
      cbind(scenario = "2b", name_var = variable, var = get(variable), parameters_nb = parameters_nb, param_name = parameters_name, table_res(S2b[[nom_simul]])),
      cbind(scenario = 3, name_var = variable, var = get(variable), parameters_nb = parameters_nb, param_name = parameters_name, table_res(S3[[nom_simul]])),
      cbind(scenario = 4, name_var = variable, var = get(variable), parameters_nb = parameters_nb, param_name = parameters_name, table_res(S4[[nom_simul]]))
    )
  }
  sens[[parameters_name]] <- as.data.frame(tab)
}
S1 <- lapply(S1, function(x) {
  x <- result2ts(x)
  data.frame(date = time(x), x)
})
S2 <- lapply(S2, function(x) {
  x <- result2ts(x)
  data.frame(date = time(x), x)
})
S2b <- lapply(S2b, function(x) {
  x <- result2ts(x)
  data.frame(date = time(x), x)
})
S3 <- lapply(S3, function(x) {
  x <- result2ts(x)
  data.frame(date = time(x), x)
})

S4 <- lapply(S4, function(x) {
  x <- result2ts(x)
  data.frame(date = time(x), x)
})

baseline <- list(S1 = S4$baseline_nb1, # budget carbone
                 S2 = S1$baseline_nb1, # Zen 2050
                 S3 = S2$baseline_nb1, # fit for 55 cappé
                 S3b = S2b$baseline_nb1, # fit for 55 non-cappé
                 S4 = S3$baseline_nb1 # Fit for 55 + Fit for 90 cappé
) 


p_bc <- pres_scenario(baseline[["S1"]])
p_zen <- pres_scenario(baseline[["S2"]])
p_fitfor55cappe <- pres_scenario(baseline[["S3"]])
p_fitfor55noncappe <- pres_scenario(baseline[["S3b"]])
p_fitfor5590cappe <- pres_scenario(baseline[["S4"]])
p_bc

saveRDS(p_bc,
        "DT/images/graph/fiche_bc.RDS")
saveRDS(p_zen,
        "DT/images/graph/fiche_zen.RDS")
saveRDS(p_fitfor55cappe,
        "DT/images/graph/fiche_fitfor55.RDS")
saveRDS(p_fitfor55noncappe,
        "DT/images/graph/fiche_fitfor55noncappe.RDS")
saveRDS(p_fitfor5590cappe,
        "DT/images/graph/fiche_fitfor5590.RDS")

#####################################################
# Graphe de comparaison des émissions par scénario ##
#####################################################
emissions <- sapply(c(baseline[c(1:3,5)] #On ne prend pas le non-cappé
), `[`, , "Em_flux1b")

size <- 1.2

emissions <- ts(emissions, start = baseline$S1[1,1])
emissions <- window(emissions, end = 2053)
em_0 <- ts(emissions[1,1], start = start(emissions),
           end = end(emissions))
emissions <- ts.union(emissions, em_0)
colnames(emissions) <- c("Budget carbone\n1.6 °C",
                         "ZEN", "Fit for 55", 
                         "Fit for 55 + 90",
                         "Sans contrainte")
bc <- c(2030, window(emissions, start = 2030, end = 2030)[,1])
zen <- c(2045, window(emissions, start = 2045, end = 2045)[,2])
fit55 <- c(2045, window(emissions, start = 2045, end = 2045)[,3])
fit5590 <- c(2045, window(emissions, start = 2045, end = 2045)[,4])
em_0_lab <- c(2045, window(emissions, start = 2045, end = 2045)[,5])

data_label <- data.frame(rbind(bc,
                               zen,
                               fit55,
                               fit5590,
                               em_0_lab),
                         label = colnames(emissions),
                         hjust = c(1, 0.5, 0.9, 0.7, 0.5),
                         vjust = c(1, -1, 1.5, 1.2, 1))
colnames(data_label)[1:2] <- c("x", "y")
p <- gg_plot(emissions, outDec = ".",linewidth = size,
             add_points = points_in_plots) +
  geom_text(aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust, color = label),
            data = data_label,inherit.aes = FALSE) +
  scale_color_grey() +
  scale_x_continuous(breaks = c(2022, seq(2025, 2050, by = 5))) +
  scale_y_continuous(breaks = seq(0.0, 0.4, by = 0.05)) +
  theme(legend.position="none",
        panel.grid.minor = element_blank()) +
  labs(y = latex2exp::TeX("GtCO${}_2$éq"))
p
saveRDS(p, file = "DT/images/graph/comparaison_emissions.RDS")

budget_carbone_temperature <- read.csv(
  "DT/data/budget_carbone_temperature.csv",
  header=TRUE, row.names = 1,
  nrows = 2, check.names = FALSE)
budget_carbone_temperature <- data.frame(
  t(budget_carbone_temperature), temp = colnames(budget_carbone_temperature),
  x = 2023
)
bc_tp <- budget_carbone_temperature
emissions_st <- sapply(c(baseline[c(1:3,5)] #On ne prend pas le non-cappé
), `[`, , "Em_stock")

size <- 1.2

emissions_st <- ts(emissions_st, start = baseline$S1[1,1])
emissions_st <- window(emissions_st, start = 2023, end = 2053)
em_0 <- ts(cumsum(c(0, rep(0.368771294121981, length(time(emissions_st))-1))),
           start = start(emissions_st))
emissions_st <- ts.union(emissions_st, em_0)
colnames(emissions_st) <- c("Budget carbone\n1.6 °C",
                            "ZEN", "Fit for 55", 
                            "Fit for 55 + 90",
                            "Sans contrainte")
bc <- c(2053, window(emissions_st, start = 2050, end = 2050)[,1])
zen <- c(2053, window(emissions_st, start = 2050, end = 2050)[,2])
fit55 <- c(2053, window(emissions_st, start = 2050, end = 2050)[,3])
fit5590 <- c(2053, window(emissions_st, start = 2050, end = 2050)[,4])
em_0_lab <- c(2041, window(emissions_st, start = 2041, end = 2041)[,5])

data_label <- data.frame(rbind(bc,
                               zen,
                               fit55,
                               fit5590, 
                               em_0_lab),
                         label = colnames(emissions_st),
                         hjust = c(0, 0, 0, 0, 1.2),
                         vjust = c(0.5, 0.5, 0.5, 0.5, 0.5))
colnames(data_label)[1:2] <- c("x", "y")
data_fleche <- data.frame(x = 2051.5, xend = 2051.5,
                          y = data_label$y[c(2,3,4)],
                          yend = data_label$y[c(3,4,1)])
data_fleche$label <- sprintf("%.f %%",(data_fleche$yend/data_fleche$y - 1)*100)
data_fleche$ylabel <- (data_fleche$y + data_fleche$yend) / 2

p <- gg_plot(emissions_st, outDec = ".",linewidth = size,
             add_points = points_in_plots) +
  geom_text(aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust, color = label),
            data = data_label,
            inherit.aes = FALSE) +
  scale_color_grey() +
  geom_hline(aes(yintercept = France), data = bc_tp, alpha = 0.2,
             linetype =1) +
  geom_text(aes(x = x, y = France, label = temp), data = bc_tp,
            vjust = 1, hjust = 0.4, inherit.aes = FALSE,
            alpha = 0.5) +
  scale_x_continuous(breaks = c(2023, seq(2025, 2050, by = 5)),
                     limits = c(2023, 2059)) +
  scale_y_continuous(breaks = c(0, seq(1, 7, by = 1))) +
  coord_cartesian(ylim = c(0, 7.2)) +
  theme(legend.position="none",
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank()) +
  labs(y = latex2exp::TeX("GtCO${}_2$éq")) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend),
               data = data_fleche,
               arrow = arrow(length = unit(0.3, "cm")),
               linetype = 1,inherit.aes = F,
               lineend = "round",
               linejoin = "bevel",
               linewidth=0.2,
               alpha = 0.7) +
  geom_text(aes(x = x, y = ylabel, label = label),
            data = data_fleche,inherit.aes = FALSE,
            hjust = 1, alpha = 0.6)
p 
saveRDS(p, file = "DT/images/graph/comparaison_emissions_stocks.RDS")

##########################################################
## Cibles régulières d'émissions, Budget Carbone 1.75°C
##########################################################
annees <- c(10, 5, 2)
zen_interm_pib_1.6 <- ZEN[sprintf("ZEN_Interm%s_Emax_3.93_S4", annees)]
names(zen_interm_pib_1.6) <- sprintf("%s ans", annees)
p_zen_interm_pib_1.6 <- fiche_pib_p_serie(zen_interm_pib_1.6)
wrap_plots(p_zen_interm_pib_1.6, ncol = 2) +
  plot_layout(guides = 'collect')  &
  labs(color = "Cibles tous les",
       linetype = "Cibles tous les")

zen_interm_pib_1.75 <- ZEN[sprintf("ZEN_Interm%s_Emax_6.23_S4", annees)]
names(zen_interm_pib_1.75) <- sprintf("%s ans", annees)
p_zen_interm_pib_1.75 <- fiche_pib_p_serie(zen_interm_pib_1.75)
wrap_plots(p_zen_interm_pib_1.75, ncol = 2) +
  plot_layout(guides = 'collect')&
  labs(color = "Cibles tous les",
       linetype = "Cibles tous les")

saveRDS(list(z1.6 = p_zen_interm_pib_1.6, z1.75=p_zen_interm_pib_1.75),
        "DT/images/graph/zen_interm_pib.RDS")


##########################################################
## Budget carbone 1.6°C avec différentes dates de démarrage de la transition
##########################################################
annees <- c(2023, 2028, 2033)
zen_bc_pib1.6 <- lapply(annees, function(y){
  ZEN[[sprintf("ZEN_BC%s_Emax_3.93_S4", y)]]
})
names(zen_bc_pib1.6) <- annees
zen_bc_pib1.75 <- lapply(annees, function(y){
  ZEN[[sprintf("ZEN_BC%s_Emax_6.23_S4", y)]]
})
names(zen_bc_pib1.75) <- annees

zen_bc_pib1.6 <- fiche_pib_p_serie(zen_bc_pib1.6)
zen_bc_pib1.75 <- fiche_pib_p_serie(zen_bc_pib1.75)

wrap_plots(zen_bc_pib1.6, ncol = 2) +
  plot_layout(guides = 'collect') &
  labs(color = "Démarrage en",
       linetype = "Démarrage en")

saveRDS(list(z1.6 = zen_bc_pib1.6, z1.75=zen_bc_pib1.75),
        "DT/images/graph/zen_bc_pib.RDS")


## Flux émissions
annees <- c(2023, 2028, 2033)
em_flux <- do.call(cbind, lapply(annees, function(y){
  ZEN[[sprintf("ZEN_BC%s_Emax_3.93_S4", y)]][,"Em_flux1b"]
}))
colnames(em_flux) <- annees
em_flux <- ts(em_flux, start = ZEN[[sprintf("ZEN_BC%s_Emax_3.93_S4", annees[1])]][1,"date"])
em_flux <- window(em_flux, end = 2053)
dataGraph <- data.frame(as.numeric(time(em_flux)), em_flux)
colnames(dataGraph) <- c("date", colnames(em_flux))

dataGraph <- reshape2::melt(dataGraph, id="date",
                            variable.name = "Démarrage en")
dataGraph$x <- as.numeric(dataGraph$value)
d2023 <- c(2025, window(em_flux, start = 2025, end = 2025)[,1])
d2028 <- c(2032, window(em_flux, start = 2032, end = 2032)[,2])
d2033 <- c(2032, window(em_flux, start = 2032, end = 2032)[,3])

data_label <- data.frame(rbind(d2023,
                               d2028,
                               d2033),
                         label = sprintf("Démarrage en\n%s", annees),
                         hjust = c(1, 1, 0),
                         vjust = c(1, 1, -0.5))
colnames(data_label)[1:2] <- c("x", "y")

p <- ggplot(data = dataGraph, aes(x = date, y = value, group = `Démarrage en`)) +
  geom_line() +
  geom_point(data = dataGraph[as.character(dataGraph$`Démarrage en`) > 2023,] , aes(shape = `Démarrage en`)) +
  scale_shape_manual(values=c(16, 0)) +
  theme_bw() +
  labs(y = latex2exp::TeX("GtCO${}_2$éq"),
       x = NULL) +
  scale_x_continuous(breaks = c(seq(2023, 2050, by = 5), 2050))  +
  geom_text(aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
            data = data_label,inherit.aes = FALSE,
            size = 3) +
  theme(legend.position="none",
        panel.grid.minor = element_blank())
p  
saveRDS(p, file = "DT/images/graph/em_flux_bcdiff.RDS")

###########################################
# Construction des analyses de sensibilité #
###########################################

size <- 1

all_var = c(
  "cuminv1_pib" = "Cumul investissements bruns", 
  "cuminv2_pib" = "Cumul investissements verts",
  "last_inv1" = "Date du dernier investissement brun", 
  "pechv1_pib" = "Cumul du capital brun échoué",
  "ave_conso_gap" = "Baisse moyenne de la consommation"
)
graph_sens_Emax <- lapply(names(all_var), function(yvar){
  d_graph <- sens[["var_Emax"]][, c("var", "scenario", yvar)]
  d_graph <- reshape2::melt(d_graph, id.vars = c("var", "scenario"))
  d_graph <- arrange(d_graph, scenario)
  d_graph$value <- as.numeric(d_graph$value)
  d_graph$var <- as.numeric(d_graph$var)
  d_graph <- d_graph[d_graph$scenario == "4",]
  
  p <- ggplot(data = d_graph, aes(x = var, y = value)) + 
    geom_line(linewidth=size) +
    theme_bw() + #ca enleve le fond gris
    labs(title=all_var[yvar], x=latex2exp::TeX("$E_{max}$"), y="")
  # + 
  #   scale_colour_manual(values = scales::hue_pal()(3)[3],
  #                       labels = c("Budget carbone")) 
  if (points_in_plots) {
    p <- p + geom_point(aes(shape = scenario))
  }
  p +
    guides (color="none", shape = "none")
})
wrap_plots(graph_sens_Emax, ncol = 2)

size <- 1.2

label_xvar = c("rho" = "$\\rho$", "delta_1" = "$\\delta$",
               "epuits" = "Puits carbone", "part_K_b" = "Part $K_b$",
               "sigma" = "$\\sigma$")
ind_sens <- c("var_rho", "var_delta", "var_puits", "var_part_K_b", "var_sigma")
graph_sens_autres <- lapply(ind_sens, function(ind){
  res <- lapply(names(all_var), function(yvar){
    d_graph <- sens[[ind]][, c("var", "scenario", yvar)]
    d_graph <- d_graph[d_graph$scenario %in% c("1", "2", "3", "4"), ]
    variable <- sens[[ind]][1, "name_var"]
    variable <- label_xvar[variable]
    d_graph <- reshape2::melt(d_graph, id.vars = c("var", "scenario"))
    d_graph <- arrange(d_graph, scenario)
    d_graph$value <- as.numeric(d_graph$value)
    d_graph$var <- as.numeric(d_graph$var)
    labels = c("4" = "Budget carbone\n1.6 °C",
               "1" = "ZEN",
               "2" = "Fit for 55",
               "3" = "Fit for 55 + 90")
    d_graph$scenario <- factor(labels[d_graph$scenario],
                               levels = labels,
                               ordered = TRUE)
    d_graph$scenario <- factor(d_graph$scenario,
                               levels = labels,
                               ordered = TRUE)
    
    # # version grise : 
    g <- ggplot(data = d_graph, aes(x = var, y = value, linetype=scenario, color = scenario)) +
      geom_line(linewidth = size) +
      theme_bw() +
      labs(title=all_var[yvar], x=latex2exp::TeX(variable), y="") +
      # scale_linetype_manual(values= c("dotted", "dashed", "solid"),
      #                       labels = labels) +
      scale_color_grey()+
      labs(colour = "Scénario",
           linetype = "Scénario",
           shape = "Scénario")
    if (points_in_plots) {
      g <- g +
        geom_point(aes(shape = scenario)) +
        scale_shape_ordinal(labels = labels)
    }
    
    g 
  })
  # Graphique avec la legende
  res <- c(res, list(guide_area()))
  res 
})
names(graph_sens_autres) <- ind_sens
all_graph_sens <- c(list(var_Emax = graph_sens_Emax),
                    graph_sens_autres)
wrap_plots(all_graph_sens$var_rho, ncol = 2) +
  plot_layout(guides = 'collect')
wrap_plots(all_graph_sens$var_delta, ncol = 2) +
  plot_layout(guides = 'collect')
(wrap_plots(all_graph_sens$var_puits, ncol = 2) +
    plot_layout(guides = 'collect')) &
  theme(plot.title = element_text(size = 12,
                                  face="bold"))
(wrap_plots(all_graph_sens$var_sigma, ncol = 2) +
    plot_layout(guides = 'collect')) &
  theme(plot.title = element_text(size = 12,
                                  face="bold"))
saveRDS(all_graph_sens, file = "DT/images/graph/analyse_sens.RDS")

####################################
##### Graph couts d'ajustement ##### 
####################################
data_ca <- list("Nuls" = S4$baseline_nb1,
                "Modérés" = S4$var_Cout_ajust_nb66,
                "Élevés" = S4$var_Cout_ajust_nb68)
p_pib_ca <- fiche_pib_p_serie(data_ca)
wrap_plots(p_pib_ca , ncol = 2)  +
  plot_layout(guides = 'collect') &
  labs(color = "Coûts d'échouage",
       linetype = "Coûts d'échouage") 
saveRDS(p_pib_ca, file = "DT/images/graph/pib_ca.RDS")

#####################
##### Tableaux ######
#####################

#Construction table baseline

name_rowtab <- c(
  "Cumul des émissions en 2050 (Gt CO2 eq.)",
  "Écart utilité intertemporelle p/r scénario sans contrainte (%)",
  "Capital brun",
  "Capital vert",
  "Investissement brun",
  "Investissement vert",
  "Investissement total",
  "Échouage brun",
  "Consommation en niveau",
  "Perte en Mds",
  "Perte en %",
  "PIB en niveau",
  "Perte en Mds",
  "Perte en %"
)

var_act <- function(Sequence,rho){
  nmax <- length(Sequence)-1
  res <- Sequence
  if (nmax>1){
    res <-  (
      sum(Sequence[2:(nmax)]/(1+rho)^(-2+2:(nmax)))+Sequence[nmax+1]/rho/(1+rho)^(nmax-2)
    )*rho/(1+rho)*1000
  }
  res
}

constr_table <- function(S,name){
  tabi <- as.matrix(c(S$Em_stock[S$date==2050],
                      S$crit[2],
                      var_act(S$capital1b,rho),
                      var_act(S$capital2,rho),
                      var_act(S$inv1b,rho),
                      var_act(S$inv2,rho),
                      var_act(S$inv1b + S$inv2,rho),
                      var_act(S$echv1,rho),
                      var_act(S$conso,rho),
                      0,
                      0,
                      var_act(S$prod,rho),
                      0,
                      0))
  colnames(tabi) <- name
  tabi
}

tab0 <- as.matrix(
  c((S1$baseline_nb1$Em_flux1b[1]-epuits)*(Tmax+1),
    S1$baseline_nb1$crit[1],
    S1$baseline_nb1$capital1b[1]*1000,
    S1$baseline_nb1$capital2[1]*1000,
    S1$baseline_nb1$inv1b[1]*1000,
    S1$baseline_nb1$inv2[1]*1000,
    (S1$baseline_nb1$inv1b[1] + S1$baseline_nb1$inv2[1])*1000,
    S1$baseline_nb1$echv1[1]*1000,
    S1$baseline_nb1$conso[1]*1000,
    0,
    0,
    S1$baseline_nb1$prod[1]*1000,
    0,
    0))
colnames(tab0) <- "Sans contrainte carbone"

tabbaseline <- cbind(
  tab0,
  constr_table(S4$baseline_nb1, "Budget Carbone 1.6°C et ZEN"),
  constr_table(ZEN$ZEN_BC2023_Emax_6.23_S4, "Budget Carbone 1.75°C et ZEN"),
  constr_table(S1$baseline_nb1, "ZEN"),
  constr_table(S2$baseline_nb1, "Fit for 55"),
  constr_table(S3$baseline_nb1, "Fit for 55 + 90"))

Ajust_table <- function(tabbaseline,tab0,name_rowtab){
  rownames(tabbaseline) <- name_rowtab
  # Ecart utilité
  tabbaseline[2,] <- (tabbaseline[2,]-tab0[2])/tab0[2] * 100
  tabbaseline[10,] <- tabbaseline[9,]-tab0[9]
  tabbaseline[11,] <- tabbaseline[10,]/tab0[9]*100
  tabbaseline[13,] <- tabbaseline[12,]-tab0[12]
  tabbaseline[14,] <- tabbaseline[13,]/tab0[12]*100
  # tabbaseline <- round(tabbaseline,digits=2)
  tabbaseline
}

tabbaseline <- Ajust_table(tabbaseline,tab0,name_rowtab)

tabInterm1.6 <- cbind(
  constr_table(ZEN$ZEN_BC2023_Emax_3.93_S4,"Cibles Annuelles"),
  constr_table(ZEN$ZEN_Interm2_Emax_3.93_S4,"Cibles Tous les 2 ans"),
  constr_table(ZEN$ZEN_Interm5_Emax_3.93_S4,"Cibles Tous les 5 ans"),
  constr_table(ZEN$ZEN_Interm10_Emax_3.93_S4,"Cibles Tous les 10 ans")
)

tabInterm1.6 <- Ajust_table(tabInterm1.6,tab0,name_rowtab)

tabInterm1.75 <- cbind(
  constr_table(ZEN$ZEN_BC2023_Emax_6.23_S4,"Cibles Annuelles"),
  constr_table(ZEN$ZEN_Interm2_Emax_6.23_S4,"Cibles Tous les 2 ans"),
  constr_table(ZEN$ZEN_Interm5_Emax_6.23_S4,"Cibles Tous les 5 ans"),
  constr_table(ZEN$ZEN_Interm10_Emax_6.23_S4,"Cibles Tous les 10 ans")
)

tabInterm1.75 <- Ajust_table(tabInterm1.75,tab0,name_rowtab)

tabRetard1.6 <- cbind(
  constr_table(ZEN$ZEN_BC2023_Emax_3.93_S4,"Budget Carbone 2023"),
  constr_table(ZEN$ZEN_BC2028_Emax_3.93_S4,"Budget Carbone 2028"),
  constr_table(ZEN$ZEN_BC2033_Emax_3.93_S4,"Budget Carbone 2033")
)

tabRetard1.6 <- Ajust_table(tabRetard1.6,tab0,name_rowtab)

tabRetard1.75 <- cbind(
  constr_table(ZEN$ZEN_BC2023_Emax_6.23_S4,"Budget Carbone 2023"),
  constr_table(ZEN$ZEN_BC2028_Emax_6.23_S4,"Budget Carbone 2028"),
  constr_table(ZEN$ZEN_BC2033_Emax_6.23_S4,"Budget Carbone 2033"),
  constr_table(ZEN$ZEN_BC2038_Emax_6.23_S4,"Budget Carbone 2038")
)

tabRetard1.75 <- Ajust_table(tabRetard1.75,tab0,name_rowtab)
all_table <- list(tabbaseline=tabbaseline,
                  tabInterm1.6=tabInterm1.6, tabInterm1.75=tabInterm1.75,
                  tabRetard1.6=tabRetard1.6, tabRetard1.75=tabRetard1.75)
saveRDS(all_table,
        "DT/images/tableaux.RDS")

all_table <- readRDS("DT/images/tableaux.RDS")
# Baseline :
data <- data.frame(rowname = rownames(all_table$tabbaseline), all_table$tabbaseline)
colnames(data)[-1] <- colnames(all_table$tabbaseline)
# colnames(data)[-1] <- gsub("Budget Carbone", "", colnames(data)[-1],
#                            ignore.case = TRUE)
data$group_row <- c(rep("", 8), 
                    rep("Consommation", 3),
                    rep("PIB", 3))
data$rowname <- 
  gsub("\\w+ en niveau", "Niveau", data$rowname)
data |> 
  gt(groupname_col = "group_row") |> 
  fmt_number(decimals = 2, sep_mark = " ", rows = c(1,2,11,14)) |> 
  fmt_number(decimals = 0, sep_mark = " ", rows = c(3, 4, 5, 6, 7, 8, 9, 10, 12, 13))
# tab_spanner(
#   label = "Budget carbone",
#   columns = 3:4
# )

# tabInterm
data <- data.frame(rowname = rownames(all_table$tabInterm1.6), all_table$tabInterm1.6)
colnames(data)[-1] <- gsub("Cibles", "", colnames(all_table$tabInterm1.6),
                           ignore.case = TRUE)
data$group_row <- c(rep("", 8), 
                    rep("Consommation", 3),
                    rep("PIB", 3))
data$rowname <- 
  gsub("\\w+ en niveau", "Niveau", data$rowname)
data |> 
  gt(groupname_col = "group_row") |> 
  fmt_number(decimals = 1, sep_mark = " ") |>
  tab_spanner(
    label = "Cibles",
    columns = 2:5
  )

# tabRetard :
data <- data.frame(rowname = rownames(all_table$tabRetard1.6), all_table$tabRetard1.6)
colnames(data)[-1] <- gsub("Budget Carbone", "", colnames(all_table$tabRetard1.6),
                           ignore.case = TRUE)
data$group_row <- c(rep("", 8), 
                    rep("Consommation", 3),
                    rep("PIB", 3))
data$rowname <- 
  gsub("\\w+ en niveau", "Niveau", data$rowname)
data |> 
  gt(groupname_col = "group_row") |> 
  fmt_number(decimals = 1, sep_mark = " ") |>
  tab_spanner(
    label = "Budget carbone",
    columns = 2:5
  )
