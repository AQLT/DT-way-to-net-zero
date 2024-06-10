#' ---
#' title: "Statistiques utilisées dans l'article"
#' author: "Abbas R., Carnot N., Lequien M., Quartier-la-Tente A. et Roux S."
#' date: "Mai 2024"
#' output: html_document
#' ---

# Rapport obtenu en utilisant la fonction rmarkdown::render("R/4_statistiques_article.R")
if (basename(getwd()) == "R"){
  # Si l'espace de travail est le sous dossier "R" on se replace à la racine du projet.
  path <- normalizePath("../")
  setwd(path)
  knitr::opts_knit$set(root.dir = path)
}

source("R/0_functions.R")
library(dplyr, warn.conflicts = FALSE)
library(stringr)
library(ggplot2)

list_parameters <- read.csv("R/parameters.csv", encoding = "utf-8")
andeb <- 2022

# Statistiques mobilisées dans l’article ----

p_fitfor55 <- readRDS("DT/images/graph/fiche_fitfor55.RDS")
p_fitfor55noncappe <- readRDS("DT/images/graph/fiche_fitfor55noncappe.RDS")
p_bc <- readRDS("DT/images/graph/fiche_bc.RDS")
p_zen <- readRDS("DT/images/graph/fiche_zen.RDS")
p_fitfor5590 <- readRDS("DT/images/graph/fiche_fitfor5590.RDS")

## bac à sable pour regarder les données: Fit for 55 cappé et non cappé
# p_fitfor55noncappe[[1]]$data
# p_fitfor55noncappe[[2]][[1]]$data
# p_fitfor55noncappe[[2]][[2]]$data
# 
# p_fitfor55[[1]]$data
# p_fitfor55[[2]][[1]]$data
# p_fitfor55[[2]][[2]]$data

## Fit for 55 ----
## émissions supplémentaires dues au fait de ne pas capper fit for 55 après 2030:
p_fitfor55noncappe[[2]][[2]]$data |>
  filter(date == 2050, variable == "Émissions\n(stock)") |> 
  pull(value) - 
  p_fitfor55[[2]][[2]]$data |>
  filter(date == 2050, variable == "Émissions\n(stock)") |> 
  pull(value)

## investissement brun (hors résiduel) nul dès 2023
p_fitfor55[[1]]$data |>
  filter(date <= 2025, variable == "Investissement brun") 

## investissement brun et vert pour fit for 55 et fit for 55 non cappé autour de 2030
p_fitfor55noncappe[[1]]$data |>
  filter(2029 <= date & date <= 2034, variable %in% c("Investissement brun", "Investissement vert")) 
p_fitfor55[[1]]$data |>
  filter(2029 <= date & date <= 2034, variable %in% c("Investissement brun", "Investissement vert")) 

## Évolution du stock de capital brun en 2031
p_fitfor55noncappe[[2]][[1]]$data|>
  filter(2030 <= date & date <= 2033, variable == "Capital brun") |> 
  pull(value) 
p_fitfor55[[2]][[1]]$data|>
  filter(2030 <= date & date <= 2033, variable == "Capital brun") |> 
  pull(value) 

## Échouage: nul sauf années de la contrainte pour ZEN, FF55 et FF5590, quelle part pour BC en 2023 et 2050 ----
p_zen[[2]][[1]]$data |>
  filter(2023 <= date & date < 2050, variable == "Capital brun\néchoué") |> 
  pull(value) |> 
  sum()

p_fitfor55[[2]][[1]]$data |>
  filter(2023 <= date & date < 2050 & date!=2030, variable == "Capital brun\néchoué") |> 
  pull(value) |> 
  sum()

p_fitfor5590[[2]][[1]]$data |>
  filter(2023 <= date & date < 2050 & date!=2030 & date != 2040, variable == "Capital brun\néchoué") |> 
  pull(value) |> 
  sum()

p_fitfor5590[[2]][[1]]$data |>
  filter(date %in% c(2030, 2040, 2050), variable == "Capital brun\néchoué")

## Part du capital brun échoué en 2023 et 2050 avec BC
round(
  (p_bc[[2]][[1]]$data |> 
     filter(date %in% c(2023, 2050), variable == "Capital brun\néchoué") |> 
     pull(value) /
     p_bc[[2]][[1]]$data |> 
     filter(date == 2022, variable == "Capital brun") |> 
     pull(value)) * 100)

## scénario de sensibilité Emax
### Note: on peut aussi regarder, dans 4_figures_articles.R, graph_sens_Emax[[1]]$data
### Emax > 6.34, égalité entre scénarios BC et ZEN
### investissement brun (hors résiduel) nul en dessous de 5.5Gt
liste_var <- c(var_fiche_pib, var_capital, var_emissions)
scen_zen <- format_result(readRDS("est_mod/simul_rev_1_baseline_S1.RDS"))

result <- readRDS(sprintf("est_mod/simul_rev_%s_%s_S4.RDS", 1, "baseline"))
epuits <- result$parameters$epuits
rho <- result$parameters$rho
voir <- format_result(result)
res <- table_res(voir)
d <- data.frame(Emax = list_parameters[list_parameters$param_nb == 1, "Emax"], 
                cuminv1_pib = res[, "cuminv1_pib"], 
                last_inv1 = res[, "last_inv1"],
                egalite_ZEN = isTRUE(all.equal(voir[liste_var],
                                               scen_zen[liste_var])))

param_name <- "var_Emax"
param_nb <- list_parameters[list_parameters$param_name == param_name, "param_nb"]
for (nb in param_nb){
  #nb <- 86
  result <- readRDS(sprintf("est_mod/simul_rev_%s_%s_S4.RDS", nb, param_name))
  epuits <- result$parameters$epuits
  rho <- result$parameters$rho
  voir <- format_result(result)
  res <- table_res(voir)
  d <- rbind(d, c(list_parameters[list_parameters$param_nb == nb, "Emax"],
                  res[, "cuminv1_pib"], 
                  res[, "last_inv1"],
                  isTRUE(all.equal(voir[liste_var],
                                   scen_zen[liste_var]))))
}
d <- d[order(d$Emax),]
d

## investissement total en 2049 supérieur à l’I initial total pour ZEN et FF55, pas les autres [égalité pour BC]
data <- p_fitfor55[[1]]$data
data <- p_fitfor55noncappe[[1]]$data
data <- p_bc[[1]]$data
data <- p_zen[[1]]$data
data <- p_fitfor5590[[1]]$data

(data[data$date==2049 & data$variable == "Investissement vert", "value"] + 
data[data$date==2049 & data$variable == "Investissement brun", "value"]) / 
(data[data$date==2022 & data$variable == "Investissement vert", "value"] + 
data[data$date==2022 & data$variable == "Investissement brun", "value"])

(data[data$date==2047 & data$variable == "Investissement vert", "value"] + 
 data[data$date==2047 & data$variable == "Investissement brun", "value"]) / 
(data[data$date==2022 & data$variable == "Investissement vert", "value"] +
 data[data$date==2022 & data$variable == "Investissement brun", "value"])

## date de dépassement du BC 1.6 = 3.93 Gt ----
### ZEN 2036
d <- p_zen
data <- d[[2]][[2]]$data
data[data$variable == "Émissions\n(stock)" & inrange(data$date, 2035, 2036), ]
cumul_em_zen <- data[data$variable == "Émissions\n(stock)" & data$date ==2050, "value"]
cumul_em_zen
### FF55 2038
d <- p_fitfor55
data <- d[[2]][[2]]$data
data[data$variable == "Émissions\n(stock)" & inrange(data$date, 2037, 2039), ]
cumul_em_fitfor55 <- data[data$variable == "Émissions\n(stock)" & data$date ==2050, "value"]
cumul_em_fitfor55
### FF5590 2039
d <- p_fitfor5590
data <- d[[2]][[2]]$data
data[data$variable == "Émissions\n(stock)" & inrange(data$date, 2038, 2040), ]
cumul_em_fitfor5590 <- data[data$variable == "Émissions\n(stock)" & data$date ==2050, "value"]
cumul_em_fitfor5590

### sans contrainte carbone 2033
d <- p_zen
data <- d[[2]][[2]]$data
#### vérif émissions (nettes du puits de carbone) cumulées entre 2023 et 2050 inclus:
emissions_annuelles_initales <- 
  data[data$date == 2022 & data$variable ==  "Émissions\n(flux)", "value"] - 
  data[data$date == 2051 & data$variable ==  "Émissions\n(flux)", "value"]
cumul_em_sanscontrainte <- (2050-2023 + 1) * emissions_annuelles_initales
(2032-2023 + 1) * emissions_annuelles_initales
(2033-2023 + 1) * emissions_annuelles_initales

### BC
d <- p_bc
data <- d[[2]][[2]]$data
cumul_em_bc <- data[data$variable == "Émissions\n(stock)" & data$date ==2050, "value"]

## comparaison des cumuls d’émissions en 2050 dans les différents scénarios ----
### ZEN vs sans contrainte -0.39
cumul_em_zen / cumul_em_sanscontrainte - 1
### FF55 vs ZEN -0.12
cumul_em_fitfor55 / cumul_em_zen - 1
### FF5590 vs FF55 -0.19
cumul_em_fitfor5590 / cumul_em_fitfor55 - 1
### BC vs FF5590 -0.13
cumul_em_bc / cumul_em_fitfor5590 - 1

## investissement brun FF55 et FF5590 et BC et ZEN ----
d <- p_fitfor5590
d <- p_fitfor55
d <- p_bc
d <- p_zen
d[[1]]$data


## Cibles tous les 10, 5, 2 ans et BC, selon le budget carbone alloué (Emax) ----
Emax <- 5.9
#### BC
param_nb <- list_parameters[list_parameters$Emax == Emax, "param_nb"]
param_name <- list_parameters[list_parameters$Emax == Emax, "param_name"]
f <- sprintf("est_mod/simul_rev_%s_%s_S4.RDS", param_nb, param_name)
voirBC <- format_result(readRDS(f))

#### cibles ponctuelles
nint <- 10
f = sprintf("est_mod/ZEN_Interm%s_Emax_%s_S4.RDS", nint, Emax)
res <- readRDS(f)
voir <- format_result(res)
all.equal(voir, voirBC)
cumul_echouage_2022_2049 <- sum(voir[voir$t <= 2050-2023, "echv1"])

### BC 3.93
Emax <- 3.93
param_nb <- 1
param_name <- list_parameters[list_parameters$param_nb == param_nb, "param_name"]
f <- sprintf("est_mod/simul_rev_%s_%s_S4.RDS", param_nb, param_name)
voirBC <- format_result(readRDS(f))
#### echouage la première et dernière année
voirBC[voirBC$t %in% c(1, 28), "echv1"]

#### écart de cumul d’émissions entre cible à 10, 5 et 2 ans et BC
nint <- 10
f = sprintf("est_mod/ZEN_Interm%s_Emax_%s_S4.RDS", nint, Emax)
res <- readRDS(f)
voir <- format_result(res)
# Em cumulées
round(max(voir[, "Em_stock"]), 2)
voir[voir$t == 28, "Em_stock"]
voir[voir$t == 28, "Em_stock"] - voirBC[voirBC$t == 28, "Em_stock"]
voir[voir$t == 28, "Em_stock"] / voirBC[voirBC$t == 28, "Em_stock"] - 1

nint <- 5
f = sprintf("est_mod/ZEN_Interm%s_Emax_%s_S4.RDS", nint, Emax)
res <- readRDS(f)
voir <- format_result(res)
round(max(voir[, "Em_stock"]), 2)
voir[voir$t == 28, "Em_stock"]
voir[voir$t == 28, "Em_stock"] - voirBC[voirBC$t == 28, "Em_stock"]
voir[voir$t == 28, "Em_stock"] / voirBC[voirBC$t == 28, "Em_stock"] - 1

nint <- 2
f = sprintf("est_mod/ZEN_Interm%s_Emax_%s_S4.RDS", nint, Emax)
res <- readRDS(f)
voir <- format_result(res)
round(max(voir[, "Em_stock"]), 2)
voir[voir$t == 28, "Em_stock"] - voirBC[voirBC$t == 28, "Em_stock"]
voir[voir$t == 28, "Em_stock"] / voirBC[voirBC$t == 28, "Em_stock"] - 1

# 5.4 cumul d'émissions
zen_interm_5.4 <- format_result(readRDS(sprintf("est_mod/ZEN_Interm%s_Emax_%s_S4.RDS", 5, 5.4)))
zen_bc_5.4 <- format_result(readRDS(sprintf("est_mod/ZEN_BC2023_Emax_%s_S4.RDS", 5.4)))

zen_interm_5.5 <- format_result(readRDS(sprintf("est_mod/ZEN_Interm%s_Emax_%s_S4.RDS", 5, 5.5)))
zen_bc_5.5 <- format_result(readRDS(sprintf("est_mod/ZEN_BC2023_Emax_%s_S4.RDS", 5.5)))

zen_bc_5.4[1:30, c("t", "echv1")]
zen_bc_5.5[1:30, c("t", "echv1")]
zen_interm_5.4[1:30, c("t", "echv1")]
zen_interm_5.5[1:30, c("t", "echv1")]

zen_interm_5.4[1:30, c("t", "inv1")]
zen_interm_5.5[1:30, c("t", "inv1")]

## BC retardé ----
#### baisse des émissions dans le scénario BC 1.6°C, démarrage en 2033
bc_2023 <- format_result(readRDS(sprintf("est_mod/ZEN_BC%s_Emax_3.93_S4.RDS", "2023")))
bc_2028 <- format_result(readRDS(sprintf("est_mod/ZEN_BC%s_Emax_3.93_S4.RDS", "2028")))
bc_2033 <- format_result(readRDS(sprintf("est_mod/ZEN_BC%s_Emax_3.93_S4.RDS", "2033")))
bc_2033[bc_2033$t == 2033-2022, "Em_flux1b"] / bc_2033[bc_2033$t == 2032-2022, "Em_flux1b"] - 1
bc_2033[c("t", "Em_flux1b")] 

### BC 3.93
#### echouage la première année
bc_2023[bc_2023$t== (2023 - 2022), "echv1"]
bc_2028[bc_2028$t== (2028 - 2022), "echv1"]
bc_2033[bc_2033$t== (2033 - 2022), "echv1"]
bc_2033[bc_2033$t== (2050 - 2022), "echv1"]


## cout d’ajustement ----
# scenarios n° 66 et 68 (85 trop extrême) 
param_nb <- 66
param_name <- list_parameters[list_parameters$param_nb == param_nb, "param_name"]
a66 <- c(0, 0, 0, 0, 0, list_parameters[list_parameters$param_nb == param_nb, "Cout_Ajust"])
f <- sprintf("est_mod/simul_rev_%s_%s_S4.RDS", param_nb, param_name)
voir66 <- format_result(readRDS(f))
voir66[voir66$t %in% c(1, 28), "echv1"]
voir66[voir66$t %in% c(1, 28), "echv1"] / voirBC[voirBC$t %in% c(1, 28), "echv1"]
voir66$Em_stock

param_nb <- 68
param_name <- list_parameters[list_parameters$param_nb == param_nb, "param_name"]
a68 <- c(0, 0, 0, 0, 0, list_parameters[list_parameters$param_nb == param_nb, "Cout_Ajust"])
f <- sprintf("est_mod/simul_rev_%s_%s_S4.RDS", param_nb, param_name)
voir68 <- format_result(readRDS(f))
voir68[voir68$t %in% c(1, 28), "echv1"]
voir68[voir68$t %in% c(1, 28), "echv1"] / voirBC[voirBC$t %in% c(1, 28), "echv1"]
1/ voir68[voir68$t %in% c(1, 28), "echv1"] * voirBC[voirBC$t %in% c(1, 28), "echv1"]
voir68$Em_stock

param_nb <- 85
param_name <- list_parameters[list_parameters$param_nb == param_nb, "param_name"]
a85 <- c(0, 0, 0, 0, 0, list_parameters[list_parameters$param_nb == param_nb, "Cout_Ajust"])
f <- sprintf("est_mod/simul_rev_%s_%s_S4.RDS", param_nb, param_name)
voir85 <- format_result(readRDS(f))
voir85[voir85$t %in% c(1, 28), "echv1"]
voir85[voir85$t %in% c(1, 28), "echv1"] / voirBC[voirBC$t %in% c(1, 28), "echv1"]
voir85$Em_stock

d_graph <- rbind(cbind(voir66[, c("t", "Em_flux1b")], scen = 66), 
                 cbind(voir85[, c("t", "Em_flux1b")], scen = 85),
                 cbind(voir68[, c("t", "Em_flux1b")], scen = 68))
g <- ggplot(d_graph[d_graph$t <= 31, ], aes(x = t, y = Em_flux1b, group=scen, color = scen)) +
  geom_line()
g


cout_aj(voir85[voir85$t == 1, "echv1"], a85)
cout_aj(voir68[voir68$t == 1, "echv1"], a68)
cout_aj(voir66[voir66$t == 1, "echv1"], a66)
util(voir85$conso)
util(voir68$conso)
util(voir66$conso)

# cout_aj(voir85[, "echv1"], a85) / util(voir85$conso) * 100
cout_aj(voir66[, "echv1"], a66) / util(voir66$conso) * 100

# coût d’ajustement, en % de l’utilité initiale
# cout_aj(voir85[, "echv1"], a85) / util(voir85[1 , "conso"]) * 100
cout_aj(voir68[, "echv1"], a68) / util(voir68[1 , "conso"]) * 100 # 1% en 2023

cout_aj(voir66[, "echv1"], a66) / util(voir66[1 , "conso"]) * 100 # 0.1% en 2023



### Émissions cumulées scénario BC
d <- p_bc
data <- d[[2]][[2]]$data
cumsum(data[data$variable == "Émissions\n(flux)" & data$date > andeb, "value"] - 0.035)
