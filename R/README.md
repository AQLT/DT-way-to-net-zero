## Description des programmes

- `R/0_functions.R` : fonctions utilisées dans les autres scripts

- `R/1_estimations.R` : estimation des modèles et analyse de sensibilité (selon paramètres définis dans le fichier `parameters.csv`) qui sont ensuite sauvegardés dans le dossier `est_mod`.

- `R/2_variantes_BC.R` : estimation de 2 variantes au budget carbone : démarrage retardé (2028, 2033 ou 2038) et avec des cibles espacées dans le temps (2, 5 ou 10 ans).
Les modèles sont sauvegardés dans le dossier `est_mod`.

- `R/3_figures_article.R` : création de l'ensemble des graphiques du document de travail

- `R/4_statistiques_articles.R` : calcul des statistiques utilisées dans l'article, il est utilisé pour générer un fichier html lisible à partir du lien <https://inseefrlab.github.io/DT-way-to-net-zero/R/4_statistiques_article.html>.
