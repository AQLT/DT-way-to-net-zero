# DT-way-to-net-zero

[![Build](https://github.com/inseefrlab/dt-way-to-net-zero/workflows/Dockerize/badge.svg)](https://hub.docker.com/repository/docker/inseefrlab/dt-way-to-net-zero)
[![Onyxia](https://img.shields.io/badge/Launch-Datalab-orange?logo=R)](https://datalab.sspcloud.fr/launcher/ide/rstudio?autoLaunch=false&service.image.custom.enabled=true&service.image.pullPolicy=%C2%ABAlways%C2%BB&service.image.custom.version=%C2%ABinseefrlab%2Fdt-way-to-net-zero%3Alatest%C2%BB&init.personalInit=%C2%ABhttps%3A%2F%2Fraw.githubusercontent.com%2Finseefrlab%2Fdt-way-to-net-zero%2Fmaster%2F.github%2Fsetup_onyxia.sh%C2%BB)

Ce dépôt contient tous les programmes du document de travail

**Abbas R., Carnot N., Lequien M., Quartier-la-Tente A. et Roux S. (2024)**, *En chemin vers la neutralité carbone. Mais quel chemin ?*,
Document de travail Insee, G2024/11.

Pour citer cet article :

    @article{inseeDTG202411,
      title={En chemin vers la neutralité carbone. Mais quel chemin ?},
      author={Abbas, Riyad and Carnot, Nicolas and Lequien, Matthieu and Quartier{-la-}Tente, Alain and Roux, S{\'e}bastien},
      journal={Document de travail Insee},
      number={G2024/11},
      year={2024},
      url={https://github.com/InseeFrLab/DT-way-to-net-zero}
    }
    
## Installation

Tous les programmes sont en R.
Pour installer les packages nécessaires dans les versions utilisées pour le document de travail, il faut tout d'abord installer le package `renv` et ensuite, une fois le projet chargé il suffit de
lancer le code `renv::restore()`.

Une [image docker](https://hub.docker.com/repository/docker/inseefrlab/dt-way-to-net-zero) a été construite pour assurer la reproductibilité complète de l’environnement utilisé pour ce document de travail. 
Elle peut également être directement utilisée avec [Onyxia](https://github.com/InseeFrLab/onyxia-web), la plateforme
*datascience* développée par l’[Insee](https://www.insee.fr/fr/accueil) en cliquant sur
[![Onyxia](https://img.shields.io/badge/Launch-Datalab-orange?logo=R)](https://datalab.sspcloud.fr/launcher/ide/rstudio?autoLaunch=false&service.image.custom.enabled=true&service.image.pullPolicy=%C2%ABAlways%C2%BB&service.image.custom.version=%C2%ABinseefrlab%2Fdt-way-to-net-zero%3Alatest%C2%BB&init.personalInit=%C2%ABhttps%3A%2F%2Fraw.githubusercontent.com%2Finseefrlab%2Fdt-way-to-net-zero%2Fmaster%2F.github%2Fsetup_onyxia.sh%C2%BB).

    
## Description des programmes et des fichiers

Les programmes suivants sont utilisés :

- `R/0_functions.R` : fonctions utilisées dans les autres scripts

- `R/1_estimations.R` : estimation des modèles et analyse de sensibilité (selon paramètres définis dans le fichier `parameters.csv`) qui sont ensuite sauvegardés dans le dossier `est_mod`.

- `R/2_variantes_BC.R` : estimation de 2 variantes au budget carbone : démarrage retardé (2028, 2033 ou 2038) et avec des cibles espacées dans le temps (2, 5 ou 10 ans).
Les modèles sont sauvegardés dans le dossier `est_mod`.

- `R/3_figures_article.R` : création de l'ensemble des graphiques du document de travail

- `R/4_statistiques_articles.R` : calcul des statistiques utilisées dans l'article, il est utilisé pour générer un fichier html lisible à partir du lien <https://inseefrlab.github.io/DT-way-to-net-zero/R/4_statistiques_article.html>.


Le dossier `DT` contient le programme pour générer le document de travail.

Les fichiers `.Rprofile`, `renv.lock` et le dossier `renv` sont associés à l'utilisation du package `renv` qui permet de contrôler les versions des packages utilisés.

Le dossier `.github`, non utile pour l'utilisateur, contient les fichiers relatifs à la création de l'image docker et l'initialisation du service Onyxia.