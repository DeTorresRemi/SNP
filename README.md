# Semi and non parametric econometrics #
Projet basé sur l'article de CIZEK qui présente un estimateur AMSTLE qu'on généralisera au cas semi-paramétrique.

## Mode d'emploi ##

ATTENTION: Avant de commencer, quelques modifictions du code s'imposent:\
  1/ Dans AMSTLE_semiparametric.R,\
      ligne 9: mettre le nom du dossier personnel où est rangé AMSTLE_parametric.R\
  2/ Dans main.R,\
      ligne 3,4: mettre le nom du dossier où sont rangés AMSTLE_parametric.R et AMSTLE_semiparametric.R
      
Pour lancer les estimations par AMSTLE, il faut se placer dans main.R:\
  0/ Exécuter les lignes 1-21 pour charger les packages, les données, les reformater et fixer les paramètres.\
  1/ Exécuter la ligne 22 pour lancer l'estimation AMSTLE dans le cas paramétrique\
  2/ Exécuter la ligne 23 pour lancer l'estimation AMSTLE dans le cas semi-paramétrique
  
NB: les résultats ont été passés en copier commentaires à partir de la ligne 45 de main.R et également sauvegardés dans Résultats_sp.RData
