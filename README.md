# TIPE2024
Recherche de balise en course d'orientation

Ce repository est mon projet de TIPE de 2024. Le but de ce projet est d'implémenter puis de comparer des stratégies de recherche de balise.

# Conditions de recherches
Il y a N agents cherchant une balise ponctuelle se trouvant dans une zone modélisé par un polygone quelconque, avec N le premier apramètre du problème. Chaque agent peut librement se déplacer dans le plan. La vitesse de déplacement des agents est majorée par une constante qui est commune à chaque agent. Cette constante est un deuxième paramètre du problème. Pour qu'une balise soit détectée, elle doit se trouver à moins d'une distance d d'un agent. Cette distance d est un troisième paramètre du problème.

# Stratégies adoptées
Nous implémonterons 3 stratégies. La première est la stratégie aléatoire qui fera office d'expérience témoin. La seconde est une stratégie inspirée de la thermodynamique, dans laquelle les N agents sont assimilés à des molécules d'un gaz. La dernière est une stratégie utilisant une courbe de remplissage pour parcourir l'ensemble du polygone et être certain de passer à moins d'une distance d de la balise.

## Stratégie aléatoire
Dans cette stratégie, les N agents adoptent une direction aléatoire à chaque instant. Lorsqu'un déplacement fait sortir l'agent en dehors du polygone, alors on applique la règle du rebond, c'est-à-dire que l'agent rebondit sur le ou les côtés du polygone concerné en appliquant la loi de la réflection de Snell-Descartes.

## Stratégie Thermodynamique
Dans cette stratégie, à chaque instant, pour chaque agent, une nouvelle direction est calculée en additionnant la direction précédente avec un vecteur aléatoire et un vecteur calculé en fonction de la position des autres agents. Le vecteur aléatoire est multiplié par un coefficient a paramètre du problème. Et le vecteur issue de la position des autres agents est la somme des forces répulsives induites par les autres agents.
Une telle stratégie devrait permettre à l'ensemble des agents de recouvrir la surface du polygone et devrait ainsi aboutir à la découverte de la balise.

## Stratégie de la courbe de remplissage
Dans cette stratégie, les N agents peuvent sortir du polygone. En début de recherche on opère un découpage du polygone en un ensemble de triangles et pour chaque triangle on calcul la courbe de Polya correspondante. Puis on répartit les courbes de Polya à parcourir parmi les N agents. De cette manière, il existe au moins un agent qui passera suffisamment près de la balise pour la trouver.

# Résultats
Les résultats sont enregistrés dans un fichier local puis lu par le script lecteur.py qui produira un graphique représentant les temps de recherches par rapport aux paramètres du problème et indiquant aussi l'écart-type.
