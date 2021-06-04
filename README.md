# Implementation of different TSA for IOTA on OMNeT++

Implémentation et comparaison de différents TSA (Tips Selection Algorithm) pour la cryptomonnaie IOTA (https://www.iota.org/) : 

- IOTA : www.descryptions.com/Iota.pdf
- G-IOTA : https://ieeexplore.ieee.org/document/8845163
- E-IOTA : https://ieeexplore.ieee.org/document/9223294

Ce projet se base sur le git suivant : https://github.com/richardg93/TangleSim et requiert donc le logiciel OMNeT++ afin de simuler le Tangle (voir https://blog.iota.org/the-tangle-an-illustrated-introduction-4d5eae6fe8d4/). Plusieurs modifications ont été apportées à celui-ci (affichage avec graphviz, fichiers log, benchmark, implémentations des différents TSA).   

Les fichiers Tangle.cc, Tangle.h, TangleModules.cc, TangleSim.ned et omnetpp.ini sont les fichiers sources nécessaires pour pouvoir simuler les différents TSA avec OMNeT++. 

Le dossier data contient les fichiers logs, les statistiques de performances des différents TSA, l'affichage en .dot (graphviz) et des scripts Python : 

- Stats.py permet d'effectuer une moyenne des performances d'un TSA via les fichiers logs du dossier log. 
- Le notebook TSAResult.ipynb permet d'effectuer un barplot afin de comparer les performances des différents TSA.
- TangleGen.py permet d'afficher le Tangle avec un code couleur dans le dossier im. Il utilise les fichiers logs du dossier Traking. 
- RTracking contient les fichiers logs du git de référence.

Pour plus de détails voir le compte rendu. 
