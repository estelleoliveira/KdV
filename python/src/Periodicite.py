# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 01:49:26 2023
Projet_simulation_numerique : Fonctions utiles
@author: Estelle Oliveira, Galaad Barraud et Rémi Bourges
"""

# Périodicité

def periodicite(j,M):
    # Calcul de jm
    jm = j-1
    if j == 0:
        jm = M-1
# Fin de calcul de jm

# Calcul de jp
    jp = j+1
    if j == M-1:
        jp = 0
# Fin de calcul de jp

# Calcul de jv
    jv = j-2
    if j == 1:
        jv = M-1
    elif j == 0:
        jv = M-2
# Fin de calcul de jv

# Calcul de jw
    jw = j+2
    if j == M-2:
        jw = 0
    elif j == M-1:
        jw = 1
# Fin de calcul de jw
    return jm, jp, jv, jw 

