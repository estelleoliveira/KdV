#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 14 16:16:41 2023
Projet_simulation_numerique : KDV modèle WANG_2008
@author: Estelle Oliveira, Galaad Barraud et Rémi Bourges
"""

# Import de la bibliothèque Numpy pour calculs numériques
import numpy as np
import Periodicite as p
import fex_kdv as fexkdv
import time 

def kdv_solver(L, T, M, N, eta, mu, x0, Uinf, U0):
    # Définition de l'espace-temps
    dx = L / M  #pas de l'espace
    dt = T / N  #pas de temps
    x = np.linspace(0, L-dx, M)     #création de l'abscisse
    t = np.linspace(0, T, N+1)      #création du vecteur temporel
    # Initialisation de u (matrice résultat) 
    u = np.zeros((M, N+1))
    # Définition des constantes
    alpha = (eta * dt) / (3 * dx)
    S = ((mu**2) * dt) / (dx**3)
    delta=0.022
    DELTA=delta/(np.sqrt((U0-Uinf)/12))
   
   
# ----------------------------------------------------------------
#Calcul de fex, la solution éxacte de KdV
# ----------------------------------------------------------------
    fex= fexkdv.fex(M,N,T,L,delta,U0,Uinf,x0,x,t)

# ----------------------------------------------------------------
# Formule d'initialisation de toute la colonne 0
# ----------------------------------------------------------------
    u[:, 0] = U0/np.cosh((x - x0)/DELTA)**2
    
# ----------------------------------------------------------------
#Calcul du ratio r et rtheo
# ----------------------------------------------------------------
    r=dt/dx
    rtheo=2/(np.max(np.abs(fex))+(4*(mu**2))/(dx**2))
    print("rW  numérique ", r)
    print("rW  théorique max  ",rtheo)

    start = time.time()
# ----------------------------------------------------------------
#Calcul de u[jm,n+1]
# ----------------------------------------------------------------

# Boucle j pour le calcul de u[j,1]
    for j in range(0, M):
        
        jm, jp, jv, jw = p.periodicite(j,M)
        
        # Calcul de u[j,1], on pose A,B,C pour plus de visibilité
        A = u[j, 0]
        B = - ((alpha/2)*(u[jp, 0] + u[j, 0] + u[jm, 0])*(u[jp, 0] - u[jm, 0]))
        C = - ((S/2)*(u[jw, 0] - 2*u[jp, 0] + 2*u[jm, 0] - u[jv, 0]))
        u[j, 1] = A + B + C
        
# Boucle n et j pour le calcul de u[j-1,n+1]
    for n in range(1, N):
        for j in range(0, M):

            jm, jp, jv, jw = p.periodicite(j,M)
            
# Calcul de u[j-1,n+1], on pose D, E, F pour plus de visibilité
            D = -alpha*(u[jp, n] + u[j, n] + u[jm, n])*(u[jp, n] - u[jm, n])
            E = -S*(u[jw, n] - (2*u[jp, n]) + (2*u[jm, n]) - u[jv, n])
            F = u[jm, n] - u[jp, n] + u[jp, n-1]
            u[jm, n+1] = D + E + F
            
        # Tests habituel et messages
        umax=np.max(np.abs(u))
        if umax>1e6:    
            print('Instabilité numérique: boum au pas '+str(n))
            break
        if (n+1)%100==0:
            print('Pas '+str(n+1)+' sur '+str(N))
        if (n+1)==N:
            print('Dernier pas atteint')

    print('temps de calcul modèle Wang = ', time.ctime(time.time() - start)[13:19],' minute')
    return u