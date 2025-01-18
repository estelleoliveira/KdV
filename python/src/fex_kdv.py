#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 14 16:16:41 2023
Projet_simulation_numerique : fex, solution éxacte de KdV
@author: Estelle Oliveira, Galaad Barraud et Rémi Bourges
"""

# Import de la bibliothèque Numpy pour calculs numériques
import numpy as np

def fex(M,N,T,L,delta,U0,Uinf,x0,x,t):

    DELTA=delta/(np.sqrt((U0-Uinf)/12))
    c=Uinf+(U0-Uinf)/3
    T1,X=np.meshgrid(t,x)

    # ----------------------------------------------------------------
    #Calcul de fex, la solution éxacte de KdV
    # ----------------------------------------------------------------
    fex= U0/(np.cosh(((X - x0- c*T1)-np.round((X-x0-c*T1)/L)*L)/DELTA)**2)
    
    return fex