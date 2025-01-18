#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 14 16:16:41 2023
Projet_simulation_numerique : KDV modèle WANG_2008, Z-K 1965 et 6-P
@author: Estelle Oliveira, Galaad Barraud et Rémi Bourges
"""


# Import de la bibliothèque Numpy pour calculs numériques
import numpy as np
# Library for creation of static, animated, and interactive visualizations
import matplotlib.pyplot as plt
# Make live animations in Matplotlib by using Animation classes
from matplotlib import animation
import KdV_W as kdvw
import KdV_ZK as kdvzk
import KdV_6P as kdv6p
import fex_kdv as fexkdv


# *********************************
# Liste des paramètres numériques *
# *********************************
eta = 1
mu = 0.022
L = 2
T = 10
Uinf = 0
U0 = 1
M = 100
N = 3000
x0 = L/2
#paramètres de la fonction ZK 
L_ZK = 2
M_ZK = 100
T_ZK = 10
N_ZK = 3000
#paramètres de la fonction 6P
L_6P = 2
M_6P = 100
T_6P = 10
N_6P = 3000

dx = L / M
dt = T / N
x = np.linspace(0, L-dx, M)
t = np.linspace(0, T, N+1)
delta=0.022

#appel de fex
fexa=fexkdv.fex(M,N,T,L,delta,U0,Uinf,x0,x,t)
# appel du solveur W
fw=kdvw.kdv_solver(L, T, M, N, eta, mu, x0, Uinf, U0)
# appel du solveur ZK
fzk=kdvzk.ZK_schema(L_ZK, T_ZK, M_ZK, N_ZK, eta, mu, x0, Uinf, U0)
# appel du solveur sixP
f6p=kdv6p.sixP_schema(L_6P, T_6P, M_6P, N_6P, eta, mu, x0, Uinf, U0)

#%%
#création d'un plot de f_w et fex
fig1=plt.figure()  
ax1 = plt.axes()  
ax1.plot(fw[:,0], '-', label='f_W')
ax1.plot(fexa[:,0], '+', label='f_ex')
ax1.set_xlabel('x')
ax1.set_ylabel('f')
ax1.legend(loc='upper right')
plt.plot()

# ******************
# Fonction makefig *
# ******************
fig = plt.figure() 
ax = plt.axes()                                        # Créer des axes

def makefig(n):
    ax.cla()
    ax.plot(x, fw[:,n], '-', label='f_W')
    ax.plot(x, fzk[:,n], '-', label='f_ZK')
    ax.plot(x, f6p[:,n], '-', label='f_6P')
    ax.plot(x, fexa[:,n], 'm:', label='f_ex')
    ax.set_xlabel('x')
    ax.set_ylabel('f')
    ax.text(0.5, 2, 't ='+str(np.around(t[n], decimals=1))+' s')
    ax.legend(loc='upper right')

# **************************
# Création d'une animation *
# **************************


G = animation.FuncAnimation(fig, makefig, interval=40, frames=range(0, N, 10), repeat=True)