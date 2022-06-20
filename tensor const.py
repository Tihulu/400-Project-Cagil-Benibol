#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 16:33:36 2022

@author: cagilbenibol
"""
import numpy as np
import math as m
import cmath as cm
from numpy.linalg import eig
n=16
Q2_4 = 29.166
k = 58.33
c = 45.84
lambd8 = 34.55
Qz = np.zeros((n,n))
Pz = np.zeros((n,n))
Qpos = np.zeros((n,n))
F = np.zeros((n,n))
Ppos = np.zeros((n,n))
for i,j in zip(range(n),range(1,n)):
        Pz[i][j] = -m.sqrt(j)
        Pz[j][i] = m.sqrt(j)
        Qz[i][j] = m.sqrt(j)
        Qz[j][i] = m.sqrt(j)
#print(Pz)

#Qpos Ppos
# for j in range(n):
#     for k in range(n):
#         F[j][k]=(1/(m.sqrt(n))) * np.exp(  ((2*np.pi*complex(0,1))/(4*n)) * (2*(j+1) - ((n+1)*(2*(k+1) - (n+1)))  )   )
#     Qpos[j][j]=( m.sqrt((2*np.pi)/(4*n)) * (2*(j+1) - (n+1)))
# Ppos = np.matmul(np.conj(F),Qpos,F)

#oscilatory basis
I=np.identity(n)
Posc = (complex(0,1)/(m.sqrt(2))) * Pz
Qosc = (1/(m.sqrt(2))) * Qz
x = np.kron(Qosc,I)
y = np.kron(I,Qosc)
Px = np.kron(Posc,I)
Py = np.kron(I,Posc)

thetax= np.arctan(np.abs(y/x))
thetap= np.arctan(np.abs(Py/Px))
'''
M1_4=29.167
M2=7.638
V = lambda x: M1_4*((1- np.exp(x/M2))**2)

V=V(Qosc)
Hmtrx = ((Posc*Posc)/2) + V
'''
V = lambda x : np.exp((-4*x)/c) * ( (Q2_4*(np.exp((-8*x)/c))) - k*(np.exp((-8*x)/c)) + lambd8*np.identity(16))
V=V(Qosc)

Hmtrx = ((Posc*Posc)/2) + V
#Hmtrx = ((Posc**2)/2) + ((Qosc*Qosc)/2) + ((0.275*(Qosc**4))/4)

Xp = [[0, 1], [1, 0]]
Yp = [[0, -complex(0,1)], [complex(0,1), 0]]
Zp = [[1, 0], [0, -1]]
Ip = [[1,0],[0,1]]
#C=np.zeros((4**4,4**4))
C=[]
Cconst=[]
#Cconst=np.zeros((4**4,4**4))
#matrix.diagonal().sum()

for i,a in zip((Xp,Yp,Zp,Ip),("X","Y","Z","I")):
    for j,b in zip((Xp,Yp,Zp,Ip),("X","Y","Z","I")):
        for k,c in zip((Xp,Yp,Zp,Ip),("X","Y","Z","I")):
            for l,d in zip((Xp,Yp,Zp,Ip),("X","Y","Z","I")):
                kl=np.kron(k,l)
                jkl=np.kron(j,kl)
                ijkl=np.kron(i,jkl)
                C=np.matmul(Hmtrx,ijkl)
                cweight=C.diagonal().sum()
                Cconst.append(cweight)
                print(a,b,c,d,cweight,)
                