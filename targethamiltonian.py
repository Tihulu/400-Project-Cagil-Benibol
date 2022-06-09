#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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
for j in range(n):
    for k in range(n):
        F[j][k]=(1/(m.sqrt(n))) * np.exp(  ((2*np.pi*complex(0,1))/(4*n)) * (2*(j+1) - ((n+1)*(2*(k+1) - (n+1)))  )   )
    Qpos[j][j]=( m.sqrt((2*np.pi)/(4*n)) * (2*(j+1) - (n+1)))
Ppos = np.matmul(np.conj(F),Qpos,F)

#oscilatory basis
I=np.identity(n)
Posc = (complex(0,1)/(m.sqrt(2))) * Pz
Qosc = (1/(m.sqrt(2))) * Qz
x = np.kron(Qpos,I)
y = np.kron(I,Qosc)
Px = np.kron(Ppos,I)
Py = np.kron(I,Posc)

thetax= np.arctan(np.abs(y/x))
thetap= np.arctan(np.abs(Py/Px))
#P2x= np.matmul(Px,Px)
#print('P2x',P2x)
#V = lambda x : np.exp((-4*x)/c) * ( (Q2_4*(np.exp((-8*x)/c))) - k*(np.exp((-8*x)/c)) + lambd8*np.identity(16))

#V = lambda x: ((x**4) * 3) + 0.5 * (x**2)  
M1_4=29.167
M2=7.638
V = lambda x: M1_4*((1- np.exp(x/M2))**2)
#print('V ','\n',V(Qosc))

H = lambda x, Px : ((Px**2)/2) + V(x)

V=V(Qosc)
#Hmtrx = H(Qosc, Posc)
Hmtrx = ((Posc*Posc)/2) + V

#print('Q \n', Qpos)

#print('Px \n', Px)

#print('x \n', x)

#print('y \n', y)

print('H ','\n',Hmtrx)

print('V ','\n',V)

e_val,e_vec = eig(np.nan_to_num(Hmtrx))
print('eig values, \n', e_val)