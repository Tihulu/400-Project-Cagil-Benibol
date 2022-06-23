#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import math as m
import cmath as cm
from numpy.linalg import eig
n=16
Qz = np.zeros((n,n))
Pz = np.zeros((n,n))
a = np.zeros((n,n))
adag = np.zeros((n,n))
for i,j in zip(range(n),range(1,n)):
        a[i][j] = m.sqrt(j)
        adag[j][i] = m.sqrt(j)

#oscilatory basis
I=np.identity(n)
Posc = (1/(m.sqrt(2))) * (-(complex(0,1)))  * (a - adag)
Qosc = (1/(m.sqrt(2))) * (a + adag)
x = np.kron(Qosc,I)
y = np.kron(I,Qosc)
Px = np.kron(Posc,I)
Py = np.kron(I,Posc)

att=(1/(m.sqrt(2))) * (Qosc + (complex(0,1) * Posc))


Q2=np.matmul(Qosc,Qosc)
Pcomplex=np.conjugate(Posc)
P2=np.matmul(Posc,Posc)

at=(np.matmul(a,adag) - Q2)*2

#V = lambda x: ((x**4) * 3) + 0.5 * (x**2)  
#V=0.5*(np.matmul(Qosc,Qosc))
Q2_4 = 29.166
k = 58.33
c = 45.84
lambd8 = 34.55

M1_4=29.167
M2=7.638

n=4
a_=1
lambdx = 0.005
lambdy = lambdx
lambdmix=0.001

x = np.kron(Qosc,I)

x2=np.matmul(x,x)
x4=np.matmul(x2,x2)
y = np.kron(I,Qosc)
y2=np.matmul(y,y)
y4=np.matmul(y2,y2)
Px = np.kron(Posc,I)
Py = np.kron(I,Posc)



#V = lambda x : np.exp((-4*x)/c) * ( (Q2_4*(np.exp((-8*x)/c))) - k*(np.exp((-8*x)/c)) + lambd8*np.identity(16))
#V = lambda x: M1_4*(np.matmul( (I - np.exp(x/M2)),(I - np.exp(x/M2))  ))
V = lambda x,y,lambdx,lambdy,lambdmix,a : 0.5*(x2 + y2) + (lambdx*(x4)) + (lambdy*(y4)) + (1/(a_**4))*lambdmix*(x4)*(y4)

#V=V(Qosc)
V=V(x,y,lambdx,lambdy,lambdmix,a)
Hmtrx =  ((np.matmul(Px,Px))/2) + (  (np.matmul(Py,Py))/2   ) + V
# 0.5*P2 +
Hmtrx1 = 0.5*P2 + 0.5*Q2
Hmtrx2 = np.matmul(a,adag) + 0.5*I
#print('Q \n', Qpos)

#print('Px \n', Px)

#print('x \n', x)

#print('y \n', y)

#print('H ','\n',Hmtrx)

#print('V ','\n',V)

e_val1,e_vec = eig(np.nan_to_num(Hmtrx))
print('eig values, \n', e_val1)
e_val2,e_vec2 = eig(np.nan_to_num(Hmtrx2))
print('eig values2, \n', e_val2)