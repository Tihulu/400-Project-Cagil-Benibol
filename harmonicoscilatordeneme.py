#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 25 17:33:36 2022

@author: cagilbenibol
"""

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

x = np.kron(Qosc,I)
x2=np.matmul(x,x)
x4=np.matmul(x2,x2)
y = np.kron(I,Qosc)
y2=np.matmul(y,y)
y4=np.matmul(y2,y2)
Px = np.kron(Posc,I)
Py = np.kron(I,Posc)




#2D Harmonic Osc.
V=0.5*(x2+y2)
Hmtrx =  ((np.matmul(Px,Px))/2) + (  (np.matmul(Py,Py))/2   ) + V

# Harmonic osc in 1D
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
e_val2,e_vec2 = eig(np.nan_to_num(Hmtrx1))
print('eig values2, \n', e_val2)




#paulization

Xp = [[0, 1], [1, 0]]
Yp = [[0, -complex(0,1)], [complex(0,1), 0]]
Zp = [[1, 0], [0, -1]]
Ip = [[1,0],[0,1]]
C=[]
Cconst=[]


Xp = [[0, 1], [1, 0]]
Yp = [[0, -complex(0,1)], [complex(0,1), 0]]
Zp = [[1, 0], [0, -1]]
Ip = [[1,0],[0,1]]
C=[]
Cconst=[]

for i,a in zip((Xp,Yp,Zp,Ip),("X","Y","Z","I")):
    for j,b in zip((Xp,Yp,Zp,Ip),("X","Y","Z","I")):
        for k,c in zip((Xp,Yp,Zp,Ip),("X","Y","Z","I")):
            for l,d in zip((Xp,Yp,Zp,Ip),("X","Y","Z","I")):
                kl=np.kron(k,l)
                jkl=np.kron(j,kl)
                ijkl=np.kron(i,jkl)
                C=np.matmul(Hmtrx1,ijkl)
                cweight=C.diagonal().sum()
                Cconst.append(cweight)
                #(-446.84995595499294 * X ^ X ^ X ^ X) + \
                print("(",cweight," *",a," ^",b," ^",c," ^",d,")"," + \ ")
                