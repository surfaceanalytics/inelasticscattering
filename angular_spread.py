# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 14:52:52 2020

@author: Mark
"""
import numpy as np
import matplotlib.pyplot as plt

def AngDist(x, L=5000):
    t1 = 1/(4*np.pi)
    t2 = L
    t3 = 1 + L * ((np.sin(x / 180 * np.pi))**2)
    t4 = np.log(1+L)
    
    y = t1*t2/t3/t4
    
    return y

def multiConv(y,z,n):
    zz = []
    zz += [z]
    for i in range(n):
        z1 = zz[i]
        z2 = np.convolve(y,z1)
        #z3 = z2[int(len(z1)/2) - 1:len(z1)+int(len(z1)/2)-1]
        zz += [z2]
    return zz


x1 = np.arange(-90,90,0.1)
y1 = np.array([AngDist(i) for i in x1])
y1 = y1/np.sum(y1)
plt.plot(x1,y1)

x2 = np.arange(-90,90,0.1)
y2 = np.array([np.cos( i / 180*np.pi ) for i in x2])
y2 = y2/np.sum(y2)

z = multiConv(y1,y2,500)


area = []
angle = 30
delta_angle = (89-angle)/500

minidx = np.where(x1 < -1* angle)[-1][-1]
maxidx = np.where(x1 > angle)[-1][0]

for i in z:
    area += [np.sum(i[minidx:maxidx])]
    angle += delta_angle
    minidx = np.where(x1 < -1* angle)[-1][-1]
    maxidx = np.where(x1 > angle)[-1][0]


    
ratios = [area[j]/area[j-1] for j in range(len(area)) if j != 0]   
ratios=ratios[1:]    
    
plt.plot(area)
plt.plot(ratios)


for i in z:
    plt.plot(x1,i)
plt.show()



    
