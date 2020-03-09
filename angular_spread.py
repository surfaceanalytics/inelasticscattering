# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 14:52:52 2020

@author: Mark
"""
import numpy as np
import matplotlib.pyplot as plt

def AngDist(x, L=1000):
    t1 = 1/(4*np.pi)
    t2 = L
    t3 = 1 + L * ((np.sin(x / 180 * np.pi))**2)
    t4 = np.log(1+L)
    
    y = t1*t2/t3/t4
    
    return y

x1 = np.arange(-90,90,0.1)
y1 = np.array([AngDist(i) for i in x1])
y1 = y1/np.sum(y1)
plt.plot(x1,y1)


x2 = np.arange(-90,90,0.1)
y2 = np.array([np.cos( i / 180*np.pi ) for i in x2])
y2 = y2/np.sum(y2)

plt.plot(x2,y2)


z = np.convolve(y2,y1)
z = np.multiply(z,1/len(y2)) 
z = z[int(len(y1)/2) - 1:len(y1)+int(len(y1)/2)-1]
    

plt.plot(x1,z)
plt.plot(x1,y2)
plt.show()
