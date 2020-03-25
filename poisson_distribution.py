# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 16:31:38 2020

@author: Mark
"""
import numpy as np
import matplotlib.pyplot as plt

n = 10
a = np.math.factorial(5)
sigma = 0.5
density = 0.5
distance = 80000

def Poisson(n, distance, sigma, density):
    p = (1/np.math.factorial(n)) * ((distance * density * sigma)**n) * np.exp(-1 * distance * density * sigma)
    return p


p = (1/np.math.factorial(n))
p1 = ((distance * density * sigma)**n)
p3 = distance * density * sigma

p2 = np.exp(-1 * distance * density * sigma)

probs = [Poisson(i, distance, sigma, density) for i in range(n)]

print(Poisson(0, distance, sigma, density))


#%%
sigma = 0.5
density = 0.5
step = 0.01
n = 10
l = np.arange(0,1000,step)
p = [Poisson(n,i,sigma,density) for i in l]
plt.plot(l,p)
area = np.sum(p)*step
