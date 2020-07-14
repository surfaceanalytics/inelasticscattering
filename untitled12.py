# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 10:19:46 2020

@author: Mark
"""

import numpy as np
from model import Model
import matplotlib.pyplot as plt

#%%

m = Model()
filename1 = r'C:\Users\Mark\ownCloud\Muelheim Group\Projects\Gas phase background\python code\gasscattering\data\test_spectrum.TXT'
m.loadSpectrum(filename1)
filename2 = r'C:\Users\Mark\ownCloud\Muelheim Group\Projects\Gas phase background\python code\gasscattering\data\scatterers.json'
m.loadScatterers(filename2)

m.setCurrentScatterer('He')

m.unscattered_spectrum = m.loaded_spectra[0]
m.scattering_medium.setPressure(4)
m.scattering_medium.scatterer.inelastic_xsect = 0.03
m.scattering_medium.scatterer.elastic_xsect = 0.04
m.scattering_medium.distance = 0.8
m.n_events = 50
m.algorithm_id = 2
m.scatterSpectrum()

pts_nm = 20
thickness = 100

sim = m.simulation
F = sim.F
A = sim.A
T = sim.T
F_sum = np.sum(F,axis=0) # elastic scattering
F_sum1 = np.sum(F,axis=1) # inelastic scattering
A_sum = np.sum(A,axis=0) # elastic scattering
A_sum1 = np.sum(A,axis=1) # inelastic scattering
T_sum = np.sum(T,axis=0) # elastic scattering
T_sum1 = np.sum(T,axis=1) # inelastic scattering
F_sums = np.zeros(len(F_sum1))
A_sums = np.zeros(len(A_sum1))
T_sums = np.zeros(len(T_sum1))
for i in range(5):
    i = i/pts_nm
    m.scattering_medium.distance = i
    m.scatterSpectrum()
    sim = m.simulation
    F_sums = np.add(F_sums,np.sum(sim.F, axis=1))
    A_sums = np.add(A_sums,np.sum(sim.A, axis=1))
    T_sums = np.add(T_sums,np.sum(sim.T, axis=1))
    plt.plot(np.sum(sim.T, axis=1))
plt.show()

a = np.sum(T, axis=1)
a_sum = np.sum(a)

plt.plot(T_sums)
plt.plot(A_sums)
plt.plot(F_sums)

test1 = T_sums
test2 = T_sums


plt.plot(test1)
plt.plot(test2)

#%%

def Poisson(n, distance, sigma, density):
    p = ((1/np.math.factorial(n)) 
        * ((distance * density * sigma)**n) 
        * np.exp(-1 * distance * density * sigma))
    return p


density = 58.6 # density of Ag in atoms/nm^3
sigma1 = 0.05
sigma2 = 0.08

# rows are number of times scattered (dimension 100)
# columns are distance (dimension 300)
sd1 = np.array([[Poisson(i,j/pts_nm,sigma1,density) for j in range(thickness*pts_nm)] for i in range(100)])
sd2 = np.array([[Poisson(i,j/pts_nm,sigma2,density) for j in range(thickness*pts_nm)] for i in range(100)])

y = np.empty([1,len(sd1[:,0])])
# iterate over distance
for i in range(len(sd1[0,:])):
    I = np.array([sd1[:,i]])
    J = np.array([sd2[:,i]])
    cross = np.dot(I.T,J)
    x = np.array([cross[0,:]])
    y = np.concatenate((y,x),axis=0)
y = y[1:,:]


#%%
# show Poisson distributions for n elastic scattering events over 
# distances of d/100
d = 100
for i in range(d):
    plt.plot(y[i,:10])
plt.show()

el_sums = np.sum(y,axis=0)
plt.plot([i/pts_nm for i in range(300)][:20], el_sums[:20])

axisfont = 10
titlefont = 10
axislabelsize = 6
fig, ax = plt.subplots(figsize = (5,4), dpi = 200)
fig.tight_layout()
title = 'Inelastic Poisson Profiles'
limit0 = 0
limit1 = 100
x_label = 'Nr. times scattered'
y_label = 'Probability'
ax.set_xlabel(x_label, fontsize = axisfont)
ax.set_ylabel(y_label, fontsize = axisfont)
ax.set_title(title, fontsize = titlefont)
ax.tick_params(direction='out', length=2, width=1, colors='black',
               grid_color='black', labelsize=axislabelsize, grid_alpha=0.5)
for i in list(range(len(sd1[0,:]))[1::5]):
    y = sd1[:,i]
    x = list(range(len(y)))
    ax.plot(x[limit0:limit1], 
        y[limit0:limit1])
plt.show()

#%%
# shows the probability of each n scattering event, as a function of thickness

axisfont = 10
titlefont = 10
axislabelsize = 6
fig, ax = plt.subplots(figsize = (5,4), dpi = 200)
fig.tight_layout()
title = 'Probability-Thickness Profiles for Inelastic Scattering'
limit0 = 0
limit1 = 100
x_label = 'Thickness [nm]'
y_label = 'Probability'
ax.set_xlabel(x_label, fontsize = axisfont)
ax.set_ylabel(y_label, fontsize = axisfont)
ax.set_title(title, fontsize = titlefont)
ax.tick_params(direction='out', length=2, width=1, colors='black',
               grid_color='black', labelsize=axislabelsize, grid_alpha=0.5)
sums_inel = []
for i in list(range(len(sd1[:,0]))[1::5]):
    y = sd1[i,:]
    sums_inel += [np.sum(y) * 1/pts_nm]
    x = list(range(len(y)))
    ax.plot(x[limit0:limit1], 
        y[limit0:limit1])
plt.show()

sd1_1 = np.sum(sd1,axis=1) * 1/pts_nm
plt.plot(sd1_1)

print(np.sum(sd1[0,:]))
plt.plot(sums_inel)

#%%

# calculate area for each of the n scattering events, integrated over a distance d

def Area(pts_nm, thickness, sigma, density, n):
    area = np.sum([Poisson(n,i/pts_nm, sigma,density)/pts_nm for i in range(thickness * pts_nm)])
    return area
pts_nm = 100
thickness = 100 # in nm
sigma = 0.08 # in nm^2
density = 58.6 # in atoms /nm^3
areas = [Area(pts_nm, thickness, sigma, density, n) for n in range(100)]

factor = 1/(sigma * density)

norm_areas = [i/factor for i in areas]
plt.plot(areas)
