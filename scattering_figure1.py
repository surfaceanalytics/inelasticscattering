# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 13:57:44 2020

@author: Mark
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 12:19:34 2020

@author: Mark
"""
import numpy as np
import json
import re
from model import *
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from itertools import cycle


#%%
# format of data: 
# data = [{'x':'','y':'','color':'', 'limit0':0, 'limit1':-1, 'x_label':'','y_label':'','title':''}]

def draw_fig(data, ncol=2, titlefont = 8, axisfont=8, axislabelsize = 6, figsize=(5,4), dpi=300, multiplot = 1):
    colors = mcolors.TABLEAU_COLORS
    keys = [key for key in colors.keys()]
    
    if len(data)==1:
        multiplot = 0
    else:
        multiplot = multiplot
    
    if multiplot == 0:
        fig, ax = plt.subplots(figsize = figsize, dpi = dpi)
        c = cycle(keys)
        legend = []
        for i, d in enumerate(data):
            linestyle = d['linestyle']
            marker = d['marker']
            markersize = d['markersize']
            linewidth = d['linewidth']
            fig.tight_layout()
            x = d['x']
            y = d['y']
            limit0 = d['limit0']
            limit1 = d['limit1']
            x_label = d['x_label']
            y_label = d['y_label']
            title = d['title']
            if d['color'] == '':
                color = next(c)
            elif type(d['color']) == int:
                color = colors[keys[d['color']]]
            elif type(d['color']) == str:
                color = d['color']
            ax.plot(x[limit0:limit1], 
                    y[limit0:limit1],
                    linestyle = linestyle, 
                    linewidth = linewidth,
                    marker = marker,
                    markersize = markersize, 
                    color = color)
            if 'text' in [k for k in d.keys()]:
                if len(d['text']) != 0:
                    legend += [d['text']] 
        ax.set_xlabel(x_label, fontsize = axisfont)
        ax.set_ylabel(y_label, fontsize = axisfont)
        ax.set_title(title, fontsize = titlefont)
        ax.tick_params(direction='out', length=2, width=1, colors='black',
                       grid_color='black', labelsize=axislabelsize, grid_alpha=0.5)
        
        if len(legend) != 0:
            ax.legend(tuple(legend))
    else:
        nrow = int(len(data)/2)
        ncol = ncol
        fig, axs = plt.subplots(nrow, ncol, figsize = figsize, dpi=dpi)
        fig.tight_layout()
        c = cycle(keys)
        for i, ax in enumerate(fig.axes):
            linestyle = data[i]['linestyle']
            markersize = data[i]['markersize']
            marker = data[i]['marker']
            linewidth = data[i]['linewidth']
            x = data[i]['x']
            y = data[i]['y']
            limit0 = data[i]['limit0']
            limit1 = data[i]['limit1']
            x_label = data[i]['x_label']
            y_label = data[i]['y_label']
            title = data[i]['title']
            if type(data[i]['color']) == int:
                color = colors[keys[data[i]['color']]]
            elif type(data[i]['color']) == str:
                if data[i]['color'] == '':
                    color = next(c)
                else:
                    color = data[i]['color']
            #color = colors[keys[data[i]['color']]]
            ax.plot(x[limit0:limit1], 
                    y[limit0:limit1], 
                    linestyle = linestyle,
                    linewidth = linewidth,
                    markersize = markersize,
                    marker = marker,
                    color = color)
            ax.set_xlabel(x_label, fontsize = axisfont)
            ax.set_ylabel(y_label, fontsize = axisfont)
            ax.set_title(title, fontsize = titlefont)
            ax.tick_params(direction='out', length=2, width=1, colors='black',
                           grid_color='black', labelsize=axislabelsize, grid_alpha=0.5)
            if 'text' in [key for key in data[i].keys()]:
                label = data[i]['text']
                text_y = 0.85 * (max(ax.get_ylim()))
                text_x = 0.05 * (max(ax.get_xlim()))
                spl = re.findall(r'(\w+?)(\d+)', label)
                if len(spl)==0:
                    pass
                else:
                    spl = re.findall(r'(\w+?)(\d+)', label)[0]
                    label = spl[0] + '$' + '_' + spl[1] + '$'
                ax.text(text_x, text_y, label, fontsize=12)
    plt.show()
    
#%% 
# Load data

model = Model()
filename1 = r'C:\Users\Mark\ownCloud\Muelheim Group\Projects\Gas phase background\python code\gasscattering\data\test_spectrum.TXT'
model.loadSpectrum(filename1)
filename2 = r'C:\Users\Mark\ownCloud\Muelheim Group\Projects\Gas phase background\python code\gasscattering\data\scatterers.json'
model.loadScatterers(filename2)

model.setCurrentScatterer('He')

model.unscattered_spectrum = model.loaded_spectra[0]
model.scattering_medium.setPressure(4)
#model.scattering_medium.d_through_gas = 50000
#model.scattering_medium.density = 0.005
model.scattering_medium.scatterer.inelastic_xsect = 0.003
model.scattering_medium.scatterer.elastic_xsect = 0.04
model.n_events = 50
model.algorithm2()
    
T = model.T

    
#%%
# Loss function figure
    
Ly = model.scattering_medium.scatterer.loss_function.lineshape
Lx = model.scattering_medium.scatterer.loss_function.x

Dx = model.loaded_spectra[0].x
Dy = model.loaded_spectra[0].lineshape / 100000

data_to_plot = []
data_to_plot += [{'x':Lx,
                  'y':Ly,
                  'color':0,
                  'limit0':0,
                  'limit1':800,
                  'x_label':'Energy loss [eV]',
                  'y_label':'Cross section [Arb.U.]',
                  'title':'Test Loss Function',
                  'linestyle':'solid',
                  'markersize':1,
                  'linewidth':1,
                  'marker':''}]
data_to_plot += [{'x':Dx,
                  'y':Dy,
                  'color':1,
                  'limit0':800,
                  'limit1':-1,
                  'x_label':'Kinetic energy [eV]',
                  'y_label':'Intensity [Arb.U.]',
                  'title':'Test Electron Energy Distribution',
                  'linestyle':'solid',
                  'markersize':1,
                  'linewidth':1,
                  'marker':''}]
   
draw_fig(data_to_plot, figsize=(5,2)) 

#%%
# Plot all the inelastically scattered spectra

model.setCurrentScatterer('He')

n_events = 5
data1 = [{'x':np.arange(model.start,model.stop,model.step),
          'y':model.I[i,:] / 100000,#np.max(model.I[i,:]),
          'color':'',
          'limit0':0,
          'limit1':-60,
          'x_label':'Kinetic energy [eV]',
          'y_label':'Intensity [Arb.U.]',
          'title':'',
          'linestyle':'solid',
          'markersize':1,
          'linewidth':1.5,
          'marker':''} for i in range(n_events)]
data1 += [{'x':np.arange(model.start,model.stop,model.step),
          'y':model.P / 100000,#np.max(model.P),
          'color':'',
          'limit0':0,
          'limit1':-60,
          'x_label':'Kinetic energy [eV]',
          'y_label':'Intensity [Arb.U.]',
          'title':'',
          'linestyle':'solid',
          'markersize':1,
          'linewidth':1.5,
          'marker':''}]
    

data1 += [{'x':np.arange(model.start,model.stop,model.step)[::6],
          'y':model.bulk_spectrum[::6] / np.max(model.bulk_spectrum),
          'color':'grey',
          'limit0':0,
          'limit1':-10,
          'x_label':'Kinetic energy [eV]',
          'y_label':'Intensity [Arb.U.]',
          'title':'',
          'linestyle':'solid',
          'markersize':1.5,
          'linewidth':1,
          'marker':''}]


draw_fig(data1, multiplot=0, figsize=(3,2))     


#%%
# Loss functions of all gases
data2 = []
for i in ['He','H2','N2','O2']:
    model.setCurrentScatterer(i)
    Ly = model.scattering_medium.scatterer.loss_function.lineshape
    Lx = model.scattering_medium.scatterer.loss_function.x
    data2 += [{'x':Lx,
          'y':Ly,
          'color':'',
          'limit0':0,
          'limit1':-1600,
          'x_label':'Kinetic energy [eV]',
          'y_label':'Intensity [Arb.U.]',
          'title':'',
          'linestyle':'solid',
          'markersize':1,
          'linewidth':1,
          'marker':'',
          'text':i}]
    
draw_fig(data2, multiplot=1, figsize=(5,4))     

    
#%%
# Plot Poisson distribution
model.setCurrentScatterer('He')
model.unscattered_spectrum = model.loaded_spectra[0]
model.scattering_medium.scatterer.inelastic_xsect = 0.003
model.scattering_medium.scatterer.elastic_xsect = 0.04
model.n_events = 50


data3 = []
for p in [4,10,20,40,80]:
    model.scattering_medium.setPressure(p)
    model.algorithm2()
    x = list(range(len(model.inel_probs)))
    y = model.inel_probs
    data3 += [{'x':x,
      'y':y,
      'color':'',
      'limit0':0,
      'limit1':-40,
      'x_label':'Number of times scattered (n)',
      'y_label':'Probability',
      'title':'',
      'linestyle':'dashed',
      'markersize':4,
      'linewidth':1.5,
      'marker':'o',
      'text':'p: ' + str(p) + ' mbar'}]
    
draw_fig(data3, multiplot=0, figsize=(4,3))     

    

#%%
# Trend of non-scattered intensity with pressure
Ip = []
P = []
for p in [0,1,2,4,8,16,32,64,128]:
    model.scattering_medium.setPressure(p)
    model.algorithm2()
    P += [p]
    Ip += [model.inel_probs[0]]
plt.plot(P, Ip)

data4 = []
data4 += [{'x':P,
  'y':Ip,
  'color':'',
  'limit0':0,
  'limit1':-1,
  'x_label':'Pressure [mbar]',
  'y_label':'Relative Intensity [Arb.U.]',
  'title':'Intensity vs. Pressure for scattering through He',
  'linestyle':'solid',
  'markersize':3,
  'linewidth':1.5,
  'marker':'o',
  'text':''}]
    
draw_fig(data4, multiplot=0, figsize=(4,3))     

#%%

model.scattering_medium.setPressure(40)
model.scattering_medium.scatterer.inel_angle_factor = 1
model.scattering_medium.scatterer.inelastic_xsect = 0.003
model.scattering_medium.scatterer.el_angle_factor = 1
model.scattering_medium.scatterer.elastic_xsect = 0

model.algorithm2()

Iy = model.I
Ix = model.loaded_spectra[-1].x
inel_probs = model.inel_probs
non_scattered = model.non
total_inel = model.inel

data5 = []
for i in range(len(Iy[:10,0])):
    data5 += [{'x':Ix,
      'y':Iy[i,:] * inel_probs[i],
      'color':'',
      'limit0':0,
      'limit1':-1,
      'x_label':'Kinetic Energy [eV]',
      'y_label':'Intensity [Cts./Sec.]',
      'title':'',
      'linestyle':'solid',
      'markersize':0,
      'linewidth':1.5,
      'marker':'o',
      'text':''}] 

data5 += [{'x':Ix,
      'y':total_inel,
      'color':'grey',
      'limit0':0,
      'limit1':-1,
      'x_label':'Kinetic Energy [eV]',
      'y_label':'Intensity [Cts./Sec.]',
      'title':'Scattered portion of test spectrum in 4 mbar He',
      'linestyle':'dashed',
      'markersize':0,
      'linewidth':1.5,
      'marker':'o',
      'text':''}]  
    
data5 += [{'x':Ix,
      'y':model.simulated_spectrum.lineshape,
      'color':'',
      'limit0':0,
      'limit1':-1,
      'x_label':'Kinetic Energy [eV]',
      'y_label':'Intensity [Cts./Sec.]',
      'title':'Test spectrum in 4 mbar He',
      'linestyle':'solid',
      'markersize':0,
      'linewidth':1.5,
      'marker':'o',
      'text':''}]   
  
draw_fig(data5, multiplot=0, figsize=(4,3))     

#sumIy = np.sum(Iy,axis=0)

#plt.plot(inel_probs)
#%%








#%%
c = [key for key in data2[0].keys()]
print('text' in  c)
a = {'a':1,'b':2}
b = [key for key in a.keys()]
print('a' in b)

print('unscatterd probability: ' + str(model.p_unscattered))
print('inelastic probability: ' + str(np.sum(model.p_inel)))
print('elastic probability: ' + str(model.p_el))

print('total primary intensity: '+str(model.p_unscattered + model.p_el))

print('sum probability: '+str(np.sum(T)))
t = np.sum(T)
t1 = np.sum(T, axis = 1)
t11 = model.inel_probs
t2 = np.sum(T,axis=0).T
t22 = model.el_probs
print(np.sum(t2[1:]))

for i in [t1,t2]:
    plt.plot(i)
plt.show()
sum_t1 = np.sum(t1)
plt.pcolor(np.array(T[:limit,:limit]))



#%%
'''    
trends = []
for i in range(50):
    f = i * 10000
    model.scattering_medium.d_through_gas = f
    model.algorithm2()
    trends += [model.inel_probs]
    plt.plot(trends[-1])
plt.show()

totals = []
for i in range(len(trends)):
    x = 0
    for j in range(len(trends[i])):
        x += trends[i][j]
    totals += [x]
        

plt.plot(totals)

s1 = np.sum(model.el_probs[1:])
s2 = np.sum(model.inel_probs[1:])
'''

