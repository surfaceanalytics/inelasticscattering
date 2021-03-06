# -*- coding: utf-8 -*-
"""
Created on Thu May  7 13:58:58 2020

@author: Mark
"""

import numpy as np
import json
import re
from model.model import *
from model.model import Model
from model.algorithms.algorithm4 import Algorithm4
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from itertools import cycle
import xlsxwriter

#%%
'''
This module is for creating figures of the simulations
format of data: 
data = [{'x':'','y':'','color':'', 'limit0':0, 'limit1':-1, 'x_label':'','y_label':'','title':''}]
'''
class Analysis:
    def __init__(self):
        self.filename1 = r'C:\Users\Mark\ownCloud\Muelheim Group\Projects\Gas phase background\python code\gasscattering\data\test_spectrum.TXT'
        self.filename2 = r'C:\Users\Mark\ownCloud\Muelheim Group\Projects\Gas phase background\python code\gasscattering\data\scatterers.json'

    def draw_fig(self, data, **kwargs):
        if 'ncol' in kwargs.keys():
            ncol = kwargs['ncol']
        else:
            ncol = 2
        if 'titlefont' in kwargs.keys():
            titlefont = kwargs['titlefont']
        else:
            titlefont = 8
        if 'axisfont' in kwargs.keys():
            axisfont = kwargs['axisfont']
        else:
            axisfont = 8
        if 'axislabelsize' in kwargs.keys():
            axislabelsize = kwargs['axislabelsize']
        else:
            axislabelsize = 6
        if 'figsize' in kwargs.keys():
            figsize = kwargs['figsize']
        else:
            figsize = (5,4)
        if 'dpi' in kwargs.keys():
            dpi = kwargs['dpi']
        else:
            dpi = 300
        if 'multiplot' in kwargs.keys():
            multiplot = kwargs['multiplot']
        else:
            multiplot = 1            
        
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
                           grid_color='black', labelsize=axislabelsize, 
                           grid_alpha=0.5)
            if len(legend) != 0:
                ax.legend(tuple(legend),fontsize = 8)
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
    
    def _checkKwargs(self, model, **kwargs):
        
        if 'filename' in kwargs.keys():
            model.loadSpectrum(kwargs['filename'])
        else:
            model.loadSpectrum(self.filename1)
            
        if 'scatterer' in kwargs.keys():
            s = kwargs['scatterer']
            model.setCurrentScatterer(s)
        else:
            model.setCurrentScatterer('default')
 
        if 'P' in kwargs.keys():
            P = kwargs['P']
            model.scattering_medium.setPressure(P)
        else:
            model.scattering_medium.setPressure(4)
            
        if 'inelastic_xsect' in kwargs.keys():
            x = kwargs['inelastic_xsect']
            model.scattering_medium.scatterer.inelastic_xsect = x
        else:
            model.scattering_medium.scatterer.inelastic_xsect = 0.003
            
        if 'norm_factor' in kwargs.keys():
            norm = kwargs['norm_factor']
            model.scattering_medium.scatterer.norm_factor = norm
        else:
            model.scattering_medium.scatterer.norm_factor = 1
            
        if 'n_events' in kwargs.keys():
            n = kwargs['n_events']
            model.n_events = n
        else:
            model.n_events = 50
            
        if 'n_iter' in kwargs.keys():
            n_iter = kwargs['n_iter']
            model.calculation.n_iter = n_iter
        else:
            model.calculation.n_iter = 10
            
        if 'inelastic_prob' in kwargs.keys():
            inelastic_prob = kwargs['inelastic_prob']
            model.scattering_medium.scatterer.inelastic_prob = inelastic_prob
        else:
            model.scattering_medium.scatterer.inelastic_prob = 0.4

        return model
        
    def _initModel(self, **kwargs):
        model = Model()
        model.loadScatterers(self.filename2)
        return model 
    
    def exportExcel(self, filename):
        filename = filename + '.xlsx'
        with xlsxwriter.Workbook(filename) as workbook:
            worksheet = workbook.add_worksheet()
            cols_per_spec = 2
            for i, d in enumerate(self.data):
                start_col = i * cols_per_spec
                worksheet.write(0,start_col, i)
                worksheet.write(1,start_col, d['x_label'])
                worksheet.write(1,start_col+1, d['y_label'])
                x = d['x']
                y = d['y'].tolist()
                for i, val in enumerate(x):
                    worksheet.write(2+i,start_col, val)
                for i, val in enumerate(y):
                    worksheet.write(2+i, start_col+1, val)
            workbook.close()
       
    def fig1(self, **kwargs):
        ''' This plots of the test spectrum and the loss function
        '''
        model = self._initModel()
        model = self._checkKwargs(model, **kwargs)
        model.algorithm_id = 0
        model.scatterSpectrum()
        # Loss function figure
        Ly = model.scattering_medium.scatterer.loss_function.lineshape
        Lx = model.scattering_medium.scatterer.loss_function.x
        
        Dx = np.array(model.loaded_spectra[0].x)
        Dy = (np.array(model.loaded_spectra[0].lineshape) 
              - np.min(model.loaded_spectra[0].lineshape) 
              / 100000)
        
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
        self.data = data_to_plot   
        self.draw_fig(self.data, **kwargs) 
        
    def fig2(self, **kwargs):
        '''This plots all the n-convolved line shapes, but does not re-scale
        them by their Poisson factors.
        '''
        model = self._initModel()
        model = self._checkKwargs(model, **kwargs)
        model.unscattered_spectrum = model.loaded_spectra[0]

        model.algorithm_id = 1
        model.scatterSpectrum()
        data = [{'x':model.unscattered_spectrum.x,
                  'y':(np.array(model.simulation.I[i,:]) 
                       - np.min(model.simulation.I[i,:])
                       / 100000),
                  'color':'',
                  'limit0':0,
                  'limit1':-60,
                  'x_label':'Kinetic energy [eV]',
                  'y_label':'Intensity [Arb.U.]',
                  'title':'',
                  'linestyle':'solid',
                  'markersize':1,
                  'linewidth':1.5,
                  'marker':''} for i in range(len(model.simulation.I[:,0]))]
        data += [{'x':model.unscattered_spectrum.x,
                  'y':np.array(model.simulation.P) / 100000,#np.max(model.P),
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
        self.data = data    
        self.draw_fig(self.data, multiplot=0, **kwargs)     

    def fig3(self, **kwargs):
        ''' This plots the loss functions of several gases'''
        model = self._initModel()
        model = self._checkKwargs(model, **kwargs)
        data = []
        for i in ['He','H2','N2','O2']:
            model.setCurrentScatterer(i)
            Ly = model.scattering_medium.scatterer.loss_function.lineshape
            Lx = model.scattering_medium.scatterer.loss_function.x
            data += [{'x':Lx,
                  'y':Ly,
                  'color':'',
                  'limit0':0,
                  'limit1':-11600,
                  'x_label':'Kinetic energy [eV]',
                  'y_label':'Intensity [Arb.U.]',
                  'title':'',
                  'linestyle':'solid',
                  'markersize':1,
                  'linewidth':1,
                  'marker':'',
                  'text':i}]
        self.data = data
        self.draw_fig(self.data, multiplot=1, figsize=(5,4), **kwargs)     
        
    def fig4(self, **kwargs):
        ''' This plots the Poisson distribution at several pressures.
        It only works for Algorithms 2 and 3
        '''
        model = self._initModel()
        model = self._checkKwargs(model, **kwargs)
        model.algorithm_id = 0
        model.algorithms[0] = Algorithm4
        data = []
        pressures = [4,10,20,40,80]
        for p in pressures:
            model.scattering_medium.setPressure(p)
            model.scatterSpectrum()
            x = list(range(len(model.simulation.poiss_inel)))
            y = model.simulation.poiss_inel
            data += [{'x':x,
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
        self.data = data    
        self.draw_fig(self.data, **kwargs)     
        
    def fig5(self, **kwargs):
        """ Plots the line shape with all contributions from gas scattering.
        The convolved spectra have been scaled by the Poisson factors and
        the angle factors.
        """
        model = self._initModel()
        model = self._checkKwargs(model, **kwargs)
        model.unscattered_spectrum = model.loaded_spectra[0]
        model.scattering_medium.setPressure(4)
        model.scattering_medium.scatterer.norm_factor = 1
        model.scattering_medium.scatterer.inelastic_xsect = 0.0038
        model.algorithm_id = 0
        model.algorithms[0] = Algorithm4
        model.scatterSpectrum()
        
        Iy = np.array(model.simulation.I) - np.min(model.simulation.I)
        Ix = model.loaded_spectra[-1].x
        inel_probs = model.simulation.poiss_inel
       
        data5 = []
        for i in range(len(Iy[1:10,0])):
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
              'y':(np.sum(model.simulation.inel, axis=1) 
                   - np.min(model.simulation.inel)),
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
              'y':(np.array(model.simulated_spectrum.lineshape) 
                   - np.min(model.simulated_spectrum.lineshape)),
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
        self.data = data5  
        self.draw_fig(self.data, **kwargs)     

    def fig6(self, **kwargs):
        ''' Plot the contributions to the total spectrum of Ag3d in He'''
        model = self._initModel()
        model = self._checkKwargs(model, **kwargs)
        model.algorithm_id = 0
        model.unscattered_spectrum = model.loaded_spectra[0]

        model.scatterSpectrum()
        
        data6 = []
        data6 += [{'x':model.simulated_spectrum.x,
              'y':model.simulated_spectrum.lineshape,
              'color':'',
              'limit0':0,
              'limit1':-1,
              'x_label':'Kinetic Energy [eV]',
              'y_label':'Intensity [Cts./Sec.]',
              'title':'Ag3d spectrum in 4 mbar He',
              'linestyle':'solid',
              'markersize':0,
              'linewidth':1.5,
              'marker':'',
              'text':'Sum'}]   
        data6 += [{'x':model.simulated_spectrum.x,
              'y':np.sum(model.simulation.inel, axis=1),
              'color':'',
              'limit0':0,
              'limit1':-1,
              'x_label':'Kinetic Energy [eV]',
              'y_label':'Intensity [Cts./Sec.]',
              'title':'Ag3d spectrum in 4 mbar He',
              'linestyle':'solid',
              'markersize':0,
              'linewidth':1.5,
              'marker':'',
              'text':'Inelastic contribution'}] 
           
        data6 += [{'x':model.simulated_spectrum.x,
              'y':model.simulation.P,
              'color':'',
              'limit0':0,
              'limit1':-1,
              'x_label':'Kinetic Energy [eV]',
              'y_label':'Intensity [Cts./Sec.]',
              'title':'Ag3d spectrum in 4 mbar He',
              'linestyle':'solid',
              'markersize':0,
              'linewidth':1.5,
              'marker':'',
              'text':'Input spectrum'}] 
        self.data = data6
        self.draw_fig(self.data, **kwargs)

    def fig7(self, **kwargs):
        '''
        This plots the angular scaling factors.
        It is only relevant for Algorithm2
        '''
        model = self._initModel()
        model.algorithm_id = 0
        filename1 = r'C:\Users\Mark\ownCloud\Muelheim Group\Projects\Gas phase background\python code\gasscattering\data\He\ARM22\Ag3d vacuum EX340 ARM22.TXT'
        model.loadSpectrum(filename1)
        filename2 = r'C:\Users\Mark\ownCloud\Muelheim Group\Projects\Gas phase background\python code\gasscattering\data\scatterers.json'
        model.loadScatterers(filename2)
        model.unscattered_spectrum = model.loaded_spectra[-1]
        model.setCurrentScatterer('He')
        model.scattering_medium.d_through_gas = 1000
        model.scattering_medium.density = 20
        model.scattering_medium.scatterer.norm_factor = 10
        model.n_events = 50
        
        model.scatterSpectrum()
        
        A = model.simulation.A
        
        A1 = np.sum(A, axis=1)
        A2 = np.sum(A, axis=0)
        
        plt.plot(A1[:10])
        plt.plot(A2[:10])
    
#%%
        
if __name__ == '__main__':
# Trend of non-scattered intensity with pressure
    Ip = []
    P = []
    model.algorithm_id = 3
    for p in [0,1,2,4,8,16,32,64,128]:
        model.scattering_medium.setPressure(p)
        
        P += [p]
        Ip += [model.simulation.inel_factor[0]]
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

