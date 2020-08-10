# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 11:27:05 2020

@author: Mark
"""

import tkinter
import tkinter.ttk as tk
from tkinter import TOP

import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backend_bases import MouseEvent, LocationEvent
from matplotlib.widgets import RectangleSelector

import numpy as np

class Figure:       
    """Parameters
        ----------
        frame : TK Frame
                The container for the figure.
        x : ARRAY or LIST of FLOATS
            The x-values that should be used for the plot.
        y : ARRAY or LIST of FLOATS
            The y-values used for the plot.
        params : DICTIONARY
            Should contain keys 'xlabel', 'ylabel', 'title'. Values should be
            STRING.
        dpi : INT, optional
            DESCRIPTION. The default is 100.

        Returns
        -------
        None.
        """
    def __init__(self, frame, x, y, params, dpi = 100):

        self.frame = frame
        self.x = x
        self.y = y
        
        self._getParams(params)
        
        self.fig, self.ax = plt.subplots(figsize=self.size, dpi=dpi)
        self.fig.patch.set_facecolor('0.9411')
        self.ax.plot(x, y)
        self._setLabels()

        self._setUp()
        
    def _getParams(self, params):

        if 'title' in params.keys():
            self.title = params['title']
        else:
            self.title = ''
        if 'xlabel' in params.keys(): 
            self.xlabel = params['xlabel']
        else:
            self.xlabel = ''
        if 'ylabel' in params.keys(): 
            self.ylabel = params['ylabel']
        else:
            self.ylabel = ''
        if 'colours' in params.keys():
            self.colours = params['colours']
        else:
            self.colours = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', 
                            '#9467bd','#8c564b', '#e377c2', '#7f7f7f', 
                            '#bcbd22', '#17becf']
        if 'size' in params.keys():
            self.size = params['size']
        else:
            self.size = (5,4)
        if 'axis_label_fontsize' in params.keys():
            self.axis_label_fontsize = params['axis_label_fontsize']
        else:
            self.axis_label_fontsize = 12
        
        if 'axis_tick_font' in params.keys():
            self.axis_tick_font = params['axis_tick_font']
        else: 
            self.axis_tick_font = 8
        
    def _setLabels(self):
        self.ax.set_xlabel(self.xlabel, fontsize=self.axis_label_fontsize)
        self.ax.set_ylabel(self.ylabel, fontsize=self.axis_label_fontsize)
        self.ax.set_title(self.title)
        
    def _setUp(self):
        self.ax.tick_params(direction='out', length=4, width=1, 
                                 colors='black', grid_color='black', 
                                 labelsize=self.axis_tick_font, grid_alpha=0.5)
        self.fig.tight_layout()
        self.chart = FigureCanvasTkAgg(self.fig, self.frame)
        self.chart.mpl_connect('button_press_event', self.doubleClkChart)
        rectprops = dict(facecolor='white', edgecolor = 'black',
                           alpha=0.2, fill=True)
        self.RS = RectangleSelector(self.ax, 
                                    self.selector,
                                       drawtype='box', useblit=True,
                                       button=[1, 3],  # don't use middle button
                                       minspanx=5, minspany=5,
                                       spancoords='pixels',
                                       interactive=True,
                                       rectprops = rectprops) 
        self.RS.set_active(True)
        
        self.chart.get_tk_widget().pack(side=TOP)
        
    def doubleClkChart(self, event): 
        if event.dblclick:
            self.zoomOut()
                     
    def rePlotFig(self, data, params):
        """
        This fuction replots the figure and re-sets the x and y limits'''

        Parameters
        ----------
        data : LIST of DICTIONARIES, where the keys in the DICTIONARIES are
            'x', and 'y', and the respective values are lists of floats.
            'idx' which is the lindex of the plotted line.
        params : DICTIONARY
            Can contain key 'normalize', with value 0 or 1. This determines
            if the lines should be plotted normalized.

        Returns
        -------
        None.

        """
        left,right = self.ax.get_xlim()
        bottom, top = self.ax.get_ylim()
        self.ax.clear()
        self._setLabels()
        
        def _getIdx(i):
            if 'idx' in i.keys():
                idx = i['idx']
            else:
                idx=0
            return idx
        
        def plot(x,y,idx):
            colour_idx = self._getColourIdx(idx)
            self.ax.plot(i['x'],i['y']/np.max(i['y']),
                                    c=self.colours[colour_idx])

        if 'normalize' in params.keys():
            if params['normalize'] == 1:
                for i in data:
                    x = i['x']
                    y = i['y']/np.max(i['y'])
                    plot(x,y,_getIdx(i))
            else:
                for i in data:
                    x = i['x']
                    y = i['y']
                    plot(x,y,_getIdx(i))         
        else:
            for i in data:
                x = i['x']
                y = i['y']
                plot(x,y,_getIdx(i))
        
        if 'rescale' in params.keys():
            if params['rescale'] == False:
                self.ax.set_xlim(left,right)
                self.ax.set_ylim(bottom,top)
        
        self.fig.tight_layout()
        self.chart.draw()
        
    def _getColourIdx(self, idx):
        colour_idx = idx
        if colour_idx >= len(self.colours):
            colour_idx = colour_idx%len(self.colours)
        return colour_idx
            
    def zoomOut(self):
        self.press = None
        self.ax.relim()
        self.ax.autoscale()
        self.chart.draw()
        
    def selector(self, eclick, erelease):
        
        if eclick.dblclick:
            '''In the case that the user double clicks, we dont want to use 
            the selector. We use zoomOut instead.
            '''
            pass
        else:
            xvals = [eclick.xdata,erelease.xdata]
            yvals = [eclick.ydata,erelease.ydata]
            self.ax.set_xlim(min(xvals),max(xvals))
            self.ax.set_ylim(min(yvals),max(yvals))
            self.chart.draw()