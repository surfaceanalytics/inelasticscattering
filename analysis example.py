# -*- coding: utf-8 -*-
"""
Created on Thu May  7 14:18:50 2020

@author: Mark
"""

from analysis import Analysis
#%%
A = Analysis()
A.fig1(figsize = (6,3), dpi = 100)
A.fig2(dpi = 100, n_events = 4)
A.fig3()
A.fig4(dpi=100,multiplot=0, inelastic_xsect=0.001)
A.fig5(dpi=100, multiplot=0)

A.exportExcel("test")

filename = r'C:\Users\Mark\ownCloud\Muelheim Group\Projects\Gas phase background\python code\gasscattering\data\He\ARM22\Ag3d vacuum EX340 ARM22.TXT'
A.fig6(dpi=100, multiplot=0, P=10, scatterer='He', filename=filename)


A.fig2(dpi = 100, n_events = 4, scatterer='N2', filename=filename)
A.exportExcel("Ag_N2")
