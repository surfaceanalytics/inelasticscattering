# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 16:39:59 2020

@author: nicholls
"""

import matplotlib.pyplot as plt
from angular_spread_lorentzian import AngularSpreadCalc

#%% Example case of calculating angular spread for 500 eV:

# set up class, define no. of iterations, acceptance angle and energy
test_500eV_nist = AngularSpreadCalc(iterations=500,
                                acceptance_angle=60,
                                energy=500)
# load nist data
cross_section_nist = test_500eV_nist.load_nist_cross_section('DCS_He_500_.csv')
test_500eV_nist.plot_nist()

# run convolution
convolution_results_nist = test_500eV_nist.run_convolution()
# plot convolution results
test_500eV_nist.plot_convolution_results()
# limit by acceptance angle
test_500eV_nist.limit_by_acceptance_angle()
# plot limit by acceptance angle
test_500eV_nist.plot_angle_limited()
# calculate accepted area under curve
area_sum_nist = test_500eV_nist.calc_area_under_curve()
# plot area per iteration
test_500eV_nist.plot_area_under_curve()
# calculate area ratio between n and n-1 iterations
area_ratio_nist = test_500eV_nist.calc_area_ratio()
# plot area ratio between n and n-1 iterations
test_500eV_nist.plot_area_ratio()


#%% Example case of calculating angular spread for 500 eV with lorentzian:

# set up class, define no. of iterations, acceptance angle and energy
test_500eV_lorentz = AngularSpreadCalc(iterations=500,
                                acceptance_angle=60,
                                energy=500,
                                width=28)

# generate lorenztian
cross_section_lorentz = test_500eV_lorentz.gen_lorentzian_cross_section()
test_500eV_lorentz.plot_cross_section()

# run convolution
convolution_results_lorentz = test_500eV_lorentz.run_convolution()
# plot convolution results
test_500eV_lorentz.plot_convolution_results()
# limit by acceptance angle
test_500eV_lorentz.limit_by_acceptance_angle()
# plot limit by acceptance angle
test_500eV_lorentz.plot_angle_limited()
# calculate accepted area under curve
area_sum_lorentz = test_500eV_lorentz.calc_area_under_curve()
# plot area per iteration
test_500eV_lorentz.plot_area_under_curve()
# calculate area ratio between n and n-1 iterations
area_ratio_lorentz = test_500eV_lorentz.calc_area_ratio()
# plot area ratio between n and n-1 iterations
test_500eV_lorentz.plot_area_ratio()


