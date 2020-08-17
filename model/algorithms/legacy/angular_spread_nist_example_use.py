# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 16:39:59 2020

@author: nicholls
"""

import matplotlib.pyplot as plt
from angular_spread_nist import AngularSpreadCalc

#%% Example case of calculating angular spread for 1000 eV:

# set up class, define no. of iterations, acceptance angle and energy
test_1000eV = AngularSpreadCalc(iterations=500,
                                acceptance_angle=60,
                                energy=1000)
# load data
cross_section = test_1000eV.load_nist_cross_section('DCS_He_1000_.csv')
# run convolution
convolution_results = test_1000eV.run_convolution()
# plot convolution results
test_1000eV.plot_convolution_results()
# limit by acceptance angle
test_1000eV.limit_by_acceptance_angle()
# plot limit by acceptance angle
test_1000eV.plot_angle_limited()
# calculate accepted area under curve
area_sum = test_1000eV.calc_area_under_curve()
# plot area per iteration
test_1000eV.plot_area_under_curve()
# calculate area ratio between n and n-1 iterations
area_ratio = test_1000eV.calc_area_ratio()
# plot area ratio between n and n-1 iterations
test_1000eV.plot_area_ratio()


#%% To run for all downloaded nist data and save in a dictionary:

dict_of_data = [{'name': '50eV', 'energy': 50, 'filename': 'DCS_He_50_eV.csv'},
                {'name': '100eV', 'energy': 100, 'filename': 'DCS_He_100_.csv'},
                {'name': '200eV', 'energy': 200, 'filename': 'DCS_He_200_.csv'},
                {'name': '250eV', 'energy': 250, 'filename': 'DCS_He_250_.csv'},
                {'name': '300eV', 'energy': 300, 'filename': 'DCS_He_300_.csv'},
                {'name': '400eV', 'energy': 400, 'filename': 'DCS_He_400_.csv'},
                {'name': '500eV', 'energy': 500, 'filename': 'DCS_He_500_.csv'},
                {'name': '600eV', 'energy': 600, 'filename': 'DCS_He_600_.csv'},
                {'name': '700eV', 'energy': 700, 'filename': 'DCS_He_700_.csv'},
                {'name': '750eV', 'energy': 750, 'filename': 'DCS_He_750_.csv'},
                {'name': '800eV', 'energy': 800, 'filename': 'DCS_He_800_.csv'},
                {'name': '900eV', 'energy': 900, 'filename': 'DCS_He_900_.csv'},
                {'name': '1000eV', 'energy': 1000, 'filename': 'DCS_He_1000_.csv'},
                {'name': '1500eV', 'energy': 1500, 'filename': 'DCS_He_1500_.csv'},
                {'name': '2000eV', 'energy': 2000, 'filename': 'DCS_He_2000_.csv'}]

# iterate through dictionary
for entry in dict_of_data:
    instance = AngularSpreadCalc(iterations=500,
                                 acceptance_angle=60,
                                 energy=entry['energy'])

    entry['cross_section'] = instance.load_nist_cross_section(entry['filename'])
    # run convolution
    entry['convolution_results'] = instance.run_convolution()
    # plot convolution results
    instance.plot_convolution_results()
    # limit by acceptance angle
    instance.limit_by_acceptance_angle()
    # plot limit by acceptance angle
    instance.plot_angle_limited()
    # calculate accepted area under curve
    entry['area_sum'] = instance.calc_area_under_curve()
    # plot area per iteration
    instance.plot_area_under_curve()
    # calculate area ratio between n and n-1 iterations
    entry['area_ratio'] = instance.calc_area_ratio()
    # plot area ratio between n and n-1 iterations
    instance.plot_area_ratio()

    entry['class'] = instance


#%% Writing functions to plot data

def plot_dif_cross_sections(data_dict, energy_list):
    """ plot area under curve, at particular energies, after convolution.
    Inputs:
        data_dict:  list of dictionaries containing analysis results
        energy_list: list of energies to be plotted
        xlim:   maximum value to plot on x axis
    """

    for i in data_dict:
        if i['energy'] in energy_list:
            plt.plot(i['class'].emitted_elctn_x, i['cross_section'],
                     label=i['name'])
    plt.legend()
    plt.xlabel('theta (degrees)')
    plt.ylabel('dsigma / d omega (a0^2/sr)')
    plt.title('Differential Elastic-Scattering Cross Section')
    plt.xticks([-90, -60, -30, 0, 30, 60, 90])
    #plt.savefig('diff elastic cross section.png', dpi=600)
    plt.show()


def plot_area_under_curves(data_dict, energy_list, xlim=100):
    """ plot area under curve, at particular energies, after convolution.
    Inputs:
        data_dict:  list of dictionaries containing analysis results
        energy_list: list of energies to be plotted
        xlim:   maximum value to plot on x axis
    """

    for i in data_dict:
        if i['energy'] in energy_list:
            plt.plot(i['area_sum'], label=str(i['name']))

    plt.title('Area under curve')
    plt.xlabel('No. of iterations')
    plt.xlim(0, xlim)
    plt.ylabel('Intensity a.u.')
    plt.legend()
    #plt.savefig('area under curve.png', dpi=100)
    plt.show()

def plot_area_ratio_change(data_dict, energy_list, xlim=100):
    """ Plot ratio of area under curve for iteration n vs n-1, at given energies.
    Inputs:
        data_dict:  list of dictionaries containing analysis results
        energy_list: list of energies to be plotted
        xlim:   maximum value to plot on x axis
    """

    for i in data_dict:
        if i['energy'] in energy_list:
            plt.plot(range(1, len(i['area_ratio'])+1), i['area_ratio'],
                     label=str(i['name']))
    plt.title('Intensity ratio change per iteration')
    plt.xlabel('Iterations')
    plt.xlim(0, xlim)
    plt.ylabel('Area Ratio between iterations')
    plt.legend()
    #plt.savefig('Ratio change per iteration compared.png', dpi=600)
    plt.show()

#%% Use of functions to create plots in manuscript
#selected energies to use:
energies = [100, 300, 500, 1000]

# Figure a
plot_dif_cross_sections(dict_of_data, energies)

# Figure c
plot_area_under_curves(dict_of_data, energies, 100)

# Figure d
plot_area_ratio_change(dict_of_data, energies, 100)


#%% To plot data for all energies:
# to get list of all available energies
energies = [i['energy'] for i in dict_of_data]

plot_dif_cross_sections(dict_of_data, energies)

plot_area_under_curves(dict_of_data, energies, 100)

plot_area_ratio_change(dict_of_data, energies, 100)

#%% To plot how the energy loss in the first iteration is dependent upon energy

for i in dict_of_data:
    plt.scatter(i['energy'], i['area_ratio'][0], c='black', s=20)
plt.title('Energy dependence of intensity ratio change in first iteration')
plt.xlabel('Energy / eV')
plt.ylabel('Area ratio change after 1st iteration')
plt.show()
