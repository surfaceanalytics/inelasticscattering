# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 11:23:30 2021

@author: pielsticker
"""
import matplotlib.pyplot as plt
import random
import numpy as np
#%% Plotting function.
def plot_spectra(model, plot_sim=True):
    fig, ax = plt.subplots(1,1)

    ax.plot(model.unscattered_spectrum.x,
            model.unscattered_spectrum.lineshape)
    ax.plot(model.scattered_spectrum.x,
            model.scattered_spectrum.lineshape)
    legend = ["unscattered", "scattered"]

    if plot_sim:
        ax.plot(model.simulated_spectrum.x,
            model.simulated_spectrum.lineshape)
        legend += ["simulated"]


        text = ""
        for i, component in enumerate(
                model.scattering_medium.scatterer.loss_function.components):
            text += f"Component {i}: " + str(type(component).__name__) + "\n"
            for key, value in component.__dict__.items():
                if key in ["width", "position", "intensity", "B", "C",
                           "D", "Eg", "edge", "exponent", "fermi_width"]:
                    text += f"     {key}: {value}\n"

            ax.text(0.1, 0.9, text, horizontalalignment='left',
                    verticalalignment='top', transform=ax.transAxes,
                    fontsize=8)

    ax.legend(legend, fontsize=8)
    plt.tight_layout()
    plt.show()

#%% Optimization functions
def rmse(predictions, targets):
    """
    Calculate root mean square error.

    Parameters
    ----------
    predictions : array
        Array with predicted values.
    targets : array
        Array with target values.

    Returns
    -------
    float
        RMSE of the two arrays.

    """
    return np.sqrt(((predictions - targets) ** 2).sum())

def optimize_peak_height(model, grad_rate, rate):
    """


    Parameters
    ----------
    model : Model
        Inelastic scattering model. Needs to have a loaded scattered
        spectrum.
    grad_rate : float
        DESCRIPTION.
    rate : float
        DESCRIPTION.

    Returns
    -------
    rms3 : TYPE
        DESCRIPTION.

    """
    rand = random.choice([-1,1])
    inelastic_xsect = model.scattering_medium.scatterer.inelastic_xsect
    plus = inelastic_xsect + rand * rate

    model.scattering_medium.scatterer.inelastic_xsect = plus
    print(model.scattering_medium.scatterer.inelastic_xsect)
    model.scatterSpectrum()
    model.simulated_spectrum.normalize_minmax()
    diff1 = np.subtract(
        model.scattered_spectrum.lineshape,
        model.simulated_spectrum.lineshape)
    rms1 = np.sqrt(np.sum(np.square(diff1)))
    rms1 = rms1 / np.sum(model.simulated_spectrum.lineshape)
    print("rms1: " + str(rms1))
    rmse1 = rmse(model.simulated_spectrum.lineshape,
                 model.scattered_spectrum.lineshape)
    print(f"rmse1: {rmse1}")

    minus = inelastic_xsect - rand * rate
    model.scattering_medium.scatterer.inelastic_xsect = minus
    print(model.scattering_medium.scatterer.inelastic_xsect)
    model.scatterSpectrum()
    model.simulated_spectrum.normalize_minmax()
    diff2 = np.subtract(
        model.scattered_spectrum.lineshape,
        model.simulated_spectrum.lineshape)
    rms2 = np.sqrt(np.sum(np.square(diff2)))
    rms2 = rms2 / np.sum(model.simulated_spectrum.lineshape)
    print("rms2: " + str(rms2))

    gradient = (rms2-rms1) / (minus-plus)
    print("gradient: " + str(gradient))

    new_val = inelastic_xsect - gradient * grad_rate
    model.scattering_medium.scatterer.inelastic_xsect = new_val
    print("new_val: " + str(new_val))
    model.scatterSpectrum()
    model.simulated_spectrum.normalize_minmax()
    diff3 = np.subtract(
        model.scattered_spectrum.lineshape,
        model.simulated_spectrum.lineshape)
    rms3 = np.sqrt(np.sum(np.square(diff3)))
    rms3 = rms3 / np.sum(model.simulated_spectrum.lineshape)

    print("rms3: " + str(rms3))

    return rms3

def optimize_factor(model):
    grad_rate = 2
    rate = 0.1
    rand = random.choice([-1,1])
    norm_factor = model.scattering_medium.scatterer.norm_factor
    plus = norm_factor + rand * rate

    model.scattering_medium.scatterer.norm_factor = plus
    print(model.scattering_medium.scatterer.norm_factor)
    model.scatterSpectrum()
    diff1 = np.subtract(
        model.scattered_spectrum.lineshape,
        model.simulated_spectrum.lineshape)
    rms1 = np.sqrt(np.sum(np.square(diff1)))
    rms1 = rms1 / np.sum(model.simulated_spectrum.lineshape)
    print("rms1: " + str(rms1))

    minus = norm_factor - rand * rate
    model.scattering_medium.scatterer.norm_factor = minus
    print(model.scattering_medium.scatterer.norm_factor)
    model.scatterSpectrum()
    diff2 = np.subtract(
        model.scattered_spectrum.lineshape,
        model.simulated_spectrum.lineshape)
    rms2 = np.sqrt(np.sum(np.square(diff2)))
    rms2 = rms2 / np.sum(model.simulated_spectrum.lineshape)
    print("rms2: " + str(rms2))

    gradient = (rms2-rms1) / (minus-plus)
    print("gradient: " + str(gradient))

    new_val = norm_factor - gradient * grad_rate
    model.scattering_medium.scatterer.norm_factor = new_val
    print("new_val: " + str(new_val))
    model.scatterSpectrum()
    diff3 = np.subtract(
        model.scattered_spectrum.lineshape,
        model.simulated_spectrum.lineshape)
    rms3 = np.sqrt(np.sum(np.square(diff3)))
    rms3 = rms3 / np.sum(model.simulated_spectrum.lineshape)

    print("rms3: " + str(rms3))
    return rms3



def iterate_once(model, learning_rate):
    """


    Parameters
    ----------
    model : Model
        Inelastic scattering model. Needs to have a loaded scattered
        spectrum.
    learning_rate : TYPE
        DESCRIPTION.

    Returns
    -------
    rms3 : TYPE
        DESCRIPTION.

    """
    L = model.scattering_medium.scatterer.loss_function
    component_variables = [comp.__dict__ for comp in L.components]
    step_sizes = {
        'width':1,
        'position':0.5,
        'intensity':0.1,
        'B':100,
        'C':100,
        'D':10,
        'Eg':0,
        "edge":0.02,
        "exponent":0.001,
        "fermi_width": 0.01
        }
    limits = {
        'width':(0.01,20),
        'position':(0,200),
        'intensity':(0,10000),
        'B':(0,1E+9),
        'C':(0,1E+9),
        'D':(0.01,1000),
        'Eg':(0,0),
        "edge":(0.0,100.0),
        "exponent":(0.001,1),
        "fermi_width": (0.1,40)
        }
    gradients = {}
    config_1 = []
    config_2 = []

    randoms = []

    for idx, comp in enumerate(component_variables):
        c_1 = {}
        c_2 = {}
        for key, value in comp.items():
            if step_sizes[key] != 0:
                rand = random.choice([-1,1])
                randoms += [rand]
                val_1 = value + step_sizes[key] * rand
                val_2 = value - step_sizes[key] * rand
                c_1[key] = val_1
                c_2[key] = val_2
        config_1 += [c_1]
        config_2 += [c_2]

    for idx, config in enumerate(config_1):
        for key, value in config.items():
            L.editComponent(idx, {key: value})

    model.scattering_medium.scatterer.loss_function = L
    model.scatterSpectrum()
    model.simulated_spectrum.normalize_minmax()
    Diff1 = np.subtract(model.scattered_spectrum.lineshape,
                        model.simulated_spectrum.lineshape)
    RMS1 = np.sqrt(np.sum(np.square(Diff1)))
    RMS1 = RMS1 / np.sum(model.simulated_spectrum.lineshape)

    for idx, config in enumerate(config_2):
        for key, value in config.items():
            L.editComponent(idx, {key: value})

    model.scattering_medium.scatterer.loss_function = L
    model.scatterSpectrum()
    model.simulated_spectrum.normalize_minmax()
    Diff2 = np.subtract(model.scattered_spectrum.lineshape,
                        model.simulated_spectrum.lineshape)
    RMS2 = np.sqrt(np.sum(np.square(Diff2)))
    RMS2 = RMS2 / np.sum(model.simulated_spectrum.lineshape)

    for idx, comp in enumerate(component_variables):
        for key in comp.keys():
            Gradient = ((RMS2-RMS1) / (config_2[idx][key]-config_1[idx][key]))
            #gradients[key] = Gradient

            old_val = component_variables[idx][key]
            new_val = old_val - Gradient * learning_rate*step_sizes[key]

            if (new_val > limits[key][0]) and (new_val < limits[key][1]):
                L.editComponent(idx, {key: new_val})
            else:
                print(f"{key}={new_val} was outside its limits.")
                L.editComponent(idx, {key: old_val})

    model.scattering_medium.scatterer.loss_function = L
    model.scatterSpectrum()
    model.simulated_spectrum.normalize_minmax()

    print([[c.__dict__ for c in model.scattering_medium.scatterer.loss_function.components]])

    plot_spectra(model, plot_sim=True)

    diff3 = np.subtract(
        model.scattered_spectrum.lineshape,
        model.simulated_spectrum.lineshape)
    rms3 = np.sqrt(np.sum(np.square(diff3))) #/ np.sum(Sim.lineshape)

    return model, rms3#gradients, rms3

#%%

def optimize_spsa(
    model,
    limits=None,
    niter=100,
    a=1.0,
    alpha=0.602,
    c=1.0,
    gamma=0.101):
    """


    Parameters
    ----------
    model : Model
        Inelastic scattering model. Needs to have a loaded scattered
        spectrum.
    learning_rate : TYPE
        DESCRIPTION.

    Returns
    -------
    rms3 : TYPE
        DESCRIPTION.

    """
    A = 0.01 * niter

    gradients = {}
    rms_hist = {}

    for k in range(niter):
        print(k)
        ak = a/(k+1.0+A)**alpha
        ck = c/(k+1.0)**gamma

        L = model.scattering_medium.scatterer.loss_function
        component_variables = [comp.__dict__ for comp in L.components]

        config_1, config_2 = _get_configs(component_variables, ck)
# =============================================================================
#         for idx, comp in enumerate(component_variables):
#             if limits:
#                 if (new_val > limits[key][0]) and (new_val < limits[key][1]):
#                     L.editComponent(idx, {key: new_val})
#                 else:
#                     print(f"{key}={new_val} was outside its limits.")
#                     L.editComponent(idx, {key: old_val})
# =============================================================================


        rms1 = rms_for_config(model, config_1)
        rms2 = rms_for_config(model, config_2)

        for idx, comp in enumerate(component_variables):
            for key in comp.keys():
                grad = ((rms2-rms1) / (config_2[idx][key]-config_1[idx][key]))
                gradients[k] = grad

                old_val = component_variables[idx][key]

                # Parameter update
                new_val = old_val - grad * ak

                if limits:
                    if (new_val > limits[key][0]) and (new_val < limits[key][1]):
                        L.editComponent(idx, {key: new_val})
                    else:
                        print(f"{key}={new_val} was outside its limits.")
                        L.editComponent(idx, {key: old_val})

            model.scattering_medium.scatterer.loss_function = L
            model.scatterSpectrum()
            model.simulated_spectrum.normalize_minmax()

        print([[c.__dict__ for c in model.scattering_medium.scatterer.loss_function.components]])

        plot_spectra(model, plot_sim=True)

        diff3 = np.subtract(
            model.scattered_spectrum.lineshape,
            model.simulated_spectrum.lineshape)
        rms3 = np.sqrt(np.sum(np.square(diff3))) #/ np.sum(Sim.lineshape)

        rms_hist[k] = rms3

    return model, gradients, rms_hist

def _get_configs(component_variables, ck):
    config_1 = []
    config_2 = []

    for idx, comp in enumerate(component_variables):
        c_1 = {}
        c_2 = {}
        for key, value in comp.items():
            rand = random.choice([-1,1])
            val_1 = value + ck * rand
            val_2 = value - ck * rand

            c_1[key] = val_1
            c_2[key] = val_2
        config_1 += [c_1]
        config_2 += [c_2]

    return config_1, config_2

def rms_for_config(model, config):
    L = model.scattering_medium.scatterer.loss_function

    for idx, c in enumerate(config):
        for key, value in c.items():
            L.editComponent(idx, {key: value})

    model.scattering_medium.scatterer.loss_function = L
    model.scatterSpectrum()
    model.simulated_spectrum.normalize_minmax()

    diff = np.subtract(model.scattered_spectrum.lineshape,
                       model.simulated_spectrum.lineshape)
    rms = np.sqrt(np.sum(np.square(diff)))
    rms /= np.sum(model.simulated_spectrum.lineshape)

    return rms


