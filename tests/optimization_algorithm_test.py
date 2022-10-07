# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 15:52:00 2020

@author: Mark
"""
import matplotlib.pyplot as plt
import os

from model.model import Model
from converters.data_converter import DataConverter
from model.base_model import MeasuredSpectrum, Gauss


from model.shirley_optimization import (plot_spectra, optimize_factor,
                                    optimize_peak_height, iterate_once,
                                    optimize_spsa)

#%% First set-up the model.
abspath = os.path.abspath(__file__)
directory = os.path.dirname(abspath)
directory = r"C:\Users\pielsticker\Lukas\MPI-CEC\Projects\Gas phase background\inelasticscattering"# os.getcwd()

spectrum_filename = "Shirley/SiO2 O1s"#"He/Au, Ag, Cu in He"

spectrum_path = os.path.join(
    *[directory, "data", spectrum_filename + ".vms"] #"",
)

converter = DataConverter()
converter.load(spectrum_path)

model = Model()

for i, data_dict in enumerate(converter.data):
    data = converter.data[i]
    spectrum =  MeasuredSpectrum(data["data"]["x"], data["data"]["y0"])
    spectrum.normalize_minmax()
    model.loaded_spectra += [spectrum]

model.scattered_spectrum = model.loaded_spectra[0]
model.unscattered_spectrum = model.loaded_spectra[1]


plot_spectra(model, plot_sim=False)


#%% Set up scatterer and loss function
scatterer_path = os.path.join(
    *[directory, "data", "scatterers2.json"]
)

model.loadScatterers(scatterer_path)
model.setCurrentScatterer('Shirley')

model.algorithm_id = 0
model.calculation.n_iter = 50
model.algorithm_option = "bulk"
model.changeAlgorithm(1)

model.scattering_medium.scatterer.inelastic_xsect = 0.003
model.scattering_medium.scatterer.norm_factor = 0.88

# Only take the first 200 eV of the loss function
model.scattering_medium.scatterer.loss_function.stop = 200

# Initial params
# =============================================================================
# init_params = {
#     "Shirley": {
#         "edge": 9.5,
#         "exponent": 0.1,
#         "fermi_width": 20,
#         "intensity": 1
#         },
#     "Additional Gaussian 1": {
#         "peak_type": Gauss,
#         "position": 25,
#         "width": 10,
#         "intensity": 7,
#     },
#     "Additional Gaussian 2": {
#         "peak_type": Gauss,
#         "position": 60,
#         "width": 10,
#         "intensity": 0,
#     },
# }
# =============================================================================
init_params = {
    "Shirley": {
        "edge": 20,
        "exponent": 0.1,
        "fermi_width": 20,#0.5,
        "intensity": 500
        },
    "Additional Gaussian 1": {
        "peak_type": Gauss,
        "position": 15,
        "width": 5,
        "intensity": 150,
    },
    "Additional Gaussian 2": {
        "peak_type": Gauss,
        "position": 60,
        "width": 10,
        "intensity": 15,
    },
}

model.scattering_medium.scatterer.loss_function.components[0].edge = init_params["Shirley"]["edge"]
model.scattering_medium.scatterer.loss_function.components[0].exponent = init_params["Shirley"]["exponent"]
model.scattering_medium.scatterer.loss_function.components[0].fermi_width = init_params["Shirley"]["fermi_width"]
model.scattering_medium.scatterer.loss_function.components[0].intensity = init_params["Shirley"]["intensity"]

for key, d in init_params.items():
    if key != "Shirley":
        if d["intensity"] != 0:
            gauss = d["peak_type"](
                position=d["position"],
                width=d["width"],
                intensity=d["intensity"]
                )
            model.scattering_medium.scatterer.loss_function.addComponent(gauss)

model.scattering_medium.scatterer.loss_function.buildLine()

model.scatterSpectrum()
model.simulated_spectrum.normalize_minmax()

plt.plot(model.scattering_medium.scatterer.loss_function.x,
         model.scattering_medium.scatterer.loss_function.lineshape)
plt.legend(["Loss function"])
plt.show()

plot_spectra(model, plot_sim=True)


#%%
# =============================================================================
# grad_rate = 0.01
# rate = 0.0001
#
# rms1 = []
# for n in range(10):
#     rms1 += [optimize_peak_height(model, grad_rate, rate)]
#
# plt.scatter(list(range(len(rms1))),rms1)
# plt.show()
# plot_spectra(model, plot_sim=True)
# =============================================================================
#%%
# =============================================================================
#
# rms2 = []
# for n in range(10):
#     rms2 += [optimize_factor(model)]
#
# plt.scatter(list(range(len(rms2))),rms2)
# plt.show()
# plot_spectra(model, plot_sim=True)
# =============================================================================

#%%
# =============================================================================
# rms = []
# components  = []
# =============================================================================
#%%
# =============================================================================
# plt.plot(model.scattered_spectrum.lineshape)#[:-100])
# I = 25
# learning_rate = 0.01
# for i in range(I):
#     print(i)
#     model, rms_1 = iterate_once(model, learning_rate)
#     rms += [rms_1]
#     components += [[c.__dict__ for c in model.scattering_medium.scatterer.loss_function.components]]
#
#
# plt.plot(model.simulated_spectrum.lineshape)#[:-100])
# plt.show()
# plt.scatter(list(range(len(rms))),rms)
# plt.show()
# =============================================================================


#%%
#plt.scatter(list(range(len(rms))),rms)

#%%
#L = model.scattering_medium.scatterer.loss_function
#component_variables = [comp.__dict__ for comp in L.components]


#model.updateScatterersDict()
#model.saveScatterers('optimizedHe')

#%%
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

model, gradients, rms_hist = optimize_spsa(
    model,
    limits=limits,
    niter=15,
    )
# =============================================================================
#     a=1.0,
#     alpha=0.602,
#     c=1.0,
#     gamma=0.101)
# =============================================================================