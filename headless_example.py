# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 16:40:20 2021

@author: Mark
"""
from converters.data_converter import DataConverter
from model.base_model import MeasuredSpectrum, ScatteringMedium
from model.algorithms import algorithm5
from model.model import Model
from model.base_model import VacuumExcitation
import matplotlib.pyplot as plt
import os
import sys
import json

#%%
class SpectrumFigure:
    def __init__(self, model, params):
        self.unscattered_spectrum = model.unscattered_spectrum
        self.scattered_spectrum = model.scattered_spectrum
        self.simulated_spectrum = model.simulated_spectrum
        self.loss_function = model.scattering_medium.scatterer.loss_function
        self.params = params
        self.shirley_params = params["shirley_params"]

        self.fig, self.axs = plt.subplots(1, 2)

        plt.subplots_adjust(
            left=0.125, bottom=0.5, right=4.8, top=2, wspace=0.2, hspace=0.2,
        )

    def plot(self):
        self.axs[0].plot(
            self.unscattered_spectrum.x, self.unscattered_spectrum.lineshape
        )
        self.axs[0].plot(
            self.scattered_spectrum.x, self.scattered_spectrum.lineshape
        )
        self.axs[0].plot(
            self.simulated_spectrum.x, self.simulated_spectrum.lineshape
        )
        ax0_text = (
            f"No. of Iterations: {self.params['n_iter']} \n"
            f"Algorithm Option: {self.params['algorithm_option']} \n\n"
            f"Norm Factor: {self.params['norm_factor']} \n"
        )
        self.axs[0].text(
            0.1,
            0.5,
            ax0_text,
            horizontalalignment="left",
            verticalalignment="top",
            transform=self.axs[0].transAxes,
            fontsize=8,
        )
        self.axs[0].legend(["unscattered", "simulated"])

        self.axs[1].plot(self.loss_function.x, self.loss_function.lineshape)
        self.axs[1].set_xlim(0, 200)
        ax1_text = "Loss function params:\n\n"
        for key, item in self.shirley_params.items():
            ax1_text += (
                f"{key} \n"
                f"   Edge : {item['edge']} \n"
                f"   Fermi Width : {item['fermi_width']} \n"
                f"   Intensity : {item['intensity']} \n"
                f"   Exponent : {item['exponent']} \n\n"
            )
        self.axs[1].text(
            0.3,
            0.9,
            ax1_text,
            horizontalalignment="left",
            verticalalignment="top",
            transform=self.axs[1].transAxes,
            fontsize=7,
        )
        self.axs[1].legend(["Shirley Loss Function"])

        plt.tight_layout()

        return self.fig, self.axs


#%%
abspath = os.path.abspath(__file__)
directory = os.path.dirname(abspath)

directory = os.getcwd()

spectrum_filename = "SiO2 O1s"

filename = os.path.join(
    *[directory, "data", "Shirley", spectrum_filename + ".vms"]
)

converter = DataConverter()
converter.load(filename)

scattered_index = 0
unscattered_index = 1

model = Model()

for i in [unscattered_index, scattered_index]:
    data = converter.data[i]
    model.loaded_spectra += [
        MeasuredSpectrum(data["data"]["x"], data["data"]["y0"])
    ]
model.unscattered_spectrum = model.loaded_spectra[0]
model.scattered_spectrum = model.loaded_spectra[1]

for spectrum in model.loaded_spectra:
    plt.plot(spectrum.x, spectrum.lineshape)
plt.legend(["unscattered", "scattered"])
plt.show()

filename2 = os.path.join(directory, r"data\scatterers2.json")
model.loadScatterers(filename2)


#%%
model.scattering_medium.d_through_gas = 10
model.scattering_medium.pressure = 40
model.setCurrentScatterer("Shirley")
model.changeAlgorithm(1)

test_params = {
    "norm_factor": 0.88,
    "n_iter": 50,
    "algorithm_option": "bulk",
    "shirley_params": {
        "loaded": {
            "edge": 11.5,
            "exponent": 0.065,
            "fermi_width": 20,
            "intensity": 1,
        }
    },
}

model.scattering_medium.scatterer.norm_factor = test_params["norm_factor"]
model.calculation.n_iter = test_params["n_iter"]
model.algorithm_option = test_params["algorithm_option"]
model.scattering_medium.scatterer.loss_function.editComponent(
    0, test_params["shirley_params"]["loaded"]
)

vac_exc_params = {
    "Additional Vac. Excitation 1": {
        "edge": 11.5,
        "exponent": 0.065,
        "fermi_width": 20,
        "intensity": 25,
    },
    "Additional Vac. Excitation 2": {
        "edge": 11.5,
        "exponent": 0.065,
        "fermi_width": 20,
        "intensity": 500,
    },
}
legend = ["Loaded Shirley"]
for key, d in vac_exc_params.items():
    vac_exc = VacuumExcitation(
        edge=d["edge"],
        fermi_width=d["fermi_width"],
        intensity=d["intensity"],
        exponent=d["exponent"],
    )
    model.scattering_medium.scatterer.loss_function.addComponent(vac_exc)
    legend.append(f"{key}")

model.scattering_medium.scatterer.loss_function.buildLine()
test_params["shirley_params"].update(vac_exc_params)


lf = model.scattering_medium.scatterer.loss_function
for component in lf.components:
    y = component.function(lf.x)
    plt.plot(lf.x, y)
plt.legend(legend)
plt.title("Loss function components")
plt.show()

model.scatterSpectrum()

spectrum_fig = SpectrumFigure(model=model, params=test_params)
fig, axs = spectrum_fig.plot()
plt.plot()


#%%
output_folder = r"C:\Users\pielsticker\Lukas\MPI-CEC\Publications\Shirley paper\test figures"
filenumber = len(
    [f for f in os.listdir(output_folder) if f.startswith(spectrum_filename)]
)
fig_filepath = os.path.join(
    output_folder, f"{spectrum_filename}_{filenumber}.png"
)
fig.savefig(fig_filepath)
