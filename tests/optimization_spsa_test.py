# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 13:25:24 2021

@author: pielsticker
"""

"""
A class to implement Simultaneous Perturbation Stochastic Approximation.
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import noisyopt
#https://noisyopt.readthedocs.io/en/latest/api.html#minimizespsa

from converters.data_converter import DataConverter
from model.model import Model
from model.base_model import MeasuredSpectrum
#%%
# First set-up the model.
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
    model.loaded_spectra += [
        MeasuredSpectrum(data["data"]["x"], data["data"]["y0"])
    ]
# =============================================================================
#     if data_dict["spectrum_type"] == "Ag3d":
#         if "Ag in vacuum" in data_dict["group_name"]:
#             model.scattered_spectrum = model.loaded_spectra[i]
#         elif "Ag in 4 mbar He" in data_dict["group_name"]:
#             model.unscattered_spectrum = model.loaded_spectra[i]
# =============================================================================

model.scattered_spectrum = model.loaded_spectra[0]
model.unscattered_spectrum = model.loaded_spectra[1]

scatterer_path = os.path.join(
    *[directory, "data", "scatterers2.json"]
)

model.loadScatterers(scatterer_path)
model.setCurrentScatterer('Shirley')
model.scattering_medium.setPressure(4)
model.scattering_medium.scatterer.inelastic_xsect = 0.003
model.scattering_medium.scatterer.norm_factor = 1
model.algorithm_id = 0

model.scattering_medium.scatterer.loss_function.reBuild()
model.scatterSpectrum()

plt.plot(model.scattered_spectrum.x,model.scattered_spectrum.lineshape)
plt.plot(model.unscattered_spectrum.x,model.unscattered_spectrum.lineshape)

#%%

def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())

def _scatter_with_variables(model, param_dict):
    plt.plot(model.simulated_spectrum.x,model.simulated_spectrum.lineshape)
    L = model.scattering_medium.scatterer.loss_function
    component_variables = [comp.__dict__ for comp in L.components]
    for idx, comp in enumerate(component_variables):
        print(comp)
        for key, value in comp.items():
            L.editComponent(idx, {key: param_dict[key]})
            model.scattering_medium.scatterer.loss_function = L
            model.scatterSpectrum()
        print(comp)
    plt.plot(model.simulated_spectrum.x,model.simulated_spectrum.lineshape)
    plt.show()

    return model

def objective_fn(param_array, model):
    param_dict = {
    "width": param_array[0],
    "position": param_array[1],
    "intensity": param_array[2],
    "B": param_array[3],
    "D": param_array[4],
    "Eg": param_array[5],
    "edge": param_array[6],
    "exponent": param_array[7],
    "fermi_width": param_array[8]
    }

    model = _scatter_with_variables(model, param_dict)
    L = model.scattering_medium.scatterer.loss_function
    component_variables = [comp.__dict__ for comp in L.components]
    predictions = model.simulated_spectrum.lineshape
    targets = model.scattered_spectrum.lineshape

    return model, rmse(predictions, targets)

initial_param_dict = {
    "width": 2.0,
    "position": 5.0,
    "intensity": 1,
    "B": 10.0,
    "D": 10.0,
    "Eg": 0.1,
    "edge": 20.0,
    "exponent": 0.1,
    "fermi_width": 20.0
    }

bounds_dict = {
    "width": np.array([1.0,20.0]),
    "position": np.array([0,30.0]),
    "intensity": np.array([1.0,20.0]),
    "B": np.array([5.0,20.0]),
    "D": np.array([5.0,20.0]),
    "Eg": np.array([1.0,20.0]),
    "edge": np.array([0.0,50.0]),
    "exponent": np.array([0.05,0.5]),
    "fermi_width": np.array([1.0,50.0])
    }
params_in = np.array(list(initial_param_dict.values()))
bounds = np.array(list(bounds_dict.values()))


(model2, rmse) = objective_fn(params_in, model)

plt.plot(model.scattered_spectrum.x,model.scattered_spectrum.lineshape)
plt.plot(model.simulated_spectrum.x,model.simulated_spectrum.lineshape)


noisyopt.minimizeSPSA(
    func=objective_fn,
    x0=params_in,
    args=(model),
    bounds=None,
    niter=100,
    paired=True,
    a=0.5,
    c=1.0,
    disp=True)


#%%
param_ranges = {
    "width": np.arange(0.01,10,0.005),
    "position": np.arange(0,1000,0.05),
    "intensity": np.arange(0,1000,1),
    "B": np.arange(0,1E+9,100),
    "D": np.arange(0,1E+9,100),
    "Eg": np.arange(0,0,0.0000000000001),
    "edge": np.arange(0.0,25.0,0.5),
    "exponent": np.arange(0.0,25.0,0.5),
    "fermi_width": np.arange(0.0,0.005,0.0001),
    }

component_variables = [comp.__dict__ for comp in model.scattering_medium.scatterer.loss_function.components]
numbers = [-1,1]
step_sizes = {
    'width':0.005,
    'position':0.05,
    'intensity':1,
    'B':100,
    'C':100,
    'D':10,
    'Eg':0,
    "edge":0.2,
    "exponent":0.005,
    "fermi_width": 0.0001
    }
limits = {
    'width':(0.01,10),
    'position':(0,1000),
    'intensity':(0,1000),
    'B':(0,1E+9),
    'C':(0,1E+9),
    'D':(0.01,1000),
    'Eg':(0,0),
    "edge":(0.0,25.0),
    "exponent":(0.01,0.07),
    "fermi_width": (0.005,0.005)
    }
gradients = []
for idx, comp in enumerate(component_variables):
    gradients += [{}]
    for key, value in comp.items():
        pass

#%%
import noisyopt
import numpy as np

p_in = [0.1, -2.6, -1.5]
noise_var = 0.3
fitfunc = lambda p, x: p[0]*x*x + p[1]*x + p[2]
errfunc = lambda p, x, y, noise_var: np.sum ( (fitfunc( p, x ) - y)**2/ \
    noise_var**2 )
# The following error function can be used by leastsq in scipy.optimize
#errfunc2 = lambda p, x, y: fitfunc( p, x ) - y
# make some data
x_arr = np.arange(100) * 0.3
obs = p_in[0] * x_arr**2 + p_in[1] * x_arr + p_in[2]
np.random.seed(76523654)
noise = np.random.normal(size=100) * noise_var  # add some noise to the obs
obs += noise

res = noisyopt.minimizeSPSA(
    func=errfunc,
    x0=p_in,
    args=(x_arr, obs, noise_var),
    bounds=None,
    niter=100,
    paired=False,
    a=1.0,
    c=1.0,
    disp=False,
    callback=None)


p2 = np.array([0.5, 2.6, 1.5])
plt.plot(x_arr,
         np.array([errfunc(p_in, x, obs, noise_var) for x in x_arr]))
#plt.show()
plt.plot(x_arr,
         np.array([errfunc(p2, x, obs, noise_var) for x in x_arr]))
plt.show()

plt.plot(x_arr,
         np.array([fitfunc(p_in, x) for x in x_arr]))
#plt.show()
plt.plot(x_arr,
         np.array([fitfunc(p2, x) for x in x_arr]))
y = res.x
plt.plot(x_arr,
         np.array([fitfunc(y, x) for x in x_arr]))
plt.show()


#%%
def obj(x):
    return (x**2).sum() + 0.1*np.random.randn()

bounds = [[-3.0, 3.0], [0.5, 5.0]]
x0 = np.array([-2.0, 2.0])

res = noisyopt.minimizeSPSA(
    obj,
    bounds=bounds,
    x0=x0, niter=1000,
    paired=False)

x = np.arange(-5,5,0.001)

plt.plot(x,
         np.array([obj(x) for x in x]))

a = np.argmin(np.array([obj(x) for x in x]))
print(x[a])



