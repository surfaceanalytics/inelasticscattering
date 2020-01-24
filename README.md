The is a script to simulate the effect of gas molecules on the measured signal in X-ray photoelectron spectra (XPS).
The script takes measured XPS spectra, a synthesized energy loss function, and simulates the effect of gas molecules on the 
measured XPS spectrum.
When measuring XPS of a sample, the sample is usually held in vacuum. 
However, in the more recently developed form of XPS known as near-ambient pressure XPS, one can measure samples when they are
not in vacuum. That means there are gas molecules in the atmosphere surrounding the sample. When measuring the electrons that 
come off the sample due to photoemission, these electrons need to travel through the atmosphere before the reach the electron detector.
Along the way through the gas contining atmosphere, the electrons can collide with gas molecules and lose some of their kinetic 
energy. The end result is a signal in the measured XPS spectrum that comes from inelastically scattered electrons.

This script has been developed to be able to simulate that portion of the measured spectrum. 

One loads in two measured spectra: One with the scattered signal and opne without the scattered signal. 
One then also load (or builds) a synthetic lineshape which represents the energy loss function of the gas molecules.
Then one provides experimental parameters such as pressure and distance from the sample to the spectrometer.
The script then calculates the effect of the gas phase scattering. 

