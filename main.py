from collections import namedtuple
import math

import numpy

import auxfuncs as aux
import plots
import routines


#fluorophore parameters
density = 5e13
phoretype = 1
#sample and setup parameters
n = 1.003
NA = 0.34
exptime = 0.02
detector_qeff = 0.9
pixel_size = 0.3e-6
field_size = 1.5
#excitation beam parameters
w0 = 1e-6
wavelength = 490e-9
Pexc = 1e-5
I0 = 2*Pexc/(math.pi*(w0**2))
#STED beam parameters
STEDwavelength = 592e-9
PSTED = 3e6*Pexc

#other parameters
rng_seed = 17

#generates positions of fluorescing molecules and illumination intensities at those points
phores = aux.generate_fluorophore_field(w0, density, phoretype, seed=42, latmultiplier=field_size, axmultiplier=field_size,volume=False)
#plots.display_3D_fluorophore_field(phores, w0, field_size, field_size, Pexc, wavelength)

#packages some parameters for brevity
Setup = namedtuple("Setup","n NA exptime detector_qeff pixel_size field_size")
Excitation = namedtuple("Excitation","w0 wavelength Pexc I0")
STED = namedtuple("STED","STEDwavelength PSTED")
setup_pars = Setup(n,NA,exptime,detector_qeff,pixel_size,field_size)
exc_pars = Excitation(w0,wavelength,Pexc,I0)
STED_pars = STED(STEDwavelength,PSTED)

#imports filter spectrum
filter_spectrum = aux.get_filter_spectrum("test_filter")

#simulates a single exposure of fluorescence microscopy with and without STED illumination
routines.CW_STED_beam_fluorescence_exposure_comparison(phores,setup_pars,exc_pars,STED_pars,filter_spectrum,rng_seed)