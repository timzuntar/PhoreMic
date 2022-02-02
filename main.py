import auxfuncs as aux
import plots
import numpy
import matplotlib.pyplot as plt

w0 = 1e-6
density = 5e13
wavelength = 490e-9
I0 = 5e-2
n = 1.003
NA = 0.34
exptime = 0.001
detector_qeff = 0.9
phoretype = 1
pixel_size = 0.3e-6

field_size = 1.5

phores = aux.generate_fluorophore_field(w0, density, phoretype, seed=42, sizemultiplier=field_size)
intensities = aux.field_add_illumination_intensities(phores, n, wavelength, w0, I0)

plots.display_2D_fluorophore_field(phores,w0)

xsection = aux.get_absorption_xsection(phoretype,wavelength*1e9)

max_photons = aux.incident_photons_per_exposure(exptime, wavelength*1e9, xsection, numpy.amax(intensities))
print("Max. number of incident photons per fluorophore is %f." % (max_photons))

rng_seed = 17
filter_spectrum = aux.get_filter_spectrum("test_filter")

photon_counts = aux.calculate_single_image(phores, intensities, filter_spectrum, NA, n, wavelength*1e9, exptime, detector_qeff, rng_seed)

plots.display_detected_photon_counts(phores,w0,photon_counts)

#placeholder; right now all types of fluorophores are output on the same histogram, potentially reducing performance 
hist,_,_ = aux.pixel_binning(phores,photon_counts,w0*field_size,pixel_size)

plots.display_detected_image(hist)