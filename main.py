import auxfuncs as aux
import plots
import math
import numpy
import matplotlib.pyplot as plt

#fluorophore parameters
density = 5e13
phoretype = 1
#sample and setup parameters
n = 1.003
NA = 0.34
exptime = 0.001
detector_qeff = 0.9
pixel_size = 0.3e-6
field_size = 1.5
#excitation beam parameters
w0 = 1e-6
wavelength = 490e-9
Pexc = 1e-4
I0 = 2*Pexc/(math.pi*(w0**2))

phores = aux.generate_fluorophore_field(w0, density, phoretype, seed=42, latmultiplier=field_size, axmultiplier=field_size)
intensities = aux.field_add_illumination_intensities(phores, n, wavelength, w0, I0)

plots.display_2D_fluorophore_field(phores,w0)

xsections = aux.get_all_xsections(phores,wavelength)
print(xsections)

incident_photons = aux.all_incident(phores, exptime, wavelength, xsections, intensities)

max_photons = numpy.amax(incident_photons)
print("Max. number of incident photons per fluorophore is %f." % (max_photons))

rng_seed = 17
filter_spectrum = aux.get_filter_spectrum("test_filter")

photon_counts = aux.calculate_single_image(phores, incident_photons, filter_spectrum, NA, n, detector_qeff, rng_seed)

plots.display_detected_photon_counts(phores,w0,photon_counts)

#placeholder; right now all types of fluorophores are output on the same histogram, potentially reducing performance 
hist,_,_ = aux.pixel_binning(phores,photon_counts,w0*field_size,pixel_size)

plots.display_detected_image(hist)