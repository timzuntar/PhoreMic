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
PSTED = 1e6*Pexc

#other parameters
rng_seed = 17

#generates positions of fluorescing molecules and illumination intensities at those points
phores = aux.generate_fluorophore_field(w0, density, phoretype, seed=42, latmultiplier=field_size, axmultiplier=field_size)
intensities = aux.field_add_illumination_intensities(phores, n, wavelength, w0, I0)

#displays fluorophore distribution
plots.display_2D_fluorophore_field(phores,w0,field_size,Pexc,wavelength)

#imports the relevant absorption spectra and calculates cross-section at excitation wavelength
xsections = aux.get_all_xsections(phores,wavelength)
STEDxsections = aux.get_all_xsections(phores,STEDwavelength)
#as well as the saturation intensities
Isats = aux.STED_get_all_Isat(phores,STEDxsections,STEDwavelength)

#calculates the mean numbers of absorbed photons per exposure  
incident_photons = aux.all_incident(phores, exptime, wavelength, xsections, intensities)

#the following needs to be calculated for STED illumination
#excitation rates of the main beam
exc_rates = incident_photons/exptime
#intensities of STED beam
STED_intensities = aux.field_STED_illumination_intensities(phores, STEDwavelength, PSTED, NA)

STED_incident_photons = aux.STED_all_incident(phores, intensities, STED_intensities, exptime, wavelength, exc_rates, xsections, Isats)

print("Max intensity\nno STED: %e STED: %e" % (numpy.amax(intensities),numpy.amax(STED_intensities)))

max_photons = numpy.amax(incident_photons)
min_photons = numpy.amin(incident_photons)
STED_max_photons = numpy.amax(STED_incident_photons)
STED_min_photons = numpy.amin(STED_incident_photons)

print("Maximum number of incident photons per fluorophore:\nwithout STED: %f with STED: %f" % (max_photons,STED_max_photons))
print("Minimum number of incident photons per fluorophore:\nwithout STED: %f with STED: %f" % (min_photons,STED_min_photons))

#imports filter spectrum
filter_spectrum = aux.get_filter_spectrum("test_filter")

#calculates the number of photons the detector collects from each fluorophore 
photon_counts = aux.calculate_single_image(phores, incident_photons, filter_spectrum, NA, n, detector_qeff, rng_seed)
STED_photon_counts = aux.calculate_single_image(phores, STED_incident_photons, filter_spectrum, NA, n, detector_qeff, rng_seed)

plots.display_photon_counts_side_by_side(phores,Pexc,wavelength,PSTED,STEDwavelength,w0,photon_counts,STEDwavelength/(NA*2.0),STED_photon_counts,field_size,alt_type="STED")

#simulates a finite detector resolution
#placeholder; right now all types of fluorophores are output on the same histogram, potentially reducing performance 
hist,_,_ = aux.pixel_binning(phores,photon_counts,w0*field_size,pixel_size)
STEDhist,_,_ = aux.pixel_binning(phores,STED_photon_counts,w0*field_size,pixel_size)

plots.display_detected_images(pixel_size,hist,STEDhist,alt_type="STED")

profile = aux.radial_signal_profile(phores,photon_counts)
STED_profile = aux.radial_signal_profile(phores,STED_photon_counts)

plots.compare_profiles(profile,STED_profile,field_size*w0)
