import auxfuncs as aux 
import plots
import math
import re
import numpy
import matplotlib.pyplot as plt

def fluorescence_exposure(phores,setup_pars,exc_pars,rng_seed=42):
    """
    Simulates a single exposure of fluorescence microscopy

    Parameters
    ----------
    phores : 2D array
        types and positions of fluorophores
    setup_pars : 1D array
        parameters of optical setup and sample medium
    exc_pars : 1D array
        parameters of excitation illumination
    rng_seed : int
        seed for random number generators
    """

    n = setup_pars[0]
    NA = setup_pars[1]
    exptime = setup_pars[2]
    detector_qeff = setup_pars[3]
    pixel_size = setup_pars[4]
    field_size = setup_pars[5]
    w0 = exc_pars[0]
    wavelength = exc_pars[1]
    Pexc = exc_pars[2]
    I0 = exc_pars[3]

    #displays fluorophore distribution
    plots.display_2D_fluorophore_field(phores,w0,field_size,Pexc,wavelength)

    intensities = aux.field_add_illumination_intensities(phores, n, wavelength, w0, I0)
    print("Max intensity: %e" % (numpy.amax(intensities)))
    
    #imports the relevant absorption spectra and calculates cross-section at excitation wavelength
    xsections = aux.get_all_xsections(phores,wavelength)

    #calculates the mean numbers of absorbed photons per exposure  
    incident_photons = aux.all_incident(phores, exptime, wavelength, xsections, intensities)

    max_photons = numpy.amax(incident_photons)
    min_photons = numpy.amin(incident_photons)
    print("Maximum/minimum number of incident photons per fluorophore:\n %f / %f" % (max_photons,min_photons))

    #imports filter spectrum
    filter_spectrum = aux.get_filter_spectrum("test_filter")

    #calculates the number of photons the detector collects from each fluorophore 
    photon_counts = aux.calculate_single_image(phores, incident_photons, filter_spectrum, NA, n, detector_qeff, rng_seed)
    plots.display_detected_photon_counts(phores,w0,photon_counts)

    #simulates a finite detector resolution
    hist,_,_ = aux.pixel_binning(phores,photon_counts,w0*field_size,pixel_size)
    plots.display_detected_image(hist)
    profile = aux.radial_signal_profile(phores,photon_counts)

    return photon_counts,hist,profile


def CW_STED_beam_fluorescence_exposure_comparison(phores,setup_pars,exc_pars,STED_pars,rng_seed=42):
    """
    Simulates a single exposure of fluorescence microscopy with and without STED illumination

    Parameters
    ----------
    phores : 2D array
        types and positions of fluorophores
    setup_pars : 1D array
        parameters of optical setup and sample medium
    exc_pars : 1D array
        parameters of excitation illumination
    STED_pars : 1D array
        parameters of depletion beam 
    rng_seed : int
        seed for random number generators
    """

    n = setup_pars[0]
    NA = setup_pars[1]
    exptime = setup_pars[2]
    detector_qeff = setup_pars[3]
    pixel_size = setup_pars[4]
    field_size = setup_pars[5]
    w0 = exc_pars[0]
    wavelength = exc_pars[1]
    Pexc = exc_pars[2]
    I0 = exc_pars[3]
    STEDwavelength = STED_pars[0]
    PSTED = STED_pars[1]

    #displays fluorophore distribution
    plots.display_2D_fluorophore_field(phores,w0,field_size,Pexc,wavelength)

    intensities = aux.field_add_illumination_intensities(phores, n, wavelength, w0, I0)

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

    FWHM,popt = aux.FWHM_calculator_Gaussian(profile,w0)
    FWHM_STED,popt_STED = aux.FWHM_calculator_lin(STED_profile,w0)

    plots.compare_profiles(profile,STED_profile,field_size*w0,res1=FWHM,res2=FWHM_STED,popt1=popt,popt2=popt_STED)

    return photon_counts,STED_photon_counts,hist,STEDhist,profile,STED_profile

def define_emission_sampler(phoretype,filename="dye_spectra/Laplace_PDFs.dat",Npts=10,maxiter=100):
    """
    Runs the distribution optimization process on a spectrum and writes the results to file

    Parameters
    ----------
    phoretype : int
        fluorophore identifier
    filename : str
        name of file containing distribution parameters
    Npts : int
        number of interpolation points per nm
    maxiter : int
        maximum number of allowed iterations
    """
    params = aux.optimize_distribution(phoretype,Npts,maxiter)
    plots.display_sampling_pdf(phoretype,params)
    quality = input("Good fit? y/n: ")

    if (quality == 'y'):
        phorename,_ = aux.read_properties(phoretype)
        stringmatch = r"\b%s\s%s\b.*" % (str(phoretype).zfill(3),phorename)
        stringreplace = r"%s\t%s\t%.8f\t%.8f\t%.8f\t%.8f" % (str(phoretype).zfill(3),phorename,params[0],params[1],params[2],params[3])
        txt = "".join(open(filename).readlines())
        txt = re.sub(stringmatch, stringreplace, txt)
        f = open(filename, 'w')
        f.write(txt)
        f.close()

        print("Fit accepted. Parameters written to file.")
        return None
    else:
        print("Fit not accepted. Parameters not written to file.")
        return None