import math
import glob
import pickle
import numpy
import scipy.stats
import scipy.special
import scipy.constants
import scipy.integrate
import scipy.interpolate
import scipy.optimize
import matplotlib.pyplot as plt


def gaussian_point_intensities(point_coords, n, wavelength, w0, I0):
    """
    Returns intensity of the TEM00 (Gaussian) beam at a point in space.

    Parameters
    ----------
    point_coords : 2D array
        x,y, and z positions relative to center of beam waist [m]
    n : float
        refractive index of propagation medium
    wavelength : float
        laser wavelength [m]
    w0 : float
        beam diameter at waist [m]
    I0 : float
        intensity in center of waist [W/m^2]
    """
    boolmask = numpy.isclose(point_coords[:,2], 0.0, rtol=1e-05, atol=1e-12, equal_nan=False)   # checks if point is very close (<1 pm) to zero plane
    zero_indices = numpy.where(boolmask)[0]
    nonzero_indices = numpy.where(numpy.logical_not(boolmask))[0]
    intensities = numpy.empty(len(boolmask),dtype=float)    # this saves a little bit of performance as width doesn't need to be computed
    intensities[zero_indices] = I0*numpy.exp(-2*(numpy.square(point_coords[zero_indices,0]) + numpy.square(point_coords[zero_indices,1]))/(w0*w0))
    beam_waists = w0 * numpy.sqrt(1 + (point_coords[nonzero_indices,2]*wavelength/(numpy.pi*n*w0*w0))**2)
    intensities[nonzero_indices] = I0*((w0/beam_waists)**2)*numpy.exp(-2*(numpy.square(point_coords[nonzero_indices,0]) + numpy.square(point_coords[nonzero_indices,1]))/numpy.square(beam_waists))
    return intensities
        
def STED_2D_approx_point_intensity(point_coords, STEDwavelength, P, NA):
    """
    Intensity of depletion beam using the standing wave approximation, with a simple analytical solution
    NOTE: out-of-focus beam spread is approximated with characteristics of a Gaussian beam with its waist diameter at the annulus maximum.

    Parameters
    ----------
    point_coords : 2D array
        x,y, and z positions relative to center of beam waist [m]
    STEDwavelength : float
        wavelength of depletion beam [m]
    P : float
        total beam power [W]
    NA : float
        numerical aperture of depletion beam
    """
    radii = numpy.sqrt(numpy.square(point_coords[:,0]) + numpy.square(point_coords[:,1]))
    boolmask = numpy.isclose(point_coords[:,2], 0.0, rtol=1e-05, atol=1e-12, equal_nan=False)   # checks if point is very close (<1 pm) to zero plane
    zero_inside_indices = numpy.where((boolmask) & (radii - STEDwavelength/NA < 0.0))[0]
    zero_outside_indices = numpy.where((boolmask) & (radii - STEDwavelength/NA >= 0.0))[0]
    intensities = numpy.empty(len(radii),dtype=float)

    intensities[zero_inside_indices] = ((2*P/numpy.pi) * numpy.square(NA/STEDwavelength)) * numpy.square(numpy.sin(numpy.pi*NA*radii[zero_inside_indices]/STEDwavelength))
    intensities[zero_outside_indices] = 0.0

    w0 = STEDwavelength/(2*NA)
    z0 = numpy.pi*w0*w0/STEDwavelength
    widening = numpy.sqrt(1+numpy.square(point_coords[:,2]/z0))
    nonzero_inside_indices = numpy.where((numpy.logical_not(boolmask)) & (radii-2*w0*widening < 0.0))[0]
    nonzero_outside_indices = numpy.where((numpy.logical_not(boolmask)) & (radii-2*w0*widening >= 0.0))[0]
    intensities[nonzero_inside_indices] = (P/(2*math.pi*numpy.square(w0*widening[nonzero_inside_indices]))) * numpy.square(numpy.sin(numpy.pi*NA*radii[nonzero_inside_indices]/(STEDwavelength*widening[nonzero_inside_indices])))  #TODO: double-check!
    intensities[nonzero_outside_indices] = 0.0
    return intensities

def STED_2D_point_intensity(point_coords, STEDwavelength, NA_eff):
    """
    Determines the intensity of depletion beam at chosen point with prediction for a coherent plane wave passing through a 2-pi vortex plate
    as described by Neupane et al. (2013)
    NOTE: does not take into account out-of-focus beam spread, should only be used for 2D distributions!

    Parameters
    ----------
    point_coords : list of float
        x,y, and z position relative to center of beam waist [m]
    STEDwavelength : float
        wavelength of depletion beam [m]
    NA_eff : float
        effective numerical aperture of depletion beam
    """
    k = 2*math.pi/STEDwavelength
    x = k*NA_eff*numpy.sqrt(numpy.square(point_coords[:,0])+numpy.square(point_coords[:,1]))

    C = (0.5*numpy.pi*k*NA_eff)**2
    return (C/(numpy.square(x)))*numpy.square(scipy.special.jv(1,x)*scipy.special.struve(0,x)-scipy.special.jv(0,x)*scipy.special.struve(1,x))

def STED_2D_intensity_normalization(n, STEDwavelength, P, NA_eff, lossfactor=0.0, integrationwidth=10):
    """
    Determines the intensity of depletion beam from total incident power (assuming losses are constant) via numerical integration

    Parameters
    ----------
    n : float
        refractive index of propagation medium
    STEDwavelength : float
        wavelength of depletion beam [m]
    P : float
        total beam power [W]
    NA_eff : float
        effective numerical aperture of depletion beam
    lossfactor : float
        relative amount of beam power not incident on the sample
    integrationwidth : float
        maximum integration radius in multiples of peak intensity radius
    """
    peakradius = 2.45*STEDwavelength/(2*math.pi*NA_eff)
    samples = 5001
    rlist = numpy.linspace(1e-9, peakradius*integrationwidth, num=samples, endpoint=True)
    ilist = numpy.empty(samples)
    ilist = STED_2D_point_intensity(numpy.column_stack((rlist,numpy.zeros_like(rlist),numpy.zeros_like(rlist))),STEDwavelength,NA_eff)*rlist*2*numpy.pi

    integral = numpy.fabs(scipy.integrate.simpson(rlist,ilist))
    return P*(1.0-lossfactor)/integral

def STED_saturation_intensity(lifetime, STxsection, STwavelength):
    """
    Saturation intensity of STED beam

    Parameters
    ----------
    lifetime : float
        fluorescence lifetime of molecule [s]
    STxsection : float
        cross-section for stimulated emission [square meters]
    STwavelength : float
        wavelength of STED beam
    """

    ks1 = 1.0/(lifetime)
    return scipy.constants.c*scipy.constants.h*ks1/(STxsection*STwavelength)

def STED_get_all_Isat(phores, STxsections, STwavelength):
    """
    Calculates saturation intensities for all fluorophore types
    
    Parameters
    ----------
    phores : 2D array
        types and positions of fluorophores
    STxsections : 1D array
        cross-section for stimulated emission [square meters]
    STEDwavelength : float
        wavelength of depletion beam [m]
    """
    phoretypes = numpy.unique(phores[:,0]).astype(int)
    Isats = numpy.empty(len(phoretypes))

    for i,type in enumerate(phoretypes):
        _,props = read_properties(type)
        lifetime = props[3]
        Isats[i] = STED_saturation_intensity(lifetime,STxsections[i],STwavelength)
    return Isats

def STED_effective_saturation(I, Isat, ks1, vibrelaxrate):
    """
    Effective saturation factor (gamma)
    
    Parameters
    ----------
    I : float
        intensity of STED beam at point [W/m^2]
    Isat : float
        saturation intensity [W/m^2]
    ks1 : float
        rate of spontaneous fluorescence (inverse of lifetime) [1/s]
    vibrelaxrate : float
        rate of vibrational relaxation into ground state [1/s]
    """
    xi = I/Isat
    return (xi*vibrelaxrate)/(xi*ks1 + vibrelaxrate)

def STED_CW_rates(I, Isat, kex, ks1, vibrelaxrate, intersystem=0, kt1=1.0):
    """
    Probability for spontaneous decay and its rate for continuous STED illumination 
    
    Parameters
    ----------
    I : float
        intensity of STED beam at point [W/m^2]
    Isat : float
        saturation intensity [W/m^2]
    kex : float
        excitation rate (illumination-dependent) [1/s]
    ks1 : 
        rate of spontaneous fluorescence (inverse of lifetime) [1/s]
    vibrelaxrate : float
        rate of vibrational relaxation into ground state [1/s]
    intersystem : float
        intersystem crossing yield (probability for decay to triplet state)
    kt1 : float
        decay rate of triplet state [1/s]
    """
    gamma = STED_effective_saturation(I,Isat,ks1,vibrelaxrate)
    kexprime = kex*vibrelaxrate/(kex+vibrelaxrate)
    spontaneous_rate = 1.0/((1+gamma)/kexprime + 1.0/ks1 + intersystem/kt1)
    spontaneous_probability = (1.0/kexprime + 1.0/ks1 + intersystem/kt1)*spontaneous_rate
    return spontaneous_rate,spontaneous_probability  

def generate_fluorophore_field(w0, density, phoretype, existing=None, seed=42, volume=False, latmultiplier=5, axmultiplier=1, maxnum=1e7):
    """
    Randomly generates positions of fluorescing molecules

    Parameters
    ----------
    w0 : float
        beam waist diameter [m]
    density : float
        number density of fluorescing molecules [per square meter in 2D / per cubic meter in 3D]
    phoretype : int
        fluorophore identifier
    existing : 2D array
        previous array; output of function consists of new entries concatenated with previous ones
    seed : int
        seed to random number generator
    volume : bool
        controls whether a 3D or 2D fluorophore field is generated
    latmultiplier : float
        width and height of created field in multiples of beam waist diameter
    axmultiplier : float
        as latmultiplier, but along the z-axis
    maxnum : int
        maximum number of generated points
    """
    if (volume == False):
        phorenum = int((w0*2*latmultiplier)**2 * density)
    else:
        phorenum = int((w0*2*latmultiplier)**2 * (w0*2*axmultiplier) * density)

    print("%d points initialized." % (phorenum))
    if (phorenum > maxnum):
        response = input("Requested number of molecules %.2e exceeds the set maximum (%.2e).\nc = continue\nr = reduce density\nother keys = abort" % (phorenum, maxnum))
        if (response == "c"):
            pass
        elif (response == "r"):
            phorenum = maxnum
        else:
            quit()

    setseed = numpy.random.SeedSequence(seed)
    rng = numpy.random.Generator(numpy.random.PCG64(setseed))
    xpositions = rng.uniform(-w0*latmultiplier,w0*latmultiplier,phorenum)
    ypositions = rng.uniform(-w0*latmultiplier,w0*latmultiplier,phorenum)
    if (volume == False):
        zpositions = numpy.full(phorenum, 0.0, dtype=float)
    else:
        zpositions = rng.uniform(-w0*axmultiplier,w0*axmultiplier,phorenum)
    types = numpy.full(phorenum, phoretype, dtype=int)
    phores = numpy.concatenate((types, xpositions, ypositions,zpositions)).reshape((-1, 4), order='F')
    if (isinstance(existing,numpy.ndarray)):
        return numpy.vstack((existing,phores))
    else:
        return phores

def get_all_xsections(phores,wavelength):
    """
    Returns absorption cross sections for all fluorophore types.    

    Parameters
    ----------
    phores : 2D array
        types and positions of fluorophores
    wavelength : float
        illumination wavelength [m]
    """
    phoretypes = numpy.unique(phores[:,0]).astype(int)
    xsections = numpy.empty(len(phoretypes))
    for i,type in enumerate(phoretypes):
        xsections[i] = get_absorption_xsection(type,wavelength)
    return xsections

def illumination_fraction(NA, n):
    """
    Fraction of solid angle over which emitted light is collected. Assumes source to be exactly in focus,
    but the errors caused by that are negligible.

    Parameters
    ----------
    NA : float
        numerical aperture
    n : float
        refractive index of propagation medium
    """
    theta = math.asin(NA/n)
    fraction = (math.sin(theta/2.0))**2
    return fraction

def incident_photons_per_exposure(exptime, wavelength, xsection, intensity):
    """
    Expected number of photons absorbed by the fluorophore during exposure time

    Parameters
    ----------
    exptime : float
        exposure time per image, assumed to be identical to illumination time [s]
    wavelength : float
        laser wavelength [nm]
    xsection : float
        absorption cross section for given wavelength [m^2]
    intensity : float
        local illumination intensity [W/m^2]
    """
    return exptime*xsection*intensity*wavelength/(scipy.constants.c*scipy.constants.h)

def all_incident(phores, exptime, wavelength, xsections, intensities):
    """
    Expected absorbed photon numbers for all fluorophores

    Parameters
    ----------
    phores : 2D array
        types and positions of fluorophores
    exptime : float
        exposure time per image, assumed to be identical to illumination time [s]
    wavelength : float
        laser wavelength [nm]
    xsections : 1D array
        absorption cross sections of fluorophore species at given wavelength [m^2]
    intensities : 1D array
        local illumination intensities for all molecules [W/m^2]
    """
    phorenum = numpy.shape(phores)[0]
    incident_counts = numpy.empty(phorenum,dtype=float)
    phoretypes = numpy.unique(phores[:,0]).astype(int)

    for t,phoretype in enumerate(phoretypes):
        xsection = xsections[t]
        where_phores = numpy.where(phores[:,0] == phoretype)[0]
        incident_counts[where_phores] = incident_photons_per_exposure(exptime,wavelength,xsection,intensities[where_phores])
    return incident_counts

def STED_all_incident(phores, intensities, STEDintensities, exptime, wavelength, exc_rates, xsections, Isats):
    """
    Expected numbers of absorbed photons that undergo spontaneous decay under STED illumination  

    Parameters
    ----------
    phores : 2D array
        types and positions of fluorophores
    intensities : 1D array
        local illumination intensities for all molecules [W/m^2]
    STEDintensities : 1D array
        local illumination intensities of depletion beam [W/m^2]
    exptime : float
        exposure time [s]
    wavelength : float
        illumination wavelength [nm]
    exc_rates : 1D array
        excitation rates of all molecules (illumination-dependent) [1/s]
    xsections : 1D array
        absorption cross sections of fluorophore species at given wavelength [m^2]
    Isats : 1D array
        saturation intensities of fluorophore species [W/m^2]
    """
    phorenum = numpy.shape(phores)[0]
    incident_counts = numpy.empty(phorenum,dtype=float)
    phoretypes = numpy.unique(phores[:,0]).astype(int)

    for t,phoretype in enumerate(phoretypes):
        xsection = xsections[t]
        Isat = Isats[t]
        _,props = read_properties(phoretype)
        lifetime = props[3]
        vibrelaxrate,intersystem,kt1 = read_STED_properties(phoretype)
        ks1 = 1.0/(lifetime)

        where_phores = numpy.where(phores[:,0] == phoretype)[0]
        _,prob = STED_CW_rates(STEDintensities[where_phores], Isat, exc_rates[where_phores], ks1, vibrelaxrate, intersystem=intersystem, kt1=kt1)
        incident_counts[where_phores] = prob*incident_photons_per_exposure(exptime,wavelength,xsection,intensities[where_phores])
    return incident_counts

def interpolate_absorption_spectrum(filename, example_xsection, example_wavelength, show=False):
    """
    Interpolates the absorption spectrum from a file, scales it per-molecule and saves the object

    Parameters
    ----------
    filename : str
        name of file containing the spectrum to be interpolated
    example_xsection : float
        absorption cross section of fluorophore at given wavelength [m^2]
    example_wavelength : float
        corresponding wavelength [nm]
    show : bool
        controls whether interpolation curve is shown
    """
    data = numpy.loadtxt(filename,dtype=float,delimiter="\t")
    ftest = scipy.interpolate.interp1d(data[:,0],data[:,1])
    c = ftest(example_wavelength)
    multiplicative_factor = example_xsection/c
    f = scipy.interpolate.interp1d(data[:,0],data[:,1]*multiplicative_factor,bounds_error=False,fill_value=0.0)
    if (show == True):
        plt.plot(data[:,0],data[:,1]*multiplicative_factor, 'o', data[:,0], f(data[:,0]), '-')
        plt.show()

    pickle.dump(f,open(filename.replace(".txt", ".pkl"),"wb"))
    return True

def interpolate_emission_spectrum(filename, log=False, show=False):
    """
    Interpolates the emission spectrum from a file and saves the object

    Parameters
    ----------
    filename : str
        name of file containing the spectrum to be interpolated
    log : bool
        controls whether function is written as-is or converted to log(p) first 
    show : bool
        controls whether interpolation curve is shown
    """
    data = numpy.loadtxt(filename,dtype=float,delimiter="\t")
    if (log == True):
        f = scipy.interpolate.interp1d(data[:,0],numpy.log(data[:,1]))
    else:
        f = scipy.interpolate.interp1d(data[:,0],data[:,1],bounds_error=False,fill_value=0.0)

    if (show == True):
        plt.plot(data[:,0], f(data[:,0]), '-')
        plt.show()

    pickle.dump(f,open(filename.replace(".txt", ".pkl"),"wb"))
    return True

def get_absorption_xsection(phoretype, wavelength):
    """
    Estimates cross section for photon absorption at given wavelength

    Parameters
    ----------
    phoretype : int
        fluorophore identifier
    wavelength : float
        illumination wavelength [nm]
    """
    identifier = str(phoretype).zfill(3)
    filepath = "dye_spectra/" + identifier + "_*_absorption.pkl"
    pkl = glob.glob(filepath)
    if (len(pkl) == 0):
        print("No absorption file found for given identifier. Stopping.")
        quit()
    elif (len(pkl) > 1):
        print("Multiple absorption files are sharing same identifier. Stopping.")
        print(pkl)
        quit()
    with open(pkl[0], 'rb') as f:
        xsection_function = pickle.load(f)
    return xsection_function(wavelength*1e9)

def get_absorption_spectrum(phoretype):
    """
    Loads absorption spectrum for chosen fluorophore type

    Parameters
    ----------
    phoretype : int
        fluorophore identifier
    """
    identifier = str(phoretype).zfill(3)
    filepath = "dye_spectra/" + identifier + "_*_absorption.pkl"
    pkl = glob.glob(filepath)
    if (len(pkl) == 0):
        print("No absorption file found for given identifier. Stopping.")
        quit()
    elif (len(pkl) != 1):
        print("Multiple absorption files are sharing same identifier. Stopping.")
        print(pkl)
        quit()
    with open(pkl[0], 'rb') as f:
        absorption_function = pickle.load(f)
    return absorption_function

def get_emission_spectrum(phoretype):
    """
    Loads emission spectrum for chosen fluorophore type

    Parameters
    ----------
    phoretype : int
        fluorophore identifier
    """
    identifier = str(phoretype).zfill(3)
    filepath = "dye_spectra/" + identifier + "_*_emission.pkl"
    pkl = glob.glob(filepath)
    if (len(pkl) == 0):
        print("No emission file found for given identifier. Stopping.")
        quit()
    elif (len(pkl) != 1):
        print("Multiple emission files are sharing same identifier. Stopping.")
        print(pkl)
        quit()
    with open(pkl[0], 'rb') as f:
        emission_function = pickle.load(f)
    return emission_function

def get_filter_spectrum(filter_name):
    """
    Loads spectrum of chosen filter

    Parameters
    ----------
    filter_name : str
        filter identifier
    """
    with open("filter_spectra/" + filter_name + ".pkl", 'rb') as f:
        filter_function = pickle.load(f)
    return filter_function

def arbitrary_bandpass_filter_spectrum(low_edge, low_edge_width, high_edge, high_edge_width, transmittivity):
    """
    Creates user-defined bandpass spectrum with linear transition at edges

    Parameters
    ----------
    low_edge : float
        wavelength at which transmittivity increases to half of maximum value [nm]
    low_edge_width : float
        width of transmittivity increase [nm]
    high_edge : float
        wavelength at which transmittivity lowers to half of maximum value [nm]
    high_edge_width : float
        width of transmittivity decrease [nm]
    transmittivity : float
        maximum transmittivity of filter
    """
    if (any((low_edge,low_edge_width,high_edge,high_edge_width,transmittivity)) < 0.0 or 0.5*(low_edge_width+high_edge_width)>(high_edge-low_edge) or low_edge > high_edge):
        print("Error in specified parameters - filter spectrum cannot be computed.\nExiting.")
        quit()
        return None
    xdata = [low_edge-0.5*low_edge_width,low_edge+0.5*low_edge_width,high_edge-0.5*high_edge_width,high_edge+0.5*high_edge_width]
    ydata = [0.0,transmittivity,transmittivity,0.0]
    spectrum = scipy.interpolate.interp1d(xdata,ydata,kind="linear")
    return spectrum

def read_properties(phoretype):
    """
    Reads fluorophore properties from main file.
    [name, absorption cross section, corresponding wavelength, quantum yield, fluorescence lifetime]

    Parameters
    ----------
    phoretype : int
        fluorophore identifier
    """
    out = numpy.genfromtxt("dye_spectra/properties.dat", usecols=(0,3,4,5,6),comments="#",skip_header=phoretype,max_rows=1)
    if (int(out[0]) != phoretype):
        print("File properties.dat either missing lines or containing bad data!")
        quit()
    name = str(numpy.genfromtxt("dye_spectra/properties.dat", dtype="U",usecols=(1),comments="#",skip_header=phoretype,max_rows=1))
    out[4] *= 1e-9  #to account for lifetime being given in nanoseconds
    return name,out[1:5]

def read_STED_properties(phoretype):
    """
    Reads fluorophore properties relevant for stimulated emission depletion from file.
    [vibrational relaxation rate, intersystem crossing yield, triplet state lifetime]

    Parameters
    ----------
    phoretype : int
        fluorophore identifier
    """
    out = numpy.genfromtxt("dye_spectra/STED_properties.dat", usecols=(0,1,2,3),comments="#",skip_header=phoretype,max_rows=1)
    out[3] *= 1e-9  #to account for lifetime being given in nanoseconds
    if (int(out[0]) != phoretype):
        print("File STED_properties.dat either missing lines or containing bad data!")
        quit()
    return out[1:4]

def read_pdf_fit(phoretype):
    """
    Reads parameters of probability distribution that most efficiently samples emission spectrum from file.

    Parameters
    ----------
    phoretype : int
        fluorophore identifier
    """
    out = numpy.genfromtxt("dye_spectra/Laplace_PDFs.dat", usecols=(0,1,2,3,4,5),comments="#",skip_header=phoretype,max_rows=1)
    if (int(out[0]) != phoretype):
        print("File Laplace_PDFs.dat either missing lines or containing bad data!")
        quit()
    return out[2:6]
        
def naive_rejection_sampler(n,spectrum,lowbound,highbound,pdf_parameters):
    """
    Keeps rejection sampling a distribution until it succeeds

    Parameters
    ----------
    n : int
        number of samples to return
    spectrum : obj
        interpolation result
    lowbound : float
        lower bound of spectrum
    highbound : float
        high bound of spectrum
    pdf_parameters : 1D array
        parameters of distribution for sampling from emission spectrum
    """
    M=pdf_parameters[3]
    wavelengths = numpy.empty(n,dtype=float)
    for i in range(n):
        while(True):
            r = scipy.stats.laplace_asymmetric.rvs(loc=pdf_parameters[0], scale=pdf_parameters[1],kappa=pdf_parameters[2])
            if (lowbound <= r <= highbound):
                envelope = M*scipy.stats.laplace_asymmetric.pdf(r,loc=pdf_parameters[0], scale=pdf_parameters[1],kappa=pdf_parameters[2])
                p = numpy.random.uniform(0, envelope)
                if (p < spectrum(r)):
                    wavelengths[i] = r
                    break
    return wavelengths

def naive_rejection_sampler_optimized(n,spectrum,lowbound,highbound,pdf_parameters):
    """
    Keeps rejection sampling a distribution until it succeeds. Optimized version of the above function.
    TODO: Rewrite the second part with a recursive loop so process doesn't fail even with pathological distributions

    Parameters
    ----------
    n : int
        number of samples to return
    spectrum : obj
        interpolation result
    lowbound : float
        lower bound of spectrum
    highbound : float
        high bound of spectrum
    pdf_parameters : 1D array
        parameters of distribution for sampling from emission spectrum
    """
    M=pdf_parameters[3]
    wavelengths = numpy.empty(n,dtype=float)
    rs = scipy.stats.laplace_asymmetric.rvs(loc=pdf_parameters[0], scale=pdf_parameters[1],kappa=pdf_parameters[2],size=(int)(1.1*n))   #generate a bit extra as it is computationally cheap
    envelopes = M*scipy.stats.laplace_asymmetric.pdf(rs,loc=pdf_parameters[0], scale=pdf_parameters[1],kappa=pdf_parameters[2])
    ps = numpy.random.uniform(numpy.zeros_like(envelopes), envelopes)
    spectralvalues = spectrum(rs)
    indices = numpy.where((lowbound <= rs) & (rs <= highbound) & (ps < spectralvalues))[0]
    if (len(indices) >= n):
        return rs[indices][0:n]
    else:
        wavelengths[0:len(indices)] = rs[indices]
        diff = n-len(indices)
        rs2 = scipy.stats.laplace_asymmetric.rvs(loc=pdf_parameters[0], scale=pdf_parameters[1],kappa=pdf_parameters[2],size=(int)(diff*2*n/(len(indices)+1)))
        envelopes2 = M*scipy.stats.laplace_asymmetric.pdf(rs2,loc=pdf_parameters[0], scale=pdf_parameters[1],kappa=pdf_parameters[2])
        ps2 = numpy.random.uniform(numpy.zeros_like(envelopes2), envelopes2)
        spectralvalues2 = spectrum(rs2)
        indices2 = numpy.where((lowbound <= rs2) & (rs2 <= highbound) & (ps2 < spectralvalues2))[0]
        wavelengths[len(indices):] = rs2[indices2][0:diff]
        return wavelengths


def collected_photons_per_exposure(emission_spectrum, filter_spectrum, incident_photons, quantum_yield, detector_qeff, illumination, pdf_parameters, rng):
    """
    Calculates mean number of photons collected by detector during the exposure time

    Parameters
    ----------
    emission_spectrum : obj
        interpolation result
    filter_spectrum : obj
        interpolation result
    incident_photons : 1D array
        mean numbers of photons incident on fluorophores per exposure time
    quantum_eff : float
        quantum efficiency of fluorophore
    detector_quantum_eff : float
        quantum efficiency of detector
    illumination : float
        probability that an emitted photon is collected by microscope optics
    pdf_parameters : 1D array
        parameters of distribution for sampling from emission spectrum
    rng : obj
        random number generator (for checking whether photons get transmitted)
    """
    emitted_photons_at_filter = incident_photons*quantum_yield*illumination
    discretized_counts = emitted_photons_at_filter.astype(int)  # unfortunately, some discretization needs to be done here
    collected_photons = numpy.zeros_like(discretized_counts)
    
    all_photons_to_generate = numpy.sum(discretized_counts)

    low_lambda = filter_spectrum.x[0]
    high_lambda = filter_spectrum.x[-1]

    for photontgt in range(len(incident_photons)):
        photon_wavelengths = naive_rejection_sampler_optimized(discretized_counts[photontgt],emission_spectrum,emission_spectrum.x[0],emission_spectrum.x[-1],pdf_parameters)
        randnums = rng.uniform(size=discretized_counts[photontgt])
        collected_indices = numpy.where((low_lambda < photon_wavelengths[:]) & (photon_wavelengths[:] < high_lambda) & (randnums[:] < filter_spectrum(photon_wavelengths[:])))[0]
        collected_photons[photontgt] = len(collected_indices)         
            
    return collected_photons*detector_qeff

def collected_photons_per_exposure_optimized(emission_spectrum, filter_spectrum, incident_photons, quantum_yield, detector_qeff, illumination, pdf_parameters, rng):
    """
    Calculates mean number of photons collected by detector during the exposure time. Optimized version of the above function.
    TODO: separate detector characteristic into a function that accounts for wavelenth-dependent Qeff

    Parameters
    ----------
    emission_spectrum : obj
        interpolation result
    filter_spectrum : obj
        interpolation result
    incident_photons : 1D array
        mean numbers of photons incident on fluorophores per exposure time
    quantum_eff : float
        quantum efficiency of fluorophore
    detector_quantum_eff : float
        quantum efficiency of detector
    illumination : float
        probability that an emitted photon is collected by microscope optics
    pdf_parameters : 1D array
        parameters of distribution for sampling from emission spectrum
    rng : obj
        random number generator (for checking whether photons get transmitted)
    """
    emitted_photons_at_filter = incident_photons*quantum_yield*illumination
    discretized_counts = emitted_photons_at_filter.astype(int)  # unfortunately, some discretization needs to be done here
    collected_photons = numpy.zeros_like(discretized_counts)
    
    all_photons_to_generate = numpy.sum(discretized_counts)

    low_lambda = filter_spectrum.x[0]
    high_lambda = filter_spectrum.x[-1]

    photon_wavelengths = naive_rejection_sampler_optimized(all_photons_to_generate,emission_spectrum,emission_spectrum.x[0],emission_spectrum.x[-1],pdf_parameters)
    randnums = rng.uniform(size=all_photons_to_generate)

    photonsum = 0
    for photontgt in range(len(incident_photons)):
        collected_indices = numpy.where((low_lambda < photon_wavelengths[photonsum:photonsum+discretized_counts[photontgt]]) & (photon_wavelengths[photonsum:photonsum+discretized_counts[photontgt]] < high_lambda) & (randnums[photonsum:photonsum+discretized_counts[photontgt]] < filter_spectrum(photon_wavelengths[photonsum:photonsum+discretized_counts[photontgt]])))[0]
        collected_photons[photontgt] = len(collected_indices)         
            
    return collected_photons*detector_qeff

def calculate_single_image(phores, incident_photons, filter_spectrum, NA, n, detector_qeff, rng_seed):
    """
    Calculates projected photon counts at detector for all fluorophores. 

    Parameters
    ----------
    phores : 2D array
        types and positions of fluorophores
    intensities : 1D array
        illumination intensities at fluorophores
    filter_spectrum : obj
        interpolation result
    NA : float
        numerical aperture
    n : float
        refractive index of propagation medium
    wavelength : float
        laser wavelength [nm]
    exptime : float
        exposure time [s]
    detector_qeff : float
        quantum efficiency of camera/detector
    rng_seed : int
        seed for random number generator
    """

    phorenum = numpy.shape(phores)[0]
    photon_counts = numpy.empty(phorenum,dtype=float)
    illumination = illumination_fraction(NA,n)

    setseed = numpy.random.SeedSequence(rng_seed)
    rng = numpy.random.Generator(numpy.random.PCG64(setseed))
        
    phoretypes = numpy.unique(phores[:,0]).astype(int)

    for phoretype in phoretypes:
        emission_spectrum = get_emission_spectrum(phoretype)
        pdf_parameters = read_pdf_fit(phoretype)
        _,props = read_properties(phoretype)
        quantum_yield = props[2]

        where_phores = numpy.where(phores[:,0] == phoretype)[0]
        photon_counts[where_phores] = collected_photons_per_exposure_optimized(emission_spectrum, filter_spectrum, incident_photons[where_phores], quantum_yield, detector_qeff, illumination, pdf_parameters, rng)
    
    return photon_counts

def pixel_binning(phores,photon_counts,field_size,pixel_size):
    """
    Bins photon counts at detector according to effective pixel size

    Parameters
    ----------
    phores : 2D array
        types and positions of fluorophores
    photon_counts : 1D array
        numbers of detected photons per fluorophore 
    field_size : float
        dimension of fluorophore field [m]
    pixel_size : float
        length corresponding to microscope resolution [m]
    """
    numbins = math.floor(2*(field_size-0.5*pixel_size)/pixel_size)+1
    range = [[-(numbins+0.5)*pixel_size, (numbins+0.5)*pixel_size], [-(numbins+0.5)*pixel_size, (numbins+0.5)*pixel_size]]
    hist,xedges,yedges = numpy.histogram2d(phores[:,1],phores[:,2],2*numbins+1,range,density=False,weights=photon_counts)
    
    return hist,xedges,yedges

def radial_signal_profile(phores,photon_counts):
    """
    Returns the detected radial profile of the illuminated area

    Parameters
    ----------
    phores : 2D array
        types and positions of fluorophores
    photon_counts : 1D array
        numbers of detected photons per fluorophore 
    """
    phorenum = numpy.shape(phores)[0]
    profile = numpy.empty((phorenum,2),dtype=float)
    profile[:,0] = numpy.sqrt(numpy.square(phores[:,1]) + numpy.square(phores[:,2]))
    profile[:,1] = photon_counts[:]

    maxvalue = numpy.amax(profile[:,1])
    if (maxvalue > 0.0):
        profile[:,1] /= maxvalue
    else:
        return [0,0]
    profile = profile[profile[:, 0].argsort()]
    return profile

def gaussfunc(r,A,w):
    return A*numpy.exp(-2*(r**2)/(w**2))

def linfunc(r,a,b):
    return a*r+b

def FWHM_calculator_Gaussian(profile,w0):
    """
    Calculates FWHM by fitting intensity profile to a Gaussian function

    Parameters
    ----------
    profile : 2D array
        radii and normalized intensities
    w0 : float
        beam diameter at waist [m]
    """
    popt,_ = scipy.optimize.curve_fit(gaussfunc,profile[:,0],profile[:,1],[1.0,w0])

    FWHM = 2*popt[1]*math.sqrt(math.log(2)/2)
    return FWHM,popt

def FWHM_calculator_lin(profile,w0):
    """
    Naively calculates FWHM by fitting intensity profile to a linear function in case the analytic form isn't known 

    Parameters
    ----------
    profile : 2D array
        radii and normalized intensities
    w0 : float
        beam diameter at waist [m]
    """
    profile = profile[profile[:, 1] >= 0.4, :]
    profile = profile[profile[:, 1] <= 0.6, :]

    try:
        popt,_ = scipy.optimize.curve_fit(linfunc,profile[:,0],profile[:,1],[-1/w0,1.0])
        FWHM = 2*(0.5-popt[1])/popt[0]
        return FWHM,popt
    except:
        print("FWHM could not be estimated by linear fitting.")

def pdf_objective_function(params,spectrum,Npts):
    """
    Returns the amount of mismatch between the emission spectrum and asymmetric Laplace partial density function for sampling 

    Parameters
    ----------
    params : 1D array
        array of parameter values (loc, scale, kappa, M)
    spectrum : obj
       interpolation result
    Npts : int
        number of interpolation points per nm
    """
    loc = params[0]
    scale = params[1]
    kappa = params[2]
    M = params[3]

    lowbound = spectrum.x[0]
    highbound = spectrum.x[-1]
    Ntotal = int(Npts*(highbound-lowbound))
    laplace = scipy.stats.laplace_asymmetric(loc=loc,scale=scale,kappa=kappa)

    xlocs = numpy.linspace(lowbound,highbound,num=Ntotal,endpoint=True)
    spectrum_values = numpy.empty((Ntotal))
    laplace_values = numpy.empty((Ntotal))

    f = 0.0
    spectrum_values = spectrum(xlocs)
    laplace_values = M*laplace.pdf(xlocs)
    indices = numpy.where(laplace_values > spectrum_values*1.001)[0]
    f += numpy.sum(numpy.square(laplace_values[indices]-spectrum_values[indices]))
    f += 1e3*(Ntotal-len(indices))

    print("| ", end="",flush=True)
    return f/Ntotal

def optimize_distribution(phoretype,Npts=10,maxiter=100):
    """
    Calculates the optimal distribution function for rejection sampling of emitted photon wavelengths,
    where the optimum distribution defined as the one with fewest wasted photon evaluations. 

    Parameters
    ----------
    phoretype : int
        fluorophore identifier
    Npts : int
        number of interpolation points per nm
    maxiter : int
        maximum number of allowed iterations
    """
    spectrum = get_emission_spectrum(phoretype)
    maxy = 0.0
    locguess = 0.0
    for t in range(len(spectrum.x)):
        if (spectrum.y[t] > maxy):
            maxy = spectrum.y[t]
            locguess = spectrum.x[t]

    guess = [locguess,25.0,0.6,80.0]
    #bounds = [(spectrum.x[0],spectrum.x[-1]),(0.0,None),(None,None),(0.0,None)]
    print("Starting minimization. This may take several minutes.")
    try:
        result = scipy.optimize.minimize(pdf_objective_function,guess,args=(spectrum,Npts),tol=1.0,method="Nelder-Mead",options={"maxiter": maxiter,'disp': True})
        return result.x
    except:
        print("Minimization failed. Please adjust starting parameters.")