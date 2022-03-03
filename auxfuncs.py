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


def gaussian_point_intensity(point_coords, n, wavelength, w0, I0):
    """
    Returns intensity of the TEM00 (Gaussian) beam at a point in space.

    Parameters
    ----------
    point_coords : list of float
        x,y, and z position relative to center of beam waist [m]
    n : float
        refractive index of propagation medium
    wavelength : float
        laser wavelength [m]
    w0 : float
        beam diameter at waist [m]
    I0 : float
        intensity in center of waist [W/m^2]
    """
    if (numpy.isclose(point_coords[2], 0.0, rtol=1e-05, atol=1e-12, equal_nan=False) == True):
        return I0*math.exp(-2*(point_coords[0]**2 + point_coords[1]**2)/(w0*w0))
    else:
        w = w0 * math.sqrt(1 + (point_coords[2]*wavelength/(math.pi*n*w0*w0))**2)
        return I0*((w0/w)**2)*math.exp(-2*(point_coords[0]**2 + point_coords[1]**2)/(w*w))

def STED_2D_approx_point_intensity(point_coords, STEDwavelength, P, NA):
    """
    Intensity of depletion beam using the standing wave approximation, with an analytical solution
    NOTE: does not take into account out-of-focus beam spread!

    Parameters
    ----------
    point_coords : float
        refractive index of propagation medium
    STEDwavelength : float
        wavelength of depletion beam [m]
    P : float
        total beam power [W]
    NA : float
        numerical aperture of depletion beam
    """
    maxradius = STEDwavelength/NA
    radius = math.sqrt(point_coords[0]**2 + point_coords[1]**2)
    if (radius < maxradius):
        square = (math.sin(math.pi*NA*radius/STEDwavelength))**2
        return (2*P/math.pi)*((NA/STEDwavelength)**2) * square
    else:
        return 0.0

def STED_2D_point_intensity(point_coords, STEDwavelength, NA_eff):
    """
    Determines the intensity of depletion beam at chosen point with prediction for a coherent plane wave passing through a 2-pi vortex plate
    NOTE: does not take into account out-of-focus beam spread!

    Parameters
    ----------
    point_coords : float
        refractive index of propagation medium
    STEDwavelength : float
        wavelength of depletion beam [m]
    NA_eff : float
        effective numerical aperture of depletion beam
    """
    k = 2*math.pi/STEDwavelength
    x = k*NA_eff*math.sqrt(point_coords[0]**2+point_coords[1]**2)

    C = (0.5*math.pi*k*NA_eff)**2
    return (C/(x**2))*(scipy.special.jv(1,x)*scipy.special.struve(0,x)-scipy.special.jv(0,x)*scipy.special.struve(1,x))**2

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

    for i in range(samples):
        ilist[i] = STED_2D_point_intensity([rlist[i],0.0,0.0],STEDwavelength,NA_eff)*rlist[i]*2*math.pi

    integral = math.fabs(scipy.integrate.simpson(rlist,ilist))
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
    
def field_add_illumination_intensities(phores, n, wavelength, w0, I0):
    """
    Calculates (expected) illumination intensities for all points

    Parameters
    ----------
    phores : 2D array
        types and positions of fluorophores
    n : float
        refractive index of propagation medium
    wavelength : float
        illumination wavelength [m]
    w0 : float
        beam diameter at waist [m]
    I0 : float
        intensity in center of waist [W/m^2]
    """
    phorenum = numpy.shape(phores)[0]
    intensities = numpy.empty(phorenum,dtype=float)
    for i in range(phorenum):
        intensities[i] = gaussian_point_intensity((phores[i][1],phores[i][2],phores[i][3]), n, wavelength, w0, I0)
    return intensities

def field_STED_illumination_intensities(phores, STEDwavelength, P, NA):
    """
    Calculates (expected) intensities of STED beam at all points

    Parameters
    ----------
    phores : 2D array
        types and positions of fluorophores
    STEDwavelength : float
        wavelength of depletion beam [m]
    P : float
        total beam power [W]
    NA : float
        numerical aperture of depletion beam
    """
    phorenum = numpy.shape(phores)[0]
    STEDintensities = numpy.empty(phorenum,dtype=float)
    for i in range(phorenum):
        STEDintensities[i] = STED_2D_approx_point_intensity((phores[i][1],phores[i][2],phores[i][3]), STEDwavelength, P, NA)
    return STEDintensities

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
        for i in range(phorenum):
            if (phores[i,0] == phoretype):
                incident_counts[i] = incident_photons_per_exposure(exptime,wavelength,xsection,intensities[i])
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
        for i in range(phorenum):
            if (phores[i,0] == phoretype):
                _,prob = STED_CW_rates(STEDintensities[i], Isat, exc_rates[i], ks1, vibrelaxrate, intersystem=intersystem, kt1=kt1)
                incident_counts[i] = prob*incident_photons_per_exposure(exptime,wavelength,xsection,intensities[i])
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
    f = scipy.interpolate.interp1d(data[:,0],data[:,1]*multiplicative_factor)
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
        f = scipy.interpolate.interp1d(data[:,0],data[:,1])

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
    if (len(pkl) != 1):
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
    if (len(pkl) != 1):
        print("Multiple emission files are sharing same identifier. Stopping.")
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
    if (len(pkl) != 1):
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

def naive_rejection_sampler(spectrum,lowbound,highbound,pdf_parameters):
    """
    Keeps rejection sampling a distribution until it succeeds

    Parameters
    ----------
    spectrum : obj
        interpolation result
    lowbound : float
        lower bound of spectrum
    highbound : float
        high bound of spectrum
    pdf_parameters : 1D array
        parameters of distribution for sampling from emission spectrum
    """
    laplace = scipy.stats.laplace_asymmetric(loc=pdf_parameters[0], scale=pdf_parameters[1],kappa=pdf_parameters[2])
    M=pdf_parameters[3]
    while(True):
        r = scipy.stats.laplace_asymmetric.rvs(loc=pdf_parameters[0], scale=pdf_parameters[1],kappa=pdf_parameters[2])
        if (r<lowbound or r>highbound):
            continue
        envelope = M*scipy.stats.laplace_asymmetric.pdf(r,loc=pdf_parameters[0], scale=pdf_parameters[1],kappa=pdf_parameters[2])
        p = numpy.random.uniform(0, envelope)
        if (p < spectrum(r)):
            return r

def collected_photons_per_exposure(emission_spectrum, filter_spectrum, incident_photons, quantum_yield, detector_qeff, illumination, pdf_parameters, rng):
    """
    Calculates mean number of photons collected by detector during the exposure time

    Parameters
    ----------
    emission_spectrum : obj
        interpolation result
    filter_spectrum : obj
        interpolation result
    incident_photons : float
        mean number of photons incident on fluorophore per exposure time
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
    collected_photons = 0
    for i in range(int(emitted_photons_at_filter)):
        photon_wavelength = naive_rejection_sampler(emission_spectrum,emission_spectrum.x[0],emission_spectrum.x[-1],pdf_parameters)
        if (photon_wavelength > filter_spectrum.x[0] and photon_wavelength < filter_spectrum.x[-1]):
            collected_photons += 1
            probability = filter_spectrum(photon_wavelength)
            if (rng.uniform() < probability == True):
                collected_photons += 1
            else:
                continue             
        else:
            continue
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

    for t,phoretype in enumerate(phoretypes):
        emission_spectrum = get_emission_spectrum(phoretype)
        pdf_parameters = read_pdf_fit(phoretype)
        _,props = read_properties(phoretype)
        quantum_yield = props[2]
        for i in range(phorenum):
            if (phores[i,0] == phoretype):
                photon_counts[i] = collected_photons_per_exposure(emission_spectrum, filter_spectrum, incident_photons[i], quantum_yield, detector_qeff, illumination, pdf_parameters, rng)
    
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
    for i in range(phorenum):
        r = math.sqrt(phores[i][1]**2 + phores[i][2]**2)
        profile[i][0] = r
        profile[i][1] = photon_counts[i]

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

    popt,_ = scipy.optimize.curve_fit(linfunc,profile[:,0],profile[:,1],[-1/w0,1.0])
    FWHM = 2*(0.5-popt[1])/popt[0]
    return FWHM,popt

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
    for i in range(Ntotal):
        spectrum_values[i] = spectrum(xlocs[i])
        laplace_values[i] = M*laplace.pdf(xlocs[i])

        if (laplace_values[i] > spectrum_values[i]*1.001):
            f += (laplace_values[i]-spectrum_values[i])**2
        else:
            f += 1e3
    print("| ", end="",flush=True)
    return f/Ntotal

def optimize_distribution(phoretype,Npts=10,maxiter=100):
    """
    Calculates the optimal (least wasted evaluations) distribution function for rejection sampling of emitted photon wavelengths

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
    result = scipy.optimize.minimize(pdf_objective_function,guess,args=(spectrum,Npts),tol=1.0,method="Nelder-Mead",options={"maxiter": maxiter,'disp': True})

    return result.x