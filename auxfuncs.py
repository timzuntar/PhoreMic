import math
import glob
import pickle
import numpy
from numba import jit
import scipy
import scipy.stats
import scipy.interpolate
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
    square = (math.sin(math.pi*NA*math.sqrt(point_coords[0]**2 + point_coords[1]**2)/STEDwavelength))**2
    return 4*P*(NA/STEDwavelength)**2 * square

def STED_2D_point_intensity(point_coords, STEDwavelength, R, f):
    """
    Determines the intensity of depletion beam from total incident power (assuming losses are constant) via numerical integration

    Parameters
    ----------
    point_coords : float
        refractive index of propagation medium
    STEDwavelength : float
        wavelength of depletion beam [m]
    R : float
        radius of lens or filter aperture [m]
    f : float
        lens focal length [m]
    """
    k = 2*math.pi/STEDwavelength
    x = k*R*math.sqrt(point_coords[0]**2+point_coords[1]**2)/f
    C = (0.5*math.pi*k*R*R/f)**2
    return (C/(x**2))*(scipy.special.jv(1,x)*scipy.special.struve(0,x)-scipy.special.jv(0,x)*scipy.special.struve(1,x))**2

def STED_2D_intensity_normalization(n, STEDwavelength, P, R, f, lossfactor=0.0, integrationwidth=5):
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
    R : float
        radius of lens or filter aperture [m]
    f : float
        lens focal length [m]
    lossfactor : float
        relative amount of beam power not incident on the sample
    integrationwidth : float
        maximum integration radius in multiples of peak intensity radius
    """
    peakradius = f*2.45*STEDwavelength/(2*math.pi*R)
    samples = 2001
    rlist = numpy.linspace(0.0, peakradius*integrationwidth, num=samples, endpoint=True)
    ilist = numpy.empty(samples)

    for i in range(samples):
        ilist[i] = STED_2D_point_intensity((rlist[i],0.0,0.0),STEDwavelength,R,f)
    integral = scipy.integrate.simpson(ilist,rlist)
    return P/integral

def generate_fluorophore_field(w0, density, phoretype, existing=None, seed=42, sizemultiplier=5, maxnum=1e7):
    """
    Randomly generates positions of fluorescing molecules

    Parameters
    ----------
    w0 : float
        beam waist diameter [m]
    density : float
        number density of fluorescing molecules [per square meter]
    phoretype : int
        fluorophore identifier
    existing : 2D array
        previous array; output of function consists of new entries concatenated with previous ones
    seed : int
        seed to random number generator
    sizemultiplier : float
        width and height of created field in multiples of beam waist diameter
    maxnum : int
        maximum number of generated points
    """
    phorenum = int((w0*2*sizemultiplier)**2 * density)
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
    xpositions = rng.uniform(-w0*sizemultiplier,w0*sizemultiplier,phorenum)
    ypositions = rng.uniform(-w0*sizemultiplier,w0*sizemultiplier,phorenum)
    zpositions = numpy.full(phorenum, 0.0, dtype=float)
    types = numpy.full(phorenum, phoretype, dtype=int)
    phores = numpy.concatenate((types, xpositions, ypositions,zpositions)).reshape((-1, 4), order='F')
    if (existing==None):
        return phores
    else:
        return numpy.vstack((existing,phores))
    
    
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

def illumination_fraction(NA, n):
    """
    Fraction of solid angle over which emitted light is collected

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
        absorption cross section for given wavelength [square meters]
    intensity : float
        local illumination intensity
    """
    return exptime*xsection*intensity*wavelength/(scipy.constants.c*scipy.constants.h)

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
    return xsection_function(wavelength)

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
    with open("filter_spectra/" + filter_name + ".pkl", 'rb') as f:
        filter_function = pickle.load(f)
    return filter_function

def read_properties(phoretype):
    """
    Reads fluorophore properties from main file.
    [absorption cross section, corresponding wavelength, quantum yield, fluorescence lifetime]

    Parameters
    ----------
    phoretype : int
        fluorophore identifier
    """
    out = numpy.genfromtxt("dye_spectra/properties.dat", usecols=(0,2,3,4,5),comments="#",skip_header=phoretype,max_rows=1)
    if (int(out[0]) != phoretype):
        print("File properties.dat either missing lines or containing bad data!")
        quit()
    return out[1:5]


def naive_rejection_sampler(spectrum,lowbound,highbound):
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
    """
    #NOTE: The following sampling function is a placeholder for one specific fluorophore type
    #and will be replaced by an existing implementation.
    #Envelope function was determined by trial and error.
    laplace = scipy.stats.laplace_asymmetric(loc=512.0, scale=25,kappa=0.6)
    M=80.0
    while(True):
        r = scipy.stats.laplace_asymmetric.rvs(loc=512.0, scale=25,kappa=0.6)
        if (r<lowbound or r>highbound):
            continue
        envelope = M*scipy.stats.laplace_asymmetric.pdf(r,loc=512.0, scale=25,kappa=0.6)
        p = numpy.random.uniform(0, envelope)
        if (p < spectrum(r)):
            return r

def collected_photons_per_exposure(emission_spectrum, filter_spectrum, incident_photons, quantum_yield, detector_qeff, illumination, rng):
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
    rng : obj
        random number generator (for checking whether photons get transmitted)
    """
    emitted_photons_at_filter = incident_photons*quantum_yield*illumination
    collected_photons = 0
    for i in range(int(emitted_photons_at_filter)):
        photon_wavelength = naive_rejection_sampler(emission_spectrum,emission_spectrum.x[0],emission_spectrum.x[-1])
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

def calculate_single_image(phores, intensities, filter_spectrum, NA, n, wavelength, exptime, detector_qeff, rng_seed):
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
        xsection = get_absorption_xsection(phoretype,wavelength)
        emission_spectrum = get_emission_spectrum(phoretype)
        _,_,quantum_yield,_ = read_properties(phoretype)
        for i in range(phorenum):
            if (phores[i,0] == phoretype):
                incident_photons = incident_photons_per_exposure(exptime, wavelength, xsection, intensities[i])
                photon_counts[i] = collected_photons_per_exposure(emission_spectrum, filter_spectrum, incident_photons, quantum_yield, detector_qeff, illumination, rng)
    
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
    numbins = math.floor(0.5*(field_size-pixel_size)/pixel_size)
    range = [[-(numbins+0.5)*pixel_size, (numbins+0.5)*pixel_size], [-(numbins+0.5)*pixel_size, (numbins+0.5)*pixel_size]]
    hist,xedges,yedges = numpy.histogram2d(phores[:,1],phores[:,2],2*numbins+1,range,density=False,weights=photon_counts)
    
    return hist,xedges,yedges