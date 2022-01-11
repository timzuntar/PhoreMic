import math
import numpy
import scipy

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
        laser wavelength [m]
    xsection : float
        absorption cross section for given wavelength [square meters]
    intensity : float
        local illumination intensity
    """
    return exptime*xsection*intensity*wavelength/(scipy.constants.c*scipy.constants.h)