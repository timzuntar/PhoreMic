import auxfuncs as aux
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy

def display_2D_fluorophore_field(phores, w0, latmultiplier,Pexc,wavelength):
    """
    Displays positions of each fluorophore species and the excitation beam waist

    Parameters
    ----------
    phores : 2D array
        types and positions of fluorophores
    w0 : float
        beam waist diameter [microns]
    latmultiplier : float
        width and height of created field in multiples of beam waist diameter
    Pexc : float
        total power of excitation beam [W]
    wavelength : float
        excitation beam wavelength [m]
    """
    phoretypes = numpy.unique(phores[:,0]).astype(int)
    beamlabel = "excitation beam waist\n$w_0 = %.2f ~\mu m$\n$P = %.0f ~ \mu W$\n$\lambda = %.0f ~nm$" % (w0*1e6,Pexc*1e6,wavelength*1e9)
    beamwaist = plt.Circle((0, 0), w0*1e6, color='r',fill=False,label=beamlabel)

    fig = plt.figure(figsize=(9, 6), constrained_layout=True)
    ax = fig.add_subplot()

    ax.add_patch(beamwaist)
    for type in phoretypes:
        typelabel,_ = aux.read_properties(type)
        ax.scatter(phores[phores[:, 0] == type,1]*1e6,phores[phores[:,0] == type,2]*1e6,label=typelabel,s=3)
    
    ax.axis("scaled")
    plt.legend(bbox_to_anchor=(1.00,0.9), borderaxespad=0)

    plt.title("Generated fluorophore map")
    plt.xlabel(r"X [$\mu$m]")
    plt.xlim([-w0*latmultiplier*1e6*1.1,w0*latmultiplier*1e6*1.1])
    plt.ylabel(r"Y [$\mu$m]")
    plt.ylim([-w0*latmultiplier*1e6*1.1,w0*latmultiplier*1e6*1.1])
    plt.show()
    return None

def display_spectral_overview(phoretype,absorption_spectrum,emission_spectrum,filter_spectrum,exc_wavelength,alt_wavelength=None,alt_type=None):
    """
    Displays the various relevant spectra to aid in filter/illumination choice

    Parameters
    ----------
    phoretype : int
        fluorophore identifier
    absorption_spectrum : obj
        fluorophore absorption spectrum
    emission_spectrum : obj
        fluorophore emission spectrum
    filter_spectrum : obj
    exc_wavelength : float
        wavelength of excitation beam [m]
    alt_wavelength : float
        wavelength of second beam, if such is present [m]
    alt_type : str
        identifier of second method
    """
    phorename,specs = aux.read_properties(phoretype)

    lowbounds = [absorption_spectrum.x[0],emission_spectrum.x[0],filter_spectrum.x[0]]
    highbounds = [absorption_spectrum.x[-1],emission_spectrum.x[-1],filter_spectrum.x[-1]]
    min_wavelength = min(lowbounds)
    max_wavelength = max(highbounds)

    fig = plt.figure(figsize=(8, 4), constrained_layout=True)
    ax = fig.add_subplot()

    absorption_x = numpy.linspace(absorption_spectrum.x[0], absorption_spectrum.x[-1], num=int((absorption_spectrum.x[-1]-absorption_spectrum.x[0])*5+1), endpoint=True)
    emission_x = numpy.linspace(emission_spectrum.x[0], emission_spectrum.x[-1], num=int((emission_spectrum.x[-1]-emission_spectrum.x[0])*5+1), endpoint=True)
    filter_x = numpy.linspace(filter_spectrum.x[0], filter_spectrum.x[-1], num=int((filter_spectrum.x[-1]-filter_spectrum.x[0])*5+1), endpoint=True)

    absorption_vec = numpy.vectorize(absorption_spectrum)
    emission_vec = numpy.vectorize(emission_spectrum)
    filter_vec = numpy.vectorize(filter_spectrum)

    absorption_y = absorption_vec(absorption_x)
    absorption_y /= numpy.max(absorption_y)
    emission_y = emission_vec(emission_x)
    filter_y = filter_vec(filter_x)

    plt.plot(absorption_x,absorption_y,"-b",label="absorption")
    plt.plot(emission_x,emission_y,"-g",label="emission")
    plt.plot(filter_x,filter_y,"-k",label="filter spectrum")

    exc_label = "excitation beam, $\lambda = %.0f$ nm" % (exc_wavelength*1e9)
    plt.axvline(x=exc_wavelength*1e9,color="r",label=exc_label,linestyle="--")
    if (alt_type == "STED"):
        alt_label = "depletion beam, $\lambda = %.0f$ nm" % (alt_wavelength*1e9)
        plt.axvline(x=alt_wavelength*1e9,color="k",label=alt_label,linestyle="--")

    plt.legend(loc=2)
    plt.title("Overview of spectra for " + phorename + " and test filter")
    plt.xlabel(r"$\lambda$ [nm]")
    plt.xlim([min_wavelength,max_wavelength])
    plt.ylabel(r"Normalized profiles")
    plt.ylim([0.0,1.05])
    plt.show()

    return None

def display_detected_photon_counts(phores,w0,photon_counts):
    """
    Displays the detected photon numbers from each fluorophore 

    Parameters
    ----------
    phores : 2D array
        types and positions of fluorophores
    w0 : float
        beam waist diameter [microns]
    photon_counts : 1D array
        numbers of emitted photons collected by detector
    """
    beamwaist = plt.Circle((0, 0), w0*1e6, color='r',fill=False)
    
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.add_patch(beamwaist)
    points = ax.scatter(x=phores[:,1]*1e6,y=phores[:,2]*1e6,c=photon_counts,vmin=numpy.min(photon_counts), vmax=numpy.max(photon_counts),s=1)
    plt.colorbar(points)
    ax.axis('equal')
    
    plt.title("Number of detected photons per fluorophore")
    plt.xlabel(r"X [$\mu$m]")
    plt.ylabel(r"Y [$\mu$m]")
    plt.show()

    return None

def display_photon_counts_side_by_side(phores,Pexc,wavelength,PSTED,STEDwavelength,w0,photon_counts,alt_w0,alt_photon_counts,latmultiplier,alt_type=None):
    """
    Displays detected photon numbers from each fluorophore for two different imaging methods or runs 

    Parameters
    ----------
    phores : 2D array
        types and positions of fluorophores
    w0 : float
        excitation beam waist diameter [microns]
    photon_counts : 1D array
        numbers of emitted photons collected by detector for first method
    alt_w0 : float
        beam waist diameter for second method [microns]
    alt_photon_counts : 1D array
        numbers of emitted photons collected by detector for second method
    alt_type : str
        identifier of second method
    """
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 6),constrained_layout=True)
    
    beamwaist1 = plt.Circle((0, 0), w0*1e6, color='r',fill=False)
    axes[0].add_patch(beamwaist1)
    axes[0].set_title("Gaussian beam only")

    if (alt_type == None):
        beamwaist2 = plt.Circle((0, 0), alt_w0*1e6, color='r',fill=False)
        axes[1].add_patch(beamwaist2)
        axes[1].set_title("Gaussian beam only")
    elif (alt_type == "STED"):
        beamwaist2 = plt.Circle((0, 0), w0*1e6, color='r',fill=False)
        beamwaist3 = plt.Circle((0, 0), alt_w0*1e6, color='k',fill=False)
        axes[1].add_patch(beamwaist2)
        axes[1].add_patch(beamwaist3)
        axes[1].set_title("STED illumination")

    min1 = numpy.min(photon_counts)
    min2 = numpy.min(alt_photon_counts)
    max1 = numpy.max(photon_counts)
    max2 = numpy.max(alt_photon_counts)

    points = axes[0].scatter(x=phores[:,1]*1e6,y=phores[:,2]*1e6,c=photon_counts,vmin=min(min1,min2), vmax=max(max1,max2),s=3)
    axes[1].scatter(x=phores[:,1]*1e6,y=phores[:,2]*1e6,c=alt_photon_counts,vmin=min(min1,min2), vmax=max(max1,max2),s=3)
    axes[0].axis("scaled")
    axes[1].axis("scaled")
    axes[0].set_xlim([-w0*latmultiplier*1e6*1.1,w0*latmultiplier*1e6*1.1])
    axes[0].set_ylim([-w0*latmultiplier*1e6*1.1,w0*latmultiplier*1e6*1.1])
    axes[1].set_xlim([-w0*latmultiplier*1e6*1.1,w0*latmultiplier*1e6*1.1])
    axes[1].set_ylim([-w0*latmultiplier*1e6*1.1,w0*latmultiplier*1e6*1.1])
    axes[0].set_xlabel(r"X [$\mu$m]")
    axes[0].set_ylabel(r"Y [$\mu$m]")
    axes[1].set_xlabel(r"X [$\mu$m]")
    axes[1].set_ylabel(r"Y [$\mu$m]")
    plt.suptitle("Number of detected photons per fluorophore\n$P_{exc} = %.2f ~\mu W, \lambda_{exc} = %.0f ~nm, P_{STED} = %.2f ~W, \lambda_{STED} = %.0f ~nm$" % (Pexc*1e6,wavelength*1e9,PSTED,STEDwavelength*1e9))
    plt.colorbar(points,shrink=0.7)
    plt.show()
    return None

def display_detected_image(hist):
    """
    Displays the detected "image" histograms

    Parameters
    ----------
    hist : 2D array
        histogram of photon numbers
    """
    fig = plt.figure()
    ax = fig.add_subplot()
    plt.imshow(hist)
    plt.title("Number of detected photons per pixel")
    ax.set_aspect('equal')
    plt.colorbar()
    plt.xlabel("X [px]")
    plt.ylabel("Y [px]")
    plt.show()

    return None

def display_detected_images(pixel_size,hist1,hist2,alt_type=None):
    """
    Displays the detected "image" histograms of two different imaging methods or runs 

    Parameters
    ----------
    pixel_size : float
        length corresponding to microscope resolution [m]
    hist1 : 2D array
        histogram of photon numbers for first method
    hist2: 2D array
        histogram of photon numbers for first method
    alt_type : str
        identifier of second method
    """
    numbins = numpy.shape(hist1)
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 6),constrained_layout=True)
    axes[0].set_title("Gaussian beam only")
    if (alt_type == None):
        axes[1].set_title("Gaussian beam only")
    elif (alt_type == "STED"):
        axes[1].set_title("STED illumination")
    
    min1 = numpy.min(hist1)
    min2 = numpy.min(hist2)
    max1 = numpy.max(hist1)
    max2 = numpy.max(hist2)

    points = axes[0].imshow(hist1,vmin=min(min1,min2), vmax=max(max1,max2))
    axes[1].imshow(hist2,vmin=min(min1,min2), vmax=max(max1,max2))

    axes[0].axis("scaled")
    axes[1].axis("scaled")

    axes[0].set_xlabel(r"X [px]")
    axes[0].set_ylabel(r"Y [px]")
    axes[1].set_xlabel(r"X [px]")
    axes[1].set_ylabel(r"Y [px]")

    plt.suptitle("Number of detected photons per pixel\n pixel size = $%.2f ~\mu m$" % (pixel_size*1e6))
    plt.colorbar(points,shrink=0.7)
    plt.show()

    return None


def compare_profiles(regular_profile,STED_profile,size,res1=None,res2=None,popt1=None,popt2=None):
    """
    Compares normalized radial profiles of an illuminated point at each individual molecule
     
    Parameters
    ----------
    regular_profile : 2D array
        profile of regular fluorescence microscopy (2 columns)
    STED_profile : 2D array
        profile of STED-assisted microscopy (2 columns)
    size : float
        radius to which the profiles are evaluated
    res1 : float
        FWHM of response for ordinary excitation
    res2 : float
        FWHM of response with depletion beam active
    popt1 : tuple
        best fit parameters of ordinary excitation
    popt2 : tuple
        best fit parameters with depletion beam active
    """
    fig = plt.figure()
    ax = fig.add_subplot()

    pts = len(regular_profile[:,0])
    if (res1 is None or popt1 is None):
        fitlabel1 = "without STED"
    else:
        fit1 = numpy.empty(pts)
        for i in range(pts):
            fit1[i] = aux.gaussfunc(regular_profile[i,0],popt1[0],popt1[1])
        fitlabel1 = "without STED\nFWHM = $%.3f ~\mu m$" % (res1*1e6)
        ax.plot(regular_profile[:,0]*1e6,fit1,color="r")
        plt.axvline(x=0.5*res1*1e6,color="r",linestyle="--")

    if (res2 is None or popt2 is None):
        fitlabel2 = "with STED"
    else:
        fit2 = numpy.empty(pts)
        for i in range(pts):
            fit2[i] = aux.linfunc(STED_profile[i,0],popt2[0],popt2[1])
        fitlabel2 = "with STED\nFWHM = $%.3f ~\mu m$" % (res2*1e6)
        ax.plot(STED_profile[:,0]*1e6,fit2,color="k")
        plt.axvline(x=0.5*res2*1e6,color="k",linestyle="--")

    ax.scatter(regular_profile[:,0]*1e6,regular_profile[:,1],color="r",label=fitlabel1,s=2)
    ax.scatter(STED_profile[:,0]*1e6,STED_profile[:,1],color="k",label=fitlabel2,s=2)

    plt.legend(bbox_to_anchor=(0.6,0.9), borderaxespad=0)
    plt.title("Per-fluorophore radial profiles")
    plt.xlabel(r"r [$\mu$m]")
    plt.ylabel(r"$I/I_{max}$")
    plt.xlim([0.0, size*1e6])
    plt.ylim([0.0,1.05])
    plt.show()
    return None