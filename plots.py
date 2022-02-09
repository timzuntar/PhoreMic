import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy

def display_2D_fluorophore_field(phores, w0):
    """
    Displays positions of each fluorophore species and the excitation beam waist

    Parameters
    ----------
    phores : 2D array
        types and positions of fluorophores
    w0 : float
        beam waist diameter [microns]
    """
    phoretypes = numpy.unique(phores[:,0]).astype(int)

    beamwaist = plt.Circle((0, 0), w0*1e6, color='r',fill=False,label="beam waist")

    fig = plt.figure()
    ax = fig.add_subplot()

    ax.add_patch(beamwaist)
    for type in phoretypes:
        typelabel = "fluorophore type " + str(type)
        ax.scatter(phores[phores[:, 0] == type,1]*1e6,phores[phores[:,0] == type,2]*1e6,label=typelabel,s=1)
    
    ax.axis('equal')
    fig.legend()

    plt.title("Generated fluorophore map")
    plt.xlabel(r"X [$\mu$m]")
    plt.ylabel(r"Y [$\mu$m]")
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

def display_detected_image(hist):
    """
    Displays the detected "image"

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

def compare_profiles(regular_profile,STED_profile,size):
    """
    Compares normalized radial profiles of an illuminated point at each individual molecule
     
    Parameters
    ----------
    regular_profile : 2D array
        profile of regular fluorescence microscopy (2 columns)
    STED_profile : 2D array
        profile of STED-assisted microscopy (2 columns)
    """
    fig = plt.figure()
    ax = fig.add_subplot()

    ax.scatter(regular_profile[:,0]*1e6,regular_profile[:,1],label="without STED",s=2)
    ax.scatter(STED_profile[:,0]*1e6,STED_profile[:,1],label="with STED",s=2)
    fig.legend()

    plt.title("Per-fluorophore radial profiles")
    plt.xlabel(r"r [$\mu$m]")
    plt.ylabel(r"$I/I_{max}$")
    plt.xlim([0.0, size*1e6])
    plt.ylim([0.0,1.05])
    plt.show()
    return None