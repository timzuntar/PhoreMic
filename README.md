## PhoreMic
### A simple fluorescence microscopy visualization tool

### Purpose

The eventual goal is to achieve (somewhat) accurate simulation and visualization of fluorescence processes and imaging for several superresolution methods, like STED [[1]](#1) [[2]](#2) and two-photon microscopy, in the hopes that it becomes useful for demonstrations and back-of-the-napkin feasibility calculations. Since many important physical processes are abstracted and performance is unoptimized, it neither aims to be an alternative to professionally developed solutions nor claims to accurately predict results achievable with *in vivo* experiments.

### Function

The code aims to eventually support a variety of commonly used dyes, but this is not yet the case. To run simulations, users first need to specify several physical properties of the desired fluorescing species:

- normalized absorption and emission spectrum
- lifetime
- quantum efficiency
- attenuation cross-section at a given wavelength

Defining a new species is straightforward: in the folder *dye_spectra*, you add an unused numerical identifier (000-999) to the *properties.dat* file and specify the relevant fluorophore properties. *STED_properties.dat* is the equivalent file for properties applicable to STED illumination. Spectra are stored as *pkl* objects to simplify handling - functions which interpolate a user-supplied file of intensity as a function of wavelength and save it as such are included. File names must follow a relatively rigid pattern: *yyy_dyename_absorption.pkl* and *yyy_dyename_emission.pkl*, yyy being the identifier. Additional information about vibrational transition lifetimes or bleaching/quenching rates may become relevant in the future, in which case the structure of the files will be updated.

Aside from fluorophore properties, information about the excitation beam and optical setup are needed to calculate the output:

- filter transmission spectrum
- excitation beam wavelength (assumed to be monochromatic)
- beam waist diameter / numerical aperture
- detector efficiency and resolution

Currently, running the code simply results in the simulation of a single image frame either for simple fluorescence microscopy or a side-by-side comparison with STED illumination enabled. By default, all molecules are assumed to be distributed in the focal plane, and molecular density is interpreted as surface density by the function *generate_fluorophore_field*; passing *volume=True* to it generates a 3D distribution, with the density interpreted as ordinary 3D density. The lateral and axial limits of the generated field are defined as the laser beam waist width times *latmultiplier*/*axmultiplier*. The preset value should provide a decent tradeoff between accuracy and speed. 

"Image" in the above sense refers to the detected PSF of a single illuminated point as computed by *calculate_single_image*, the size of which is bounded by the beam focusing and pixel size (aberrations are not taken into account). To simulate imaging of an entire sample area, a nonuniform fluorophore distribution needs to be imaged at a large number of focal points (planned functionality).

### Workflow

When you specify the desired type and density of fluorophores, the code loads their relevant properties and spectra and first randomly generates a planar/volumetric distribution of molecules. For chosen illumination parameters, the mean number of absorbed photons per exposure time is calculated for each molecule from the absorption spectrum and local intensity.

<figure>
  <img
  align="middle"
  src="https://github.com/timzuntar/PhoreMic/blob/main/readme_images/ex_phore_map_N=4500.png?raw=true"
  alt="Generated fluorophores"
  width="600">
  <figcaption></figcaption>
</figure>


From this calculation, the code estimates the number of emission photons collected by the objective (arriving at the filter) from each fluorophore. Their wavelengths are not fixed at generation time, but are determined *post facto* by rejection sampling from the emission spectrum (this is by far the most computationally intensive part of the simulation). Another rejection sampler decides if the photon makes it through filters to the detector. Performance improvements will hopefully come, but to ensure that the model is compatible with the physical reality the sampling cannot be done in a single shot as the generated-but-rejected photons and those that weren't generated in the first place could not be distinguished in this case.

The next step depends on whether a secondary beam is specified. In the case of STED microscopy, first the appropriate beam and saturation intensities are calculated. The mean number of absorbed photons is then recomputed as the expected number to undergo *spontaneous* decay, based on equations from [[3]](#3). The rest of the process remains unchanged. Handling of STED is still provisional- for example, cross-sections for stimulated emission are assumed to be equal to those of spontaneous emission at the same wavelength.

### Results

An overview of the spectra and beam wavelengths used in this test case is visible below and can be generated for quick visualization of e.g. filter suitability. Note that the y-axis label is a slight misnomer; the filter spectrum is not normalized, since the transmittance values do not necessarily reach unity. Besides predefined filter profiles, simple bandpass profiles can be generated on the fly.

<figure>
  <img
  align="middle"
  src="https://github.com/timzuntar/PhoreMic/blob/main/readme_images/test_spectra.png?raw=true"
  alt="Spectral overview"
  width="600">
  <figcaption></figcaption>
</figure>

The below image shows the number of detected photons emitted by each of the 4500 fluorescing molecules with absorption and emission spectra corresponding to Alexa Fluor 488 (but some placeholder values) simultaneously illuminated by a 490 nm 0.01 mW excitation beam (left) with a 1-micron waist and a 592 nm 10 W "donut-shaped" depletion beam (right).

<figure>
  <img
  align="middle"
  src="https://github.com/timzuntar/PhoreMic/blob/main/readme_images/ex_photon_map_N=4500_PSTED=1000Pex.png?raw=true"
  alt="Photon maps"
  width="700">
  <figcaption></figcaption>
</figure>

Finally, an "image" of the single observed point is generated as a histogram of pixel intensity values. Shown here without (left) and with STED illumination (right). These can be quickly recomputed from the detected photon counts if you wish to visualize the impact of different detector pixel sizes.

<figure>
  <img
  align="middle"
  src="https://github.com/timzuntar/PhoreMic/blob/main/readme_images/ex_pixel_map_N=4500_PSTED=1000Pex.png?raw=true"
  alt="Pixel maps at detector"
  width="700">
  <figcaption></figcaption>
</figure>


A comparison between radial intensity profiles of both methods ("intensity" in this case being the number of detected photons from each fluorophore normalized to the maximum) reveals a strong narrowing of the profile when STED illumination is used, as one would expect.

<figure>
  <img
  align="middle"
  src="https://github.com/timzuntar/PhoreMic/blob/main/readme_images/ex_radial_comparison_N=4500_PSTED=3e6Pex.png?raw=true"
  alt="Comparison of radial profiles"
  width="500">
  <figcaption></figcaption>
</figure>

### GUI

A simple interface is eventually planned; it will be implemented either in tkinter or PyQT. 

### Note on beam profiles

The excitation light is currently assumed to be an ideal Gaussian beam propagating along the z-axis, with the central sample plane defined as perpendicular to the beam at its waist (x,y,0). The depletion beam is either assumed to have an idealized sine-squared radial profile or a somewhat more involved approximation of the shape a coherent plane wave takes on when passing through a vortex phase plate. The latter is based on the derivation from [[4]](#4), but is currently WIP. It uses an effective numerical aperture value that can differ from the actual value by a non-insignificant amount, as the intensity is otherwise overestimated by up to a factor of 2 in parts of the beam compared with the rough approximation.

<figure>
  <img
  align="middle"
  src="https://github.com/timzuntar/PhoreMic/blob/main/readme_images/STED_approximation_profiles.png?raw=true"
  alt="Comparison of approximated STED beam profiles"
  width="400">
  <figcaption></figcaption>
</figure>


### Determination of probability distribution for sampling

To minimize the number of wasted dice rolls, the probability distribution used to rejection sample wavelengths of emitted photons needs to follow the emission spectrum as closely as possible. Asymmetric Laplace distributions are a good fit for the typical profile shape. The optimal distribution is computed once with a minimization process (this can take a while). The resulting parameters are then stored in *Laplace_PDFs.dat* for use in all future simulations of that specific fluorophore.

<figure>
  <img
  align="middle"
  src="https://github.com/timzuntar/PhoreMic/blob/main/readme_images/test_emission_pdf.png?raw=true"
  alt="Emission probability distributions"
  width="400">
  <figcaption></figcaption>
</figure>

### References

<a id="1">[1]</a> 
Hell, S. W. and Wichmann, J. (1994). 
Breaking the diffraction resolution limit by stimulated emission: stimulated-emission-depletion fluorescence microscopy.
Optics Letters, 19(11), 780-782.

<a id="2">[2]</a> 
Klar, T. A. and Hell, S. W. (1999). 
Subdiffraction resolution in far-field fluorescence microscopy.
Optics Letters, 24(14), 954-956.

<a id="3">[3]</a> 
Leutenegger, M., Eggeling, C. and Hell, S. W. (2010). 
Analytical description of STED microscopy performance.
Optics Express, 18(25), 26417-26429.

<a id="4">[4]</a> 
Neupane, B., Chen, F., Sun, W., Chiu, D.T. and Wang, G. (2013)
Tuning donut profile for spatial resolution in stimulated emission depletion microscopy.
Rev Sci Instrum. 84(4), 043701

### Requirements

(conservative, you may get away with older versions of all except SciPy)

- Python 3.7.8
- SciPy 1.6.0
- NumPy 1.19.1
- Matplotlib 3.2.2