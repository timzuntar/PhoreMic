## A simple fluorescence microscopy visualization tool

### Purpose

The eventual goal is to achieve (somewhat) accurate simulation and visualization of fluorescence imaging for several methods, like STED [[1]](#1) [[2]](#2) and two-photon microscopy, in the hopes that it becomes useful for demonstrations and back-of-the-napkin feasibility calculations. Since many important physical processes are abstracted and performance is unoptimized, it does not aim to be an alternative to professionally developed solutions.

### Function

Knowledge of the physical properties of the desired fluorescing species (normalised absorption and emission spectrum, lifetime, quantum efficiency, attenuation cross-section at a give wavelength) is a prerequisite for simulation. The code aims to eventually support a variety of commonly used dyes. Adding them yourself is straightforward: in the folder *dye_spectra*, you add an unused numerical identifier (000-999) to the *properties.dat* file and specify the relevant fluorophore properties. Spectra are stored as *pkl* objects to simplify - functions which interpolate a user-supplied intensity(wavelength) file and save it as such are included. File names must follow the pattern *yyy_something_absorption.pkl* and *yyy_something_emission.pkl", yyy being the identifier. Additional information about vibrational transition lifetimes or bleaching/quenching rates may become relevant in the future, in which case property handling will be updated.

Aside from fluorophore properties, information about the excitation beam and optical setup (filter transmission spectrum, numerical aperture, detector efficiency and resolution) are needed to calculate the output. Currently the assumption is that all molecules exist in the focal plane and running the code simply results in the simulation of a single image frame.

### Current progress

When you specify the desired type and density of fluorophores, the code loads the relevant properties and spectra and first randomly generates a planar distribution of molecules. Depending on the illumination (a monochromatic laser beam is assumed), the mean number of absorbed photons per exposure time is calculated for each molecule from the absorption spectrum and local intensity.

<figure>
  <img
  align="middle"
  src="https://github.com/timzuntar/PhoreMic/blob/main/readme_images/ex_phore_map_N=4500.png?raw=true"
  alt="Generated fluorophores"
  width="600">
  <figcaption></figcaption>
</figure>


From that, the code estimates the number of emission photons collected by the objective (arriving at the filter) from each fluorophore. Their wavelengths are then determined *post facto* by rejection sampling from the emission spectrum (this is by far the most computationally intensive part of the simulation). Another rejection sampler decides if the photon makes it through the filter to the detector.

If a depletion beam is specified (for STED microscopy), first the appropriate beam and saturation intensities are calculated. The mean number of absorbed photons is then recomputed as the expected number to undergo *spontaneous* decay, based on equations from [[3]](#3). The rest of the process remains unchanged. Handling of STED is still provisional- for example, cross-sections for stimulated emission and vibrational relaxation rates need to be implemented on a per-species basis.

### Results

An overview of the spectra and beam wavelengths used in this test case should help with figuring out what's going on. Note that the y-axis label is a slight misnomer; the filter spectrum is not normalized, since it consists of transmittance values.

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

Finally, an "image" of the single observed point is generated as a histogram of pixel intensity values. Shown here without (left) and with STED illumination (right).

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
  src="https://github.com/timzuntar/PhoreMic/blob/main/readme_images/ex_radial_comparison_N=4500_PSTED=1000Pex.png?raw=true"
  alt="Comparison of radial profiles"
  width="500">
  <figcaption></figcaption>
</figure>

### Note on beam profiles
The excitation light is currently assumed to be an ideal Gaussian beam. The depletion beam is either assumed to have an idealized sine-squared radial profile or a somewhat more involved approximation of the shape a coherent plane wave takes on when passing through a vortex phase plate. The latter is based on the derivation from [[4]](#4), but is currently WIP. It uses an effective numerical aperture value that can differ from the actual value by a non-insignificant amount, and otherwise overestimates the intensity by up to a factor of 2 in parts of the beam in comparison with the rough approximation.

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
