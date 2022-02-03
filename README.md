## A simple fluorescence microscopy visualization tool
### NOTE: 

### Purpose

The eventual goal is to achieve (somewhat) accurate simulation and visualization of fluorescence imaging for several methods, like STED and two-photon microscopy, in the hopes that it becomes useful for demonstrations and back-of-the-napkin feasibility calculations. Since many important physical processes are abstracted and performance is unoptimized, it does not aim to be an alternative to professionally developed solutions.

### Function

Knowledge of the physical properties of the desired fluorescing species (normalised absorption and emission spectrum, lifetime, quantum efficiency, attenuation cross-section at a give wavelength) is a prerequisite for simulation. The code aims to eventually support a variety of commonly used dyes. Adding them yourself is straightforward: in the folder *dye_spectra*, you add an unused numerical identifier (000-999) to the *properties.dat* file and specify the relevant fluorophore properties. Spectra are stored as *pkl* objects to simplify - functions which interpolate a user-supplied intensity(wavelength) file and save it as such are included. File names must follow the pattern *yyy_something_absorption.pkl* and *yyy_something_emission.pkl", yyy being the identifier. Additional information about vibrational transition lifetimes or bleaching/quenching rates may become relevant in the future, in which case property handling will be updated.

Aside from fluorophore properties, information about the excitation beam and optical setup (filter transmission spectrum, numerical aperture, detector efficiency and resolution) are needed to calculate the output. Currently the assumption is that all molecules exist in the focal plane and running the code simply results in the simulation of a single image frame.

### Current progress

When you specify the desired type and density of fluorophores, the code loads the relevant properties and spectra and first randomly generates a planar distribution of molecules. Depending on the illumination (a monochromatic laser beam is assumed), the mean number of absorbed photons per exposure time is calculated for each molecule from the absorption spectrum and local intensity.

<figure>
  <img
  src="https://github.com/timzuntar/PhoreMic/blob/main/readme_images/example_phore_map.png?raw=true"
  alt="Generated fluorophores"
  width="500">
  <figcaption></figcaption>
</figure>


From that, the code estimates the number of emission photons collected by the objective (arriving at the filter) from each fluorophore. Their wavelengths are then determined *post facto* by rejection sampling from the emission spectrum (this is by far the most computationally intensive part of the simulation). Another rejection sampler decides if the photon makes it through the filter to the detector.

<figure>
  <img
  src="https://github.com/timzuntar/PhoreMic/blob/main/readme_images/example_photon_map.png?raw=true"
  alt="Photon detection numbers"
  width="500">
  <figcaption></figcaption>
</figure>

An "image" of the single observed point is then generated as a histogram of pixel intensity values.

<figure>
  <img
  src="https://github.com/timzuntar/PhoreMic/blob/main/readme_images/example_pixel_map.png?raw=true"
  alt="Pixel map"
  width="500">
  <figcaption></figcaption>
</figure>
