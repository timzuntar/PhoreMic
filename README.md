### WORK IN PROGRESS, CODE NOT YET USABLE

## A simple fluorescence microscopy visualization tool

### Purpose

The eventual goal is to enable a (somewhat) accurate simulation and visualization of several experimental methods, like STED and two-photon fluorescence microscopy.

### Function

At minimum, properties of the chosen fluorescent species, their absorption and emission spectra, information about the excitation beam and a filter transmission spectrum are needed to calculate the output. Decisions on the wavelength of emitted photons are made by rejection sampling, which is a major computation bottleneck and will need to be optimized sooner rather than later.

### Current progress

So far, the code can randomly generate a set of positions for the desired type of fluorophores...

<figure>
  <img
  src="https://github.com/timzuntar/PhoreMic/blob/master/readme_images/example_phore_map.png?raw=true"
  alt="Generated fluorophores"
  width="500">
  <figcaption>Fluorophore locations</figcaption>
</figure>

...calculate the absorbed and emitted numbers of photons in a specified amount of time as well as the number of photons reaching the detector...

<figure>
  <img
  src="https://github.com/timzuntar/PhoreMic/blob/master/readme_images/example_photon_map.png?raw=true"
  alt="Photon detection numbers"
  width="500">
  <figcaption>Detected photon numbers</figcaption>
</figure>

...and generate an "image" of a single observed point.

<figure>
  <img
  src="https://github.com/timzuntar/PhoreMic/blob/master/readme_images/example_pixel_map.png?raw=true"
  alt="Pixel map"
  width="500">
  <figcaption>Detected image</figcaption>
</figure>
