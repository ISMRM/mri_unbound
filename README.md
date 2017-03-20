# MRI Unbound
This site is meant to provide the code from the old mri_unbound website:
http://www.ismrm.org/mri_unbound/

The projects listed below are directly reposted from the site and can be a
starting point for future development.

## PULSE SEQUENCE DESIGN

### Spiral Waveform Generation
 
#### Spiral Gen (spiralgen_jgp_12oct)
Author: Jim Pipe
License: None

A support routine that numerically generates a single (variable density) spiral
waveform, which may be rotated by the user for multi-shot spiral methods.

## RECONSTRUCTION

### Sampling Density Calculation
 
#### SDC 3D (sdc3_nrz_11aug.zip)
Author: Nick Zwart
License: None

Calculates the 3D sampling density weights for an arbitrary set of sampling
points.
 
### Gridding
 
#### Grid 3D (grid3_dct_11aug.zip)
Authors: Nick Zwart, Dallas Turley
License: None

Performs 3D gridding with an arbitrary coordinate set, data, and sampling
density weights.  The code is threaded, but can be altered slightly to run
non-threaded if desired.

### Parallel Imaging

#### Through Time Radial Grappa (throughtimeRadialGrappa_nes_13oct.zip)
Author: Nicole E. Seiberlich
License: Creative Commons

Reconstructs undersampled radial data using through-time radial GRAPPA.
