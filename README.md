# MRI Unbound
This site is meant to provide the code from the old mri_unbound website:
http://www.ismrm.org/mri_unbound/

The projects listed below are directly reposted from the site and can be a
starting point for future development.

## PULSE SEQUENCE DESIGN

### Spiral Waveform Generation
 
#### spiralgen_jgp_12oct
A support routine that numerically generates a single (variable density) spiral
waveform, which may be rotated by the user for multi-shot spiral methods.

#### spiralgen_jgp_11apr (older version)
A support routine that numerically generates a single (variable density) spiral
waveform, which may be rotated by the user for multi-shot spiral methods.
 
## RECONSTRUCTION

### Sampling Density Calculation
 
#### sdc3_nrz_11aug.zip
Calculates the 3D sampling density weights for an arbitrary set of sampling
points.
 
### Gridding
 
#### grid3_dct_11aug.zip 
Performs 3D gridding with an arbitrary coordinate set, data, and sampling
density weights.  The code is threaded, but can be altered slightly to run
non-threaded if desired.

### Parallel Imaging

#### throughtimeRadialGrappa_nes_13oct.zip
Reconstructs undersampled radial data using through-time radial GRAPPA.
