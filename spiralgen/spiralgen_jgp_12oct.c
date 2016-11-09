/*********************************************
// Spiral Generation code
**********************************************
// Author: Jim Pipe 
// Date: May 2011
// Rev: Oct 2012
*********************************************/
// A Subset of Relevant Literature
//
// Spiral Invented:
// High-speed spiral-scan echo planar NMR imaging-I.
// Ahn, C.B., Kim, J.H. & Cho, Z.H., IEEE Transactions on Medical Imaging, 5(1) 1986.
//
// Spiral Improved:
// Fast Spiral Coronary Artery Imaging.
// Meyer CH, Hu BS, Nishimura DG, Macovski A, Magnetic Resonance in Medicine, 28(2) 1992.
//
// Variable Density Spiral
// Reduced aliasing artifacts using variable-density k-space sampling trajectories.
// Tsai CM, Nishimura DG, Magnetic Resonance in Medicine, 43(3), 2000
//
// "SLOPPY" SPIRAL
// Faster Imaging with Randomly Perturbed Undersampled Spirals and L_1 Reconstruction
// M. Lustig, J.H. Lee, D.L. Donoho, J.M. Pauly, Proc. of the ISMRM '05
//
// FLORET
// A new design and rationale for 3D orthogonally oversampled k-space trajectories
// Pipe JG, Zwart NR, Aboussouan EA, Robison RK, Devaraj A, Johnson KO, Mag Res Med 66(5) 2011
//
// Distributed Spirals
// Distributed Spirals: A New Class of 3D k-Space Trajectories
// Turley D, Pipe JG, Magnetic Resonance in Medicine, in press (also proc of ISMRM '12)

//#define DEBUG_FIELD
#ifndef MAX
#define UNDEFMAX
#define MAX(a,b) (((a)<(b))?(b):(a))
#endif

#ifndef MIN
#define UNDEFMIN
#define MIN(a,b) (((a)>(b))?(b):(a))
#endif

#define GRAST    0.005 /* Desired gradient raster time (ms) */
#define subrast    5      /* number of numerical cycles per gradient raster time */ 


// Data passed in spparams array which is spARRSIZE large

#define spARRSIZE  20

#define spGAMMA     0
#define spGMAX      1
#define spSLEWMAX   2

#define spGTYPE     3

#define spFOVXY     4
#define spFOVZ      5
#define spRESXY     6
#define spRESZ      7
#define spARMS      8

#define spSTYPE     9
#define spUSTYPE   10
#define spUS0      11
#define spUS1      12
#define spUSR      13

#define spDWELL    14
#define spREADPTS  15

#define spSLOP_PER 16

void bnispiralgen(float* spparams, int maxarray, float *gxarray, float *gyarray, float *gzarray, 
                  int *spgrad_na, int *spgrad_nb, int *spgrad_nc, int *spgrad_nd)
{
/************************************************************
************************************************************

  This function takes parameters passed in spparams array and
  returns a single spiral arm calculated numerically

  The corresponding gradient waveforms are in gxarray and gyarray
  spgrad_na reflects the number of gradient points to reach the end of k-space
  spgrad_nb = spgrad_na + the number of gradient points to ramp G to zero
  spgrad_nc = spgrad_nb + the number of gradient points to rewind k to zero
  spgrad_nd = spgrad_nc + the number of gradient points for first moment compensation

  Assignments below indicate units of input parameters
  All units input using kHz, msec, mT, and m!

  grad = gm exp(i theta) i.e. gm, theta are magnitude and angle of gradient
  kloc = kr exp(i phi)   i.e. kr, phi are magnitude and angle of k-space
  alpha = theta - phi    the angle of the gradient relative to that of k-space
                         (alpha = Pi/2, you go in a circle
                          alpha = 0, you go out radially)

  The variable rad_spacing determines the radial spacing
  in units of the Nyquist distance.
  rad_spacing = 1 gives critical sampling
  rad_spacing > 1 gives undersampling
  rad_spacing can vary throughout spiral generation to create variable density spirals

  KEY EQUATIONS:
  (1) dkr/dphi = rad_spacing*Nyquist/(2 pi)
  (2) dphi/dt = gamma gm Sin(alpha)/kr
  (3) dkr/dt = gamma gm Cos(alpha)

  Solving (1)*(2) = (3) gives
  (4) Tan(alpha) = (2*pi*kr)/(rad_spacing*Nyquist)

*************************************************************/
/* Initializations */
/************************************************************/
  double rast     = GRAST / (double)(subrast);   /* calculation "raster time" in msec */

  double gamma    = spparams[spGAMMA];   /* typically 42.577 kHz/mT */
  double gmax     = spparams[spGMAX];    /* max gradient amplitude in mT/m */
  double slewmax  = spparams[spSLEWMAX]; /* max slew rate, in mT/m/msec*/
  int    gtype    = spparams[spGTYPE];   /* 0 = calculate through readout
                                            1 = include grad ramp-down
                                            2 = include rewinder to end at k=0
                                            3 = include first moment comp */
  double fovxy    = spparams[spFOVXY];   /* enter in m */
  double resxy    = spparams[spRESXY];   /* enter in m : this should be true resolution */
  double fovz     = spparams[spFOVZ];    /* enter in m */ 
  double resz     = spparams[spRESZ];    /* enter in m : this should be true resolution */
  double arms     = spparams[spARMS];    /* number of spiral interleaves*/
  int   sptype    = spparams[spSTYPE];   /* 0 = Archimedean
                                            1 = Cylinder DST 
                                            2 = Spherical DST
                                            3 = Fermat:Floret */
  /* the next 4 variables are for variable density spirals */
  /* they create a transition in the radial spacing as the k-space radius goes from 0 to 1, i.e.*/
  /*    0 < kr < us_0 : spacing = Nyquist distance */
  /* us_0 < kr < us_1 : spacing increases to us_r (affected by ustype)*/
  /* us_1 < kr < 1    : spacing = us_r*/
  int   ustype   = spparams[spUSTYPE]; /* rate of change in undersampling
                                          0 = linear
                                          1 = quadratic
                                          2 = hanning */
  double us_0    = spparams[spUS0];
  double us_1    = spparams[spUS1];
  double us_r    = spparams[spUSR];

  /* For sloppy spirals, this lets us define periodicity in units of iteration loop time */
  /* set this to zero if you do not want sloppy spirals */
  double slop_per = spparams[spSLOP_PER];

  double nyquist = (float)(arms)/fovxy; /* radial distance per arm to meet the Nyquist limit*/
  double gamrast = gamma*rast; /* gamrast*g = dk*/
  double dgc     = slewmax*rast; /* the most the gradients can change in 1 raster period*/
  double sub_gamrast = (double)(subrast)*gamrast;
  double sub_dgc     = (double)(subrast)*dgc;

  double *kx    = NULL;
  double *ky    = NULL;
  double *kz    = NULL;
  double *gsign    = NULL;

  double kr, kmx, kmy, kmz, kmr, rnorm;
  double rad_spacing=1;
  double alpha, phi, theta;
  double ux=0,uy=0,uz=0, umag;
  double gx=0,gy=0,gz=0;
  double us_i;
  double gm=0,term;
  double gsum_ramp, gz_sum_ramp;
  double gsum, gsum0, gradtweak, gxsum, gysum, gzsum;
  double krmax, kzmax, krmax2, kzmax2;
  double krlim;
  int i, i0, i1, i_end;
  int j;

  kx    = (double*) malloc(subrast*maxarray*sizeof(double));
  ky    = (double*) malloc(subrast*maxarray*sizeof(double));
  kz    = (double*) malloc(subrast*maxarray*sizeof(double));
  gsign = (double*) malloc(subrast*maxarray*sizeof(double));

  if (kx == NULL || ky == NULL || gsign == NULL) printf ("cant allocate memory\n"); 

  for (i=0;i<subrast*maxarray;i++) gsign[i] = 1.;
  for (i=0;i<subrast*maxarray;i++) kx[i] = 0.;
  for (i=0;i<subrast*maxarray;i++) ky[i] = 0.;
  for (i=0;i<subrast*maxarray;i++) kz[i] = 0.;
  for (i=0;i<maxarray;i++) gxarray[i] = 0.;
  for (i=0;i<maxarray;i++) gyarray[i] = 0.;
  for (i=0;i<maxarray;i++) gzarray[i] = 0.;

  krmax = 0.5/resxy;
  kzmax = 0.5/resz;
  krmax2 = krmax*krmax;
  kzmax2 = kzmax*kzmax;
  krlim = krmax*(1.-(resxy/fovxy));

/* start out spiral going radially at max slew-rate for 2 time-points */
  kx[0] = 0.;
  ky[0] = 0.;
  kx[1] = gamrast*dgc;
  ky[1] = 0.;
  kx[2] = 3.*gamrast*dgc;
  ky[2] = 0.;

// IF SPHERE
  if (sptype == 2) {
    kz[0] = kzmax;
    kz[1] = sqrt(kzmax2*(1.-((kx[1]*kx[1]+ky[1]*ky[1])/krmax2))); // stay on surface of ellipsoid
    kz[2] = sqrt(kzmax2*(1.-((kx[2]*kx[2]+ky[2]*ky[2])/krmax2))); // stay on surface of ellipsoid
    }

  i = 2;
  kr = kx[2];

/******************************/
/* LOOP UNTIL YOU HIT MAX RES */
/******************************/
  while ((kr <= krlim) && (i < subrast*maxarray-1) ) {

/**************************/
/*** STEP 1:  Determine the direction (ux,uy) of the gradient at ~(i+0.5) */
/**************************/
   /* calculate dk/rast = gamma G*/
    kmx = 1.5*kx[i] - 0.5*kx[i-1];
    kmy = 1.5*ky[i] - 0.5*ky[i-1];
    kmr = sqrt(kmx*kmx + kmy*kmy);

/////////////////////////////
// Start rad_spacing logic //
/////////////////////////////
    rnorm = 2.*resxy*kmr; /* the k-space radius, normalized to go from 0 to 1 */
    
    /* determine the undersample factor */
    if (rnorm <= us_0)
      rad_spacing = 1;
    else if (rnorm < us_1) {
      us_i = (rnorm-us_0)/(us_1 - us_0); /* goes from 0 to 1 as rnorm goes from us_0 to us_1*/
      if (ustype == 0) {
/* linearly changing undersampling*/
        rad_spacing = 1. + (us_r - 1.)*us_i;
        }
      else if (ustype == 1) {
/* quadratically changing undersampling*/
        rad_spacing = 1. + (us_r - 1.)*us_i*us_i;
        }
      else if (ustype == 2) {
/* Hanning-type change in undersampling */
        rad_spacing = 1. + (us_r - 1.)*0.5*(1.-cos(us_i*M_PI));
        }
      } // if (rnorm < us_1)
    else {
      rad_spacing = us_r;
      } // rnorm > us_1

/* Undersample spiral for Spherical-Distributed Spiral */
    if (sptype == 2) {
      if(rnorm < 1.0)
        rad_spacing = MIN(fovz/resz, rad_spacing/sqrt(1.0 - (rnorm*rnorm)));
      else
        rad_spacing = fovz/resz;
      } // SDST

/* MAKE FERMAT SPIRAL FOR FLORET*/
    if (sptype == 3 && rnorm > 0.) rad_spacing *= 1./rnorm;

/* Sloppy Spirals - add variability to rad_spacing for reduced aliasing coherence */ 
// A couple different options here are commented out
// Lots of ways to be sloppy
    if (slop_per > 0) {
//      rad_spacing = MAX(1., (rad_spacing + ((rad_spacing-1.)*sin(2.*M_PI*(double)(i)/slop_per))));
//      rad_spacing += (rad_spacing-1.)*sin(2.*M_PI*slop_per*atan2(ky[i],kx[i]));
      rad_spacing += (rad_spacing-1.)*sin(2.*M_PI*slop_per*rnorm);
      }

///////////////////////////
// End rad_spacing logic //
///////////////////////////

/* See the Key Equation 4 at the beginning of the code */
    alpha = atan(2.*M_PI*kmr/(rad_spacing*nyquist));
    phi = atan2(kmy,kmx);
    theta = phi + alpha;

    ux = cos(theta);
    uy = sin(theta);

// IF SPHERICAL DST
// u dot km is zero if moving on a sphere (km is radial, u is tangential,
// thus km stays on the sphere)
// We are on an ellipsoid, but can normalize u and km by krmax and kzmax to make this work
// The final gradient vector (ux uy uz) will be tangential to the sphere
    if (sptype == 2) {
      kmz = 1.5*kz[i] - 0.5*kz[i-1];
      uz = -((ux*kmx + uy*kmy)/krmax2)*(kzmax2/kmz);
      umag = sqrt(ux*ux + uy*uy + uz*uz);
      ux = ux/umag;
      uy = uy/umag;
      uz = uz/umag;
      gz = (kz[i] - kz[i-1])/gamrast;
      }

/**************************/
/*** STEP 2: Find largest gradient magnitude with available slew */
/**************************/

/* Current gradient*/
    gx = (kx[i] - kx[i-1])/gamrast;
    gy = (ky[i] - ky[i-1])/gamrast;

/*
// solve for gm using the quadratic equation |gm u - g| = dgc
// which is
//   (gm u - g)(gm u* - g*) = dgc^2
// which gives
//   gm^2 (u u*) - gm (g u* + u g*) + g g* - dgc^2 = 0

// Replacing u u* with 1 (i.e. u is a unit vector) and
// replacing (g u* + u g*) with 2 Real[g u*]
// this is
//   gm^2 + gm (2 b) + c = 0
// giving
//   gm = -b +/- Sqrt(b^2 - c)
// The variable "term" = (b^2 - c) will be positive if we can meet the desired new gradient 
*/
    term = dgc*dgc - (gx*gx + gy*gy + gz*gz) + (ux*gx + uy*gy + uz*gz)*(ux*gx + uy*gy + uz*gz);

    if (term >= 0) {
// Slew constraint is met! Now assign next gradient and then next k value
// NOTE gsign is +1 or -1
//   if gsign is positive, we are using slew to speed up (increase gm) as much as possible
//   if gsign is negative, we are using slew to slow down (decrease gm) as much as possible
      gm  = MIN((ux*gx + uy*gy + uz*gz) + gsign[i]*sqrt(term),gmax);
      gx = gm*ux;
      gy = gm*uy;

      kx[i+1] = kx[i] + gx*gamrast;
      ky[i+1] = ky[i] + gy*gamrast;

// If SPHERE
      if (sptype == 2)
        kz[i+1] = sqrt(kzmax2*(1.-((kx[i+1]*kx[i+1]+ky[i+1]*ky[i+1])/krmax2))); // stay on surface of ellipsoid

      i++;
      } // term >= 0
    else {
// We can't go further without violating the slew rate
// This means that we've sped up too fast to turn here at the desired curvature
// We are going to iteratively go back in time and slow down, rather than speed up, at max slew
// Here we'll keep looking back until gsign is positive, then add another negative gsign, just far enough to make the current corner
      while ((i>3) && (gsign[i-1] == -1)) i--;
      gsign[i-1] = -1;
      i = i-2;
      } // term < 0

    kr = sqrt(kx[i]*kx[i] + ky[i]*ky[i]);

    } // MAIN kr loop

  i_end = i;

//********************************************
// DONE LOOPING FOR SAMPLING PORTION
// recast k to g while subsampling by subrast  
//********************************************
  gxarray[0] = 0.;
  gyarray[0] = 0.; 
  gzarray[0] = 0.; 
  gxsum = 0.;
  gysum = 0.;
  gzsum = 0.;
  for (j=1;j<=(i_end/subrast);j++) {
    i1 = j*subrast;
    i0 = (j-1)*subrast;
    gxarray[j] = (kx[i1]-kx[i0])/sub_gamrast;
    gyarray[j] = (ky[i1]-ky[i0])/sub_gamrast;
    gzarray[j] = (kz[i1]-kz[i0])/sub_gamrast;
    gxsum = gxsum + gxarray[j];
    gysum = gysum + gyarray[j];
    gzsum = gzsum + gzarray[j];
    }
  (*spgrad_na) = j;

// recalculate these ending gradient points
  gm = sqrt(gxarray[(*spgrad_na)-1]*gxarray[(*spgrad_na)-1] +
            gyarray[(*spgrad_na)-1]*gyarray[(*spgrad_na)-1] +
            gzarray[(*spgrad_na)-1]*gzarray[(*spgrad_na)-1]);
  ux = gxarray[(*spgrad_na)-1]/gm;
  uy = gyarray[(*spgrad_na)-1]/gm;
  uz = gzarray[(*spgrad_na)-1]/gm;

//**************************************************
// NOW, if requested via gtype, go to g=0 and k=0
// I've tried other ways to be faster, can't find them
//**************************************************
  
// first we'll ramp gradients to zero
// note {ux,uy} is still pointing in the gradient direction
  if (gtype > 0) {
    gz_sum_ramp = 0;
    while ((gm > 0) && (j < maxarray-1)) {
      gm = MAX(0,gm - sub_dgc);
      gxarray[j] = gm*ux;
      gyarray[j] = gm*uy;
      gzarray[j] = gm*uz;
      gxsum = gxsum + gxarray[j];
      gysum = gysum + gyarray[j];
      gzsum = gzsum + gzarray[j];
      gz_sum_ramp += gzarray[j];
      j++;
      }
    }
  *spgrad_nb = j;

// now point gradient towards the k-space origin
// {ux,uy} will be a unit vector in that direction
  if (gtype > 1) {

    /* NOTE: spherical needs a prephaser not a rewinder 
     * so just rewind x and y in that case */
    gsum = sqrt(gxsum*gxsum + gysum*gysum + gzsum*gzsum);
    if (sptype == 2 ) gsum = sqrt(gxsum*gxsum + gysum*gysum + gz_sum_ramp*gz_sum_ramp);
    gsum0 = gsum;
    ux = -gxsum/gsum;
    uy = -gysum/gsum;
    uz = -gzsum/gsum;
    if (sptype == 2) uz = -gz_sum_ramp/gsum;
    gsum_ramp = 0.5*gm*(gm/sub_dgc); /* this is *roughly* how much the area changes if we ramp down the gradient NOW*/
                                     /* this value is zero right now (gm = 0), but it will make sense below */
// increase gm while we can
    while ((gsum_ramp < gsum) && (j < maxarray-1)) {
      gm = MIN(gmax,gm+sub_dgc);
      gxarray[j] = gm*ux;
      gyarray[j] = gm*uy;
      gzarray[j] = gm*uz;
      gsum = gsum - gm;
      j++;
      gsum_ramp = 0.5*gm*(gm/sub_dgc); /* see - now this makes sense; this tells us when to start ramping down */
      }

// We've overshot it by a tiny bit, but we'll fix that later
// Ramp down for now
    while ((gm > 0) && (j < maxarray-1)) {
      gm = MAX(0,gm-sub_dgc);
      gxarray[j] = gm*ux;
      gyarray[j] = gm*uy;
      gzarray[j] = gm*uz;
      gsum = gsum - gm;
      j++;
      }
    *spgrad_nc = j;

// OK - gm is zero, but gsum is probably not EXACTLY zero. Now scale the rewinder to make the sum exactly zero
    gradtweak = gsum0/(gsum0-gsum);
    for (j=(*spgrad_nb); j<(*spgrad_nc); j++) {
      gxarray[j] = gradtweak*gxarray[j];
      gyarray[j] = gradtweak*gyarray[j];
      gzarray[j] = gradtweak*gzarray[j];
      }
  }
  free (kx);
  free (ky);
  free (kz);
  free (gsign);
}

/* undo common macro names */
#ifdef UNDEFMAX
#undef MAX
#endif

#ifdef UNDEFMIN
#undef MIN
#endif

