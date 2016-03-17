//*********************************************
 *  spiralgen_JGP.c
 *
 *  Author: Jim Pipe
 *  Date: 2011 mar 15
 *  Rev: 2011 mar 22
 *
 *  Summary: spiralgen_JGP is a support routine that takes in
 *           several design parameters in the array spparams
 *           and outputs k-space and gradient waveforms corresponding
 *           to a single 2D spiral interleaf.  The calling function 
 *           may then take this interleaf and rotate it for 
 *           multi-shot spiral trajectories.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  This code 
 *  is for research and academic purposes and is not intended for 
 *  clinical use.
 *
***********************************************/
// Author: Jim Pipe 
// Date: May 2011
// Please submit questions, comments, modifications to MRI_UNBOUND

#define MIN(a,b) (((a)>(b))?(b):(a))
#define MAX(a,b) (((a)<(b))?(b):(a))

int spiralgen(float* spparams, float **karray, float **garray, int *karrlen, int *garrlen)
{
/************************************************************/
/************************************************************

  This function takes parameters passed in spparams array and
   returns a single spiral arm calculated numerically

  The k-space trajectory for the arm is in karray
  The corresponding gradient waveform is in garray
  The parameter garrlen are the lengths of karray and garray
  The parameter karrlen indicates when sampling stops

  Units are in Hz, sec, gauss, and cm!

  grad = gm exp(i theta) i.e. gm, theta are magnitude and angle of gradient
  kloc = kr exp(i phi)   i.e. kr, phi are magnitude and angle of k-space
  alpha = theta - phi    the angle of the gradient relative to that of k-space
                         (alpha = Pi/2, you go in a circle
                          alpha = 0, you go out radially)

  The variable rfunc determines the radial spacing
  in units of the nyquist distance.
  rfunc = 1 gives critical sampling
  rfunc > 1 gives undersampling
  rfunc can vary throughout spiral generation to create variable density spirals

  KEY EQUATIONS:
  (1) dkr/dphi = rfunc*nyquist/(2 pi)
  (2) dphi/dt = gamma gm Sin(alpha)/kr
  (3) dkr/dt = gamma gm Cos(alpha)

  Solving (1)*(2) = (3) gives
  (4) Tan(alpha) = (2*pi*kr)/(rfunc*nyquist)

 ************************************************************/
/************************************************************/

  /*Initializations*/
  double d_PI = 3.141592653589793238462643383279;

  double gamma    = spparams[0]; // typically 4257.
  double fov      = spparams[1]; // in cm
  double res      = spparams[2]; // in cm
  double slewmax  = spparams[3]; // max slew rate, in gauss/cm/sec
  double greadmax = spparams[4]; // max gradient amplitude (system or bw-limited), in guass/cm
  double gsysmax  = spparams[5]; // total max gradient amplitude, in guass/cm
  double rast     = spparams[6]; // "raster time" between samples in sec, i.e. 0.000004 gives 4 usec gradient update
  int    arms     = spparams[7]; // number of spiral interleafs

  // the next 4 variables are for variable density spirals
  // they create a transition in the radial spacing as the kspace radius goes from 0 to 1, i.e.
  //    0 < kr < us_0, spacing = nyquist distance
  // us_0 < kr < us_1, spacing increases to us_r (affected by ustype)
  // us_1 < kr < 1, spacing = us_r
  double us_0    = spparams[8];
  double us_1    = spparams[9];
  double us_r    = spparams[10];
  int   ustype   = spparams[11]; // rate of change in undersampling (see code below)
  // SET us_0 = 1 or us_r = 1 to avoid variable sampling

  // some more variables determining the type of waveform
  int    gtype    = spparams[12]; // 0 = calculate through readout, 1 = include grad rampdown, 2 = include rewinder to end at k=0
  int   sptype    = spparams[13]; // 0 = Archimedean, 1 = Fermat

  double nyquist = (float)(arms)/fov; // radial distance per arm to meet the Nyquist limit
  double gamrast = gamma*rast; // gamrast*g = dk
  double dgc     = slewmax*rast; // the most the gradients can change in 1 raster period

  float *kout    = NULL;
  float *gout    = NULL;
  double *kx    = NULL;
  double *ky    = NULL;
  double *gsign    = NULL;

/* nyq_arc is for FLORET only */
/* this is the length of the arc that this set of arms will cover */
/* it is the angle, in radians, times the radius */
/* units of the length are in "nyquist distance" */
  double arc = 0.25*d_PI; /* hardcode for each hub to go to 45 degrees */
  double nyq_arc = arc*fov/res;

  double kr,kmx,kmy,kmr,rnorm,rfunc;
  double alpha, phi, theta;
  double ux,uy,gx,gy;
  double us_i;
  double gm,term;
  double kr_ramp, newslew;
  int i;

  int maxarray = 50000; /* maximum size of array  - change if necessary */
     kx = (double*) malloc(maxarray*sizeof(double));
     ky = (double*) malloc(maxarray*sizeof(double));
  gsign = (double*) malloc(maxarray*sizeof(double));

  if (kx == NULL || ky == NULL || gsign == NULL) printf ("cant allocate memory\n"); 

  for (i=0;i<maxarray;i++) gsign[i] = 1.;
  for (i=0;i<maxarray;i++) kx[i] = 0.;
  for (i=0;i<maxarray;i++) ky[i] = 0.;

/* start out spiral going radially at max slewrate for 2 timepoints */
  kx[0] = gamrast*dgc;
  ky[0] = 0.;
  kx[1] = 3.*gamrast*dgc;
  ky[1] = 0.;

  i = 1;
  kr = kx[1];

/******************************/
/* LOOP UNTIL YOU HIT MAX RES */
/******************************/
  while ((kr <= 0.5/res) && (i < maxarray-1)) {

/**************************/
/*** STEP 1:  Determine the direction (ux,uy) of the gradient at ~(i+0.5) */
/**************************/
// Dxtrapolate k at (i+0.5)
    kmx = 1.5*kx[i] - 0.5*kx[i-1];
    kmy = 1.5*ky[i] - 0.5*ky[i-1];
    kmr = sqrt(kmx*kmx + kmy*kmy);
    
// Determine the correct value for rfunc
    rnorm = 2.*res*kmr; /* the k-space radius, normalized to go from 0 to 1 */
    if (rnorm <= us_0)
      rfunc = 1;
    else if (rnorm < us_1) {
      us_i = (rnorm-us_0)/(us_1 - us_0); // goes from 0 to 1 as rnorm goes from us_0 to us_1
      if (ustype == 0) // linearly changing undersampling
        rfunc = 1. + (us_r - 1.)*us_i;
      else if (ustype == 1) // quadratically changing undersampling
        rfunc = 1. + (us_r - 1.)*us_i*us_i;
      else // assume ustype == 2, Hanning-type change in undersampling - you can add more types if you want!!!
        rfunc = 1. + (us_r - 1.)*0.5*(1.-cos(us_i*d_PI));
      }
    else
      rfunc = us_r;

  if (sptype == 1) rfunc = rfunc/(nyq_arc*rnorm); // MAKE FERMAT SPIRAL FOR FLORET

// See the Key Equation 4 at the beginning of the code
// Now estimate theta at (i+0.5)
    alpha = atan(2.*d_PI*kmr/(rfunc*nyquist));
    phi = atan2(kmy,kmx);
    theta = phi + alpha;

// This gives a unit vector in the desired direction of the next gradient point to give the right dkr/dphi
    ux = cos(theta);
    uy = sin(theta);

/**************************/
/*** STEP 2: Find largest gradient magnitude with available slew */
/**************************/

// Current gradient at (i)
    gx = (kx[i] - kx[i-1])/gamrast;
    gy = (ky[i] - ky[i-1])/gamrast;

// solve for next value of gm using the equation |gm u - g| = dgc
// which is
//   (gm u - g)(gm u* - g*) = dgc^2
// which gives
//   gm^2 (u u*) - gm (g u* + u g*) + g g* - dgc^2 = 0

// Replacing u u* with 1 (i.e. u is a unit vector) and
// replacing (g u* + u g*) with 2 Real[g u*],
// assigning b = Real[g u*]
//       and c = g g* - dgc^2,
// this is
//   gm^2 + gm (2 b) + c = 0
// giving
//   gm = -b +/- Sqrt(b^2 - c)
// The variable "term" = (b^2 - c) will be positive if we can meet the desired new gradient 

    term = dgc*dgc - (gx*gx + gy*gy) + (ux*gx + uy*gy)*(ux*gx + uy*gy);

    if (term >= 0) {

// Slew constraint can be met! Now assign next gradient and then next k value
// NOTE gsign is +1 or -1
//   if gsign is positive, we are using remaining slew to speed up (increase gm) as much as possible
//   if gsign is negative, we are using remaining slew to slow down (decrease gm) as much as possible

      gm  = MIN((ux*gx + uy*gy) + gsign[i]*sqrt(term),greadmax);
      gx = gm*ux;
      gy = gm*uy;

      kx[i+1] = kx[i] + gx*gamrast;
      ky[i+1] = ky[i] + gy*gamrast;
      i++;
      }
    else {
// We can't go further without violating the slew rate
// This means that we've sped up too fast recently to turn here at the desired curvature
// We are going to iteratively go back in time and slow down, rather than speed up, at max slew
// Thus we keep looking back for the most recent positive gsign and make it negative
// and repeat until we can make the current corner
      while ((i>3) && (gsign[i-1] == -1)) i--;
      gsign[i-1] = -1;
      i = i-2;
      }

    kr = sqrt(kx[i]*kx[i] + ky[i]*ky[i]);
    }
/***************************************************/
/******* END LOOPING FOR SAMPLING PORTION **********/
/***************************************************/

  *karrlen = i; // Done Sampling

/*****************************************************************/
/******* NOW, if requested via gtype, go to g=0 and k=0  *********/
/*****************************************************************/
  
// first we'll ramp gradients to zero
// note {ux,uy} is still pointing in the gradient direction
  if (gtype > 0) {
    while ((gm > 0) && (i < maxarray-1)) {
      gm = MAX(0,gm - dgc);
      kx[i+1] = kx[i] + gm*ux*gamrast;
      ky[i+1] = ky[i] + gm*uy*gamrast;
      i++;
      }
    }

// now, if requested via type, point gradient towards the k-space origin
// {ux,uy} will be a unit vector in that direction

  if (gtype > 1) {
    kr = sqrt(kx[i]*kx[i] + ky[i]*ky[i]);
    ux = -kx[i]/kr;
    uy = -ky[i]/kr;
    kr_ramp = 0.; // this is how far we go in k-space if we ramp down the gradient NOW

// increase gm while we can
    while ((kr_ramp < kr) && (i < maxarray-1)) {
      gm = MIN(gsysmax,gm+dgc);
      kx[i+1] = kx[i] + ux*gm*gamrast;
      ky[i+1] = ky[i] + uy*gm*gamrast;
      i++;
      kr = sqrt(kx[i]*kx[i] + ky[i]*ky[i]);
      kr_ramp = 0.5*gamma*gm*(gm/slewmax);
      }

      i--; // we went too far to technically reach k=0 with g=0, back up by 1 iteration

      ux = kx[i]/kr; // reset these to correct for numerical errors
      uy = ky[i]/kr;

// same equation as above for kr_ramp
      newslew = 0.5*gamma*gm*gm/kr;

// now decrease gm to zero
// I am forcing kr to be a function of gm (same as kr_ramp above) to make gm and kr hit 0 together
// otherwise numerical errors can make this a bit off.
    while ((kr > 0.) && (i < maxarray-1)) {
      gm = MAX(0.,gm - dgc);
      kr = 0.5*gamma*gm*(gm/newslew);
      kx[i+1] = kr*ux;
      ky[i+1] = kr*uy;
      i++;
      }
     }  

// Make sure number of gradient points is always even - take this out if you don't care
  if (i%2 == 1) i++;

// Set number of points
  *garrlen = i;

  kout    = (float*) malloc((*garrlen)*2*sizeof(float));
  gout    = (float*) malloc((*garrlen)*2*sizeof(float));

// FINAL ARRAYS are {kx(0), ky(0), kx(1), ky(1), ...} and {gx(0), gy(0), gx(1), gy(1), ...}
// g = dk / gamrast

  kout[0] = kx[0];
  kout[1] = ky[0];
  gout[0] = kx[0]/gamrast;
  gout[1] = 0;
  for (i=1;i<(*garrlen);i++) {
    kout[2*i]   = kx[i];
    kout[2*i+1] = ky[i]; 
    gout[2*i]   = (kx[i]-kx[i-1])/gamrast;
    gout[2*i+1] = (ky[i]-ky[i-1])/gamrast;
    }

  *karray = kout;
  *garray = gout;

  free (kx);
  free (ky);
  free (gsign);

  /*Thats it*/
  return(1);
}
