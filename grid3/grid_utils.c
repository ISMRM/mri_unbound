/**************************************************************************
 *  GRID_UTILS.C
 *
 *  Author: Nick Zwart, Dallas Turley, Ken Johnson, Jim Pipe 
 *  Date: 2011 apr 11
 *  Rev: 2011 aug 21

 * Summary: Interpolates non-cartesian k-space trajectories onto a cubic grid based on
 *          the optimized kaiser-bessel function described by Beatty et al. This code 
 *          has been threaded for speed and optimized for an 8 core architecture by 
 *          sorting the input coordinates into octants, which are independently processed 
 *          in threaded operations.  
 *
 *          A look-up table is created to speed up the kernel convolution for both threaded
 *          and non-threaded operations.
 *
 *          For roll-off correction, grid a single point at k(0,0,0), transform to image
 *          space and divide by the resulting magnitude.
 *
 *          Associated works:
 *              Jackson et al., IEEE TMI, Vol #10, pp. 473-478, 1991.
 *              Beatty et al., IEEE TMI, Vol #24, pp. 799-808, 2005
 *
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
 * This code is for research and academic purposes only and is not intended for clinical use. 
 *          
 *
 *************************************************************************************/

/*************************************************************************************
 * USAGE:
 *
 * There are two significant options for this code: threaded or non threaded. The function to 
 * implement a threaded gridding routine is grid3_threaded() as shown below. If a non-threaded 
 * grid routine is required, implement grid3() directly. The threaded function is described in 
 * detail below.  For non-threaded grid operations, the function to implement is grid3() which 
 * does not require the use of threadMask or centerPoints arrays.
 *
 * void grid3_threaded(dataArray_double *data_in,        // double nD, last dim == 2
 *                     dataArray_double *coords_in,      // nD, last dim == 3 
 *                     dataArray_double *weights_in,     // scalar, 1D
 *                     dataArray_double *grid_out,       // 4D, last dim == 2
 *                     dataArray_double *kernel_in,      // scalar, 1D, (spherical kernel)
 *                     float radiusFOVproduct,   
 *                     float windowLength,
 *                     int num_threads,                  // currently optimized for 8 threads
 *                     dataArray_int *threadMask_in,     // nD
 *                     dataArray_int *centerPts_in,      // nD last dim ==3
 *                     int partsize, // right now the partition matrices are square
 *                     int verbose)
 *
 * grid3_threaded takes in 3 data field pointers of the 'dataArray_double' type, which is 
 * a struct defined in the ARRAY UTILS section below, and returns a pointer to a dataArray_double.
 * There are checks to ensure that coords, data and weights contain the same number of elements.
 * 
 * The data inputs are defined as follows:
 *      
 *      data_in -    N-dimensional data of type double with the last dimension (or fastest varying
 *                   dimension) equal to 2, for complex data. In C, the array for a 3D spiral sequence
 *                   might have dimensions corresponding to: number of spiral arms (m), number of 
 *                   planes (n), ordered pairs (a+bi). If NULL data are passed to grid3_threaded,
 *                   the function will grid unity (useful for roll-off correction).
 *      
 *      coords_in -  REQUIRED.  N-dimensional data of type double with the last dimension equal 
 *                   to 3 for ordered tripplets corresponding to k-space coordinate points. By 
 *                   convention, it is assumed that the k-space coordinate points are between +/- 0.5
 *      
 *      weights_in - A 1-D array calculated by sample-density correction to appropriately scale the data
 *                   points for heavily over- or under- sampled k-space locations (such as the center of
 *                   k-space). If NULL weights are passed to grid3_threaded, the function will weight
 *                   all data points equally.
 *      
 *      kernel_in -  REQUIRED. The optimized kaiser-bessell function developed by Beatty et al. which will be 
 *                   convolved with the coordinate points. These are calculated by passing a 
 *                   blank dataArray_double structure with DEFAULT_KERNEL_TABLE_SIZE elements 
 *                   to loadGrid3Kernel().
 *      
 *      radiusFOVproduct - REQUIRED. Used to determine the chunk size accounting for overlapping 
 *                   chunk seams. rfp is in pixels and represents the radius of the grid kernel. 
 *                   The returned chunk will be of the determined size plus rfp+1
 *      
 *      windowLength - REQUIRED. The desired size of the resulting grid.
 *      
 *      num_threads - REQUIRED. The number of threaded processes to spawn, currently optimized 
 *                   for an 8 core processor.
 *      
 *      threadMask_in - REQUIRED.  A pointer to a 'dataArray_double' containing a bitmask corresponding
 *                   to the thread number. This is calculated by partitionGrid_octants().
 *      
 *      centerPts_in - REQUIRED.  A 3-element array contining the center points (kx,ky,kz) 
 *                   of each grid partition. This is calculated by partitionGrid_octants().
 *
 * Return:
 *      grid_out - A N-dimensional array with the fastest moving dimension equal to 2 for complex data.
 * 
 *                   
 * Arrays:
 *
 *      The dataArray_double is a struct that contains information specific to multi-dimensional
 *      arrays, such as:
 *              
 *              nd: (int) The number of dimensions.
 *
 *              dimensions: (unsigned long *) The number of points spanned by each dimension contained 
 *                          in a 1D array. The first element of the dimensions array (dimensions[0]) is 
 *                          reserved for designating if the data is scalar, complex, or ordered tripplets.
 *                          
 *              num_elem: (unsigned long) The total number of points in the array. (cumprod(dimensions))
 *
 *              data: (double *) A 1D double (or for dataArray_int, (int *) ) that contains all of the 
 *                          values for each element. 
 *
 *      The multi-dimensionality information is purely meta data that tags along with a 1D array containing 
 *      the actual data. Different macros have been defined below for appropriately accessing the 1D array,
 *      which are defined in the MACROS section.
 *
 *      A dataArray_double is allocated using new_dataArray_double(int nd, unsigned long *dimensions).
 *      A dataArray_double is deallocated using free_dataArray_double(dataArray_double *in).
 *      
 */



#ifndef GRID_UTILS_CPP
#define GRID_UTILS_CPP

/************************************************************************** EXTERNAL LIBS */
/* To turn asserts on or off */

#ifdef NDEBUG
#undef NDEBUG
#endif
/*define NDEBUG */
#include <assert.h>
#include <complex.h>
#include "threads.c"

/************************************************************************** ARRAY UTILS */
/* a data array */
typedef struct 
{
	double *data;                 /* a 1D data array of type double */
	int nd;                       /* the number of dimensions */
	unsigned long *dimensions;    /* an array that contains the size of each dimension */
    unsigned long num_elem;       /* the number of elements in this array (cumprod(dimensions)) */
} dataArray_double;

/* A method to allocate the memory for new "dataArrays". */
#define new_dataArray_double(_nd,_dim) _new_dataArray_double(_nd,_dim,__LINE__,__FILE__)
dataArray_double *_new_dataArray_double(int nd, unsigned long *dimensions,int line, const char *file)
{
    /* local vars */
    int i=0;

    /* return */
    dataArray_double *out = NULL;

    /* allocate output array struct */
	out = (dataArray_double *) malloc (sizeof(dataArray_double));
    if(out == NULL) printf("Error allocating memory: line %d in file %s\n",line,file);
	    assert (out != NULL); /* this probably means you don't have enough memory */

    /* load array information */
    out->nd = nd; /* copy number of dims */
    out->dimensions = (unsigned long *) malloc(nd*sizeof(unsigned long)); /* make a new dimensions array */
    if(out->dimensions == NULL) printf("Error allocating memory: line %d in file %s\n",line,file);
	    assert (out->dimensions != NULL); /* this probably means you don't have enough memory */
	for(i=0;i<nd;i++) out->dimensions[i] = dimensions[i]; /* copy input dimensions to struct */

    /* get the total number of elements in the array */
    out->num_elem = 1;
    for(i=0;i<nd;i++) out->num_elem *= dimensions[i];
    
    /* allocate the 1D data array */
	out->data = (double *) calloc(out->num_elem, sizeof(double)); 
    if(out->data == NULL) printf("Error allocating memory: line %d in file %s\n",line,file);
	    assert (out->data != NULL); /* this probably means you don't have enough memory */

    return(out);
}

/* A method to de-allocate the memory of a "dataArray" struct */
#define free_dataArray_double(_in) _free_dataArray_double(_in,__LINE__,__FILE__)
void _free_dataArray_double(dataArray_double *in,int line, const char *file)
{
    /* check that each input ptr is pointing to some memory */
    if(in == NULL) 
        printf("Error: invalid dataArray_double ptr: line %d in file %s\n",line,file);
    assert (in != NULL); /* struct */

    if(in->data == NULL) 
        printf("Error: invalid dataArray_double->data ptr: line %d in file %s\n",line,file);
    assert (in->data != NULL); /* data */

    if(in->dimensions == NULL) 
        printf("Error: invalid dataArray_double->dimensions ptr: line %d in file %s\n",line,file);
    assert (in->dimensions != NULL); /* dimensions */

    /* free each set of allocated data in order */
    free (in->dimensions);
    free (in->data);
    free (in);
}


/************************************************************************** KERNEL */

/* 
 *	Summary: Allocates the 3D spherically symmetric kaiser-bessel function 
 *	         for kernel table lookup.
 *  
 *	         This lookup table is with respect to the radius squared.
 *	         and is based on the work described in Beatty et al. MRM 24, 2005
 */
#define OVERSAMPLING_RATIO				1.5
#define KERNEL_WIDTH					5.0  
#define DEFAULT_RADIUS_FOV_PRODUCT		((KERNEL_WIDTH) / 2.0)
#define DEFAULT_KERNEL_TABLE_SIZE		800
#define DEFAULT_WINDOW_LENGTH			1.0
static double i0( double x )
{
	double ax = fabs(x);
	double ans;
	double y;

	if (ax < 3.75) 
    {
		y=x/3.75,y=y*y;
		ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
			   +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
	} 
    else 
    {
		y=3.75/ax;
		ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
				+y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
				+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
				+y*0.392377e-2))))))));
	}
	return (ans);
}


/* LOADGRID3KERNEL()
 * Loads a radius of the circularly symmetric kernel into a 1-D array, with
 * respect to the kernel radius squared.
 */
#define sqr(__se) ((__se)*(__se))
#define BETA (M_PI*sqrt(sqr(KERNEL_WIDTH/OVERSAMPLING_RATIO*(OVERSAMPLING_RATIO-0.5))-0.8))
#define I0_BETA	(i0(BETA))
#define kernel(__radius) (i0 (BETA * sqrt (1 - sqr(__radius))) / I0_BETA)

void loadGrid3Kernel(dataArray_double *kernTab)	
{
    /* check input data */
    assert( kernTab != NULL );
	long i;
    long size = kernTab->num_elem;

    /* load table */
	for (i=1; i<size-1; i++)	
    {
		kernTab->data[i] = kernel(sqrt(i/(double)(size-1))); /* kernel table for radius squared */
		assert (!isnan((kernTab->data[i])));
		assert (!isinf((kernTab->data[i])));
	}

    /* ensure center point is 1 */
    kernTab->data[0] = 1.0;

    /* ensure last point is zero */
    kernTab->data[size-1] = 0.0;
} /* end loadGrid3Kernel() */



/************************************************************************** GRID SUPPORT ROUTINES */

/* macros to retrieve complex data from multi-dimensional data arrays */
#define getas1p(_in,_i) ((_in)->data + (_i * (_in)->dimensions[0]))		/* returns a pointer  */
#define getas1c(_in, _i) (*(double complex *) ((_in)->data + (_i * 2)))		/* returns a complex double */
#define get3c(_in, _i, _j, _k) (*(double complex *) (((_in)->data) + \
        ((2) * (((_in)->dimensions[1] * ((_in)->dimensions[2]*(_k) + (_j))) + (_i)) ))) /* access the 1D array as if it was 3D, returns a complex double  */

/* set bounds for current data point based on kernel and grid limits */
void set_minmax (double x, int *min, int *max, int maximum, double radius)	
{
	*min = (int) ceil (x - radius);
	*max = (int) floor (x + radius);
	if (*min < 0) *min = 0;
	if (*max >= maximum) *max = maximum-1;
}


/************************************************************************** MAIN GRID ROUTINE */

/* GRID3()
 * This is the main grid routine.  It can be called directly in order to perform
 * gridding as a single process, or it is called by grid3_threaded() for 
 * parallel operation.
 */
void grid3( int *num_threads,			/* total number of threads */
            int *cur_thread,			/* this thread */
            dataArray_double *data,
            dataArray_double *coords,
            dataArray_double *weight,		/* sample weights */
            dataArray_double **out_list,	/* grid partitions for each thread */
            dataArray_double *kernel_in,	/* always in double */
            dataArray_double *threadMask_in,	/* mask of which points belong to which thread */
            dataArray_double *centerPts_in,	/* center pts of each grid partition */
            int *gridsize,			/* the final grid size in pixels */
            double *radiusFOVproduct_p,		/* kernel radius in pixels */
            double *windowLength_p ) 		/* coord scale */
{
    /* check input */
    assert(*cur_thread < *num_threads); 	/* check that num threads is valid */
    assert( coords != NULL ); 			/* check that coords exist */
    assert(*num_threads <= (coords->num_elem/3));/* check that num threads is valid */
    assert( (out_list[*cur_thread]) != NULL );  /* check that grid exists */
	assert(out_list[*cur_thread]->nd == 4); /* check that output dims are 3D */
    assert((out_list[*cur_thread])->dimensions[1] > 0); /* check for valid grid dim */
	assert(*radiusFOVproduct_p > 0.0); 	/* valid kernel radius in pix */
	assert(*windowLength_p > 0.0); 		/* valid coord scale */

    /* get rid of pointers */
	double radiusFOVproduct = *radiusFOVproduct_p, windowLength = *windowLength_p;

    /* grid size this has already been multiplied by osf */
	int width = *gridsize; /* square matrix */
    int max_x = out_list[*cur_thread]->dimensions[1];
    int max_y = out_list[*cur_thread]->dimensions[2];
    int max_z = out_list[*cur_thread]->dimensions[3];
    
    /* Calculate kernel radius in units of pixels and 1/FOV 
     *      rfp should also be appropriately scaled by oversampling factor*/
	double kernelRadius           = radiusFOVproduct / (double)width;
	double kernelRadius_sqr       = kernelRadius * kernelRadius;
	double kernelRadius_invSqr    = 1.0 / (kernelRadius * kernelRadius);

    /* access kernel table from kr */
	double dist_multiplier, dist_sqr;
	double width_inv = 1.0 / width;

    /* scale grid point in pixels to match table size */
	dist_multiplier = (kernel_in->dimensions[0] - 1) * kernelRadius_invSqr;

    /* grid points to check */
	int imin, imax, jmin, jmax, kmin, kmax, i, j, k;
	double x, y, z, ix, jy, kz;

    /* kr */
	double dx_sqr, dy_sqr, dz_sqr, dz_sqr_PLUS_dy_sqr;

    /* zero output field */
    for(unsigned long ii=0;ii<out_list[*cur_thread]->num_elem;ii++)
        out_list[*cur_thread]->data[ii] = 0.;

    /* get center pt for current partition */
    int center_x = width/2;
    int center_y = width/2;
    int center_z = width/2;

    /* use partition specific center points for multi-thread */
    if(centerPts_in != NULL)
    {
        center_x = centerPts_in->data[(*cur_thread)*3  ];
        center_y = centerPts_in->data[(*cur_thread)*3+1];
        center_z = centerPts_in->data[(*cur_thread)*3+2];
    }

    /* grid main loop */
    for (unsigned long p=0; p<(coords->num_elem/3); p++)
    {
        /* check if this point is on the current thread */
        if( threadMask_in != NULL )
        if( (threadMask_in->data[p]) != (*cur_thread) ) continue;

        /* get data */
        double complex dat = getas1c(data,p);

        /* weight the data */
        if (weight != NULL) dat *= weight->data[p];

		/* get the coordinates of the datapoint to grid */
		x = (coords->data[p*3  ]) / windowLength; /* these usually vary between -.5 -- +.5 */
		y = (coords->data[p*3+1]) / windowLength; /* these usually vary between -.5 -- +.5 */
		z = (coords->data[p*3+2]) / windowLength; /* these usually vary between -.5 -- +.5 */

		/* set the boundaries of final dataset for gridding this point */
        ix = x * width + center_x;
        set_minmax(ix, &imin, &imax, max_x, radiusFOVproduct);
        jy = y * width + center_y;
        set_minmax(jy, &jmin, &jmax, max_y, radiusFOVproduct);
        kz = z * width + center_z;
        set_minmax(kz, &kmin, &kmax, max_z, radiusFOVproduct);

		/* grid this point onto the neighboring cartesian points */
        for (k=kmin; k<=kmax; k++)	
        {
            kz = (k - center_z) * width_inv;
            dz_sqr = kz - z;
            dz_sqr *= dz_sqr;
            for (j=jmin; j<=jmax; j++)	
            {
                jy = (j - center_y) * width_inv;
                dy_sqr = jy - y;
                dy_sqr *= dy_sqr;
                dz_sqr_PLUS_dy_sqr = dz_sqr + dy_sqr;
                if (dz_sqr_PLUS_dy_sqr < kernelRadius_sqr)	
                {
                    for (i=imin; i<=imax; i++)	
                    {
                        ix = (i - center_x) * width_inv;
                        dx_sqr = ix - x;
                        dx_sqr *= dx_sqr;
                        dist_sqr = dx_sqr + dz_sqr_PLUS_dy_sqr;
                        if (dist_sqr < kernelRadius_sqr)	
                        {
                            /* get kernel value */
                            double complex val = kernel_in->data[(int) round (dist_sqr * dist_multiplier)];

                            /* multiply data by current kernel val */
                            val *= dat;

                            /* grid complex or scalar */
                            get3c(out_list[*cur_thread], i,j,k) += val;
                        } /* kernel bounds check, spherical support */
                    } /* x 	 */
                } /* kernel bounds check, spherical support */
            } /* y */
        } /* z */

    } /* data point */

} /* end grid3() */



/************************************************************************** THREADED GRID SUPPORT ROUTINES */

/* Split grid coordinates into octants for threaded processing
 * and return an array containing start and end locations for 
 * each axis (x,y,z) for each chunk.
 *
 * Returns a chunk mask (map) for each thread relating the coordinate to the thread.
 *
 * The radius-field-fov product is used to determine the chunk size accounting
 * for overlapping chunk seams.  rfp is in pixels and represents the radius of 
 * the grid kernel.  The returned partition will be of the an eigth of the grid size 
 * plus rfp + 1.
 *
 * gridsize = len(x) = len(y) = len(z)
 */

/* GETOCTANTFROMCOORD()
 * For coordinate points that fall exactly on zero.
 * octant boundaries for <x,y,z>, (greater than) or 
 * (less than or equal to) to zero:
 *  0) > ,> ,>    bitmask: 0,0,0
 *  1) <=,> ,>    bitmask: 1,0,0
 *  2) > ,<=,>    etc...
 *  3) <=,<=,>
 *  4) > ,> ,<=
 *  5) <=,> ,<=
 *  6) > ,<=,<=
 *  7) <=,<=,<=
 */
int getOctantFromCoord(double *crd) /* crd must be an ordered tripplet */
{
    /* bitmask for octant number (eg. 111 == 7) */
    int bitmask = 0;
    
    /* x, bit 2^2 = 4 */
    if( crd[0] <= 0. ) bitmask += 4;

    /* y, bit 2^1 = 2 */
    if( crd[1] <= 0. ) bitmask += 2;

    /* z, bit 2^0 = 1 */
    if( crd[2] <= 0. ) bitmask += 1;

    return(bitmask); 
}


/* PARTITIONGRID_OCTANTS()
 * Sorts coordinate points for threaded gridding.
 */
void partitionGrid_octants(dataArray_double *crds_in, 
                           double rfp_in, 
                           int gridsize_in,
                           dataArray_double **threadMask_out, 
                           dataArray_double **partitionsCenterPixel_out, 
                           int *partsize_out)
{
    int numPartitions = 8; /* octants     */
    int rfp = (int)ceil(rfp_in);
    unsigned long dim[2];

    /** label each crd with a corresponding thead */
    /* allocate memory */
    if(*threadMask_out != NULL) free_dataArray_double(*threadMask_out);
    dim[0] = 3;
    dim[1] = (crds_in->num_elem/3);
    *threadMask_out = new_dataArray_double(2,dim);

    /* assign thread numbers based on which partition they fall into */
    for(unsigned long ii=0;ii<(crds_in->num_elem/3);ii++)
        (*threadMask_out)->data[ii] = (double) getOctantFromCoord(crds_in->data + ii*3);

    /** get list of center pixels for each partition */
    /* allocate memory */
    if(*partitionsCenterPixel_out != NULL) free_dataArray_double(*partitionsCenterPixel_out);
    dim[0] = 3;
    dim[1] = numPartitions;
    *partitionsCenterPixel_out = new_dataArray_double(2,dim); 

    /* Assign center pixel for each octant. 
     * If the coords for that octant are '>', then the center
     * pixel is rfp+1.  If the coords for that octant are '<=',
     * then the center pt is effmtx/2.
     */
    for(int i=0;i<numPartitions;i++)         
    {
        /* x, bit 2^2 = 4 */
        if( (i & 4) == 4 ) getas1p(*partitionsCenterPixel_out,i)[0] = gridsize_in/2; /* <= */
        else               getas1p(*partitionsCenterPixel_out,i)[0] = rfp + 1; /* > */

        /* y, bit 2^1 = 2 */
        if( (i & 2) == 2 ) getas1p(*partitionsCenterPixel_out,i)[1] = gridsize_in/2; /* <= */
        else               getas1p(*partitionsCenterPixel_out,i)[1] = rfp + 1; /* > */

        /* z, bit 2^0 = 1 */
        if( (i & 1) == 1 ) getas1p(*partitionsCenterPixel_out,i)[2] = gridsize_in/2; /* <= */
        else               getas1p(*partitionsCenterPixel_out,i)[2] = rfp + 1; /* > */
    }

    /* assign parition size */
    *partsize_out = gridsize_in/2 + rfp+1;
}



/* _ADDPARTITIN3_THREAD()
 * Threaded partition add for summing individual 3D partitions/octants
 * to the final output grid mtx.
 */
void _addPartition3_thread( int *num_threads,
                            int *cur_thread,
                            dataArray_double *part,
                            dataArray_double *fullgrid,
                            dataArray_double *centerPts_in,
                            int *part_num)
{
    /* grid sizes */
    int width = fullgrid->dimensions[1];
    int wd2 = width/2;
    int part_zwidth = part->dimensions[3];
    int partsize = part_zwidth; /* cubic grid */

    /* split threads up by even chunks of data_out memory */
    int k_start =  *cur_thread    * part_zwidth / *num_threads;
    int k_end   = (*cur_thread+1) * part_zwidth / *num_threads;

    /* start threaded loop */
    for(int k=k_start;k<k_end;k++)
    for(int j=0;j<partsize;j++)
    for(int i=0;i<partsize;i++)
    {
        int i_off = (int) centerPts_in->data[(*part_num)*3  ];
        int j_off = (int) centerPts_in->data[(*part_num)*3+1];
        int k_off = (int) centerPts_in->data[(*part_num)*3+2];

        int ii = i + wd2 - i_off;
        int jj = j + wd2 - j_off;
        int kk = k + wd2 - k_off;

        get3c(fullgrid,ii,jj,kk) += get3c(part,i,j,k);
    }
} /* end - _addPartition3_thread() */


/* GRID3_THREADED()
 * The calling function for a threaded grid.
 */
void grid3_threaded ( dataArray_double *data_in,        /* NULL, nD, last dim == 2 */
                      dataArray_double *coords_in,      /* nD, last dim == 3 */
                      dataArray_double *weights_in,     /* NULL or scalar */
                      dataArray_double *grid_out,       /* 4D, last dim == 2 (for complex data */
                      dataArray_double *kernel_in,      /* 1D (spherical kernel) */
                      double radiusFOVproduct,   
                      double windowLength,		/* desired grid size  */
                      int num_threads,
                      dataArray_double *threadMask_in,  /* nD */
                      dataArray_double *centerPts_in,   /* nD last dim == 3 */
                      int partsize) /* right now the partition matrices are square */
{
    /* check that input arrays are initialized */
    assert (coords_in != NULL);
    assert (grid_out != NULL);
    assert (kernel_in != NULL);
    assert (threadMask_in != NULL);
    assert (centerPts_in != NULL);

    /* Create a new volume to grid to
     * for each thread.
     *
     * Allocate an array of dataArray_double ptrs.
     */
    dataArray_double **out_list = (dataArray_double **) malloc( sizeof(dataArray_double *)*num_threads);

    /* allocate a partition for each thread */
    unsigned long dim[4];
    dim[0] = 2; /* complex */
    dim[1] = partsize; /* x */
    dim[2] = partsize; /* y */
    dim[3] = partsize; /* z */
    for(int i=0;i<num_threads;i++) 
        out_list[i] = (dataArray_double *) new_dataArray_double(4,dim);

    /* start threads */
    int gridsize = grid_out->dimensions[1];
    create_threads10( num_threads, 
                      grid3, 
                      data_in, 
                      coords_in, 
                      weights_in, 
                      out_list, 
                      kernel_in,
                      threadMask_in,
                      centerPts_in, 
                      &gridsize,
                      &radiusFOVproduct, 
                      &windowLength);

    /* zero out the grid  */
    for(unsigned long ii=0;ii<grid_out->num_elem;ii++)
        grid_out->data[ii] = 0.;

    /* add each parition to the grid */
    for(int p=0;p<num_threads;p++)
    {
        /* threaded-add each partition */
        create_threads4( num_threads,
                        _addPartition3_thread,
                        out_list[p], 	/* current partition */
                        grid_out,
                        centerPts_in,
                        &p );		 	/* current thread */
    }

    /* deallocate a field for each ptr */
    for(int i=0;i<num_threads;i++) free_dataArray_double(out_list[i]);
    free(out_list);

} /* end grid3_threaded() */



#endif /* GRID_UTILS_CPP -EOF */
