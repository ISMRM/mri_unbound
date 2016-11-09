/**************************************************************************
 *  GRID3_MAT.C
 *
 *  Author: Nick Zwart
 *  Date: 2011 apr 11
 *  Rev: 2011 apr 11
 *
 *  Summary: A MATLAB mex wrapper that implements 3D grid code.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  This code 
 *  is for research and academic purposes and is not intended for 
 *  clinical use.
 *
 **************************************************************************/

/************************************************************************** EXTERNAL LIBS */
#include "mex.h"

#include <math.h>
#include <stdlib.h>

#include "grid_utils.c" 


/************************************************************************** MEX Routine 
 * grid_volume = grid3_MAT(0:data, 1:crds, 2:weights, 3:effMtx, 4:numThreads) 
 * 
 * REQUIRED:
 *  data: N-D double-precision matrix >=2D, fastest varying dimension is length 2
 *
 *  coords: N-D double-precision matrix >=2D, fastest varying dimension is length 3
 *              trajectory coordinate points scaled between -0.5 to 0.5
 *
 *  weights: N-D double-precision matrix >=1D, size(coords)/3 == size(weights) == size(data)/2
 *          if no weighting is desired then create a matrix of the appropriate size and set values
 *          to unity.
 *
 *  effMtx: (integer) the length of one side of the grid matrix, range >=1
 *
 *  numThreads: (integer) 1 or 8, the code will run either non-threaded or threading on grid octants
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    unsigned long i;

    /* help text in the event of an input error */
    const char help_txt[] = "USAGE: (see grid3_MAT.c) \n\
 * grid_volume = grid3_MAT(0:data, 1:crds, 2:weights, 3:effMtx, 4:numThreads) \n\
 * \n\
 * REQUIRED: \n\
 *  data: N-D double-precision matrix >=2D, fastest varying dimension is length 2 \n\
 * \n\
 *  coords: N-D double-precision matrix >=2D, fastest varying dimension is length 3 \n\
 *              trajectory coordinate points scaled between -0.5 to 0.5 \n\
 * \n\
 *  weights: N-D double-precision matrix >=1D, size(coords)/3 == size(weights) == size(data)/2 \n\
 *          if no weighting is desired then create a matrix of the appropriate size and set values \n\
 *          to unity. \n\
 * \n\
 *  effMtx: (integer) the length of one side of the grid matrix, range >=1 \n\
 *          The kernel used here is designed for a 1.5 times oversampled grid (min). \n\
 * \n\
 *  numThreads: (integer) 1 or 8, the code will run either non-threaded or threading on grid octants \n";


    /* Check for proper number of arguments */
    /* input: 
    *    REQUIRED: data, coords, weight, effMtx, numThreads
    */
    if (nrhs < 3 && nrhs > 6) 
    {
        printf("%s",help_txt);
        mexErrMsgTxt("3 required inputs are: coords, numIter, effMtx.");
    }

    /* ouput: weights_out */
    if (nlhs != 1) 
    {
        printf("%s",help_txt);
        mexErrMsgTxt("Only 1 output arg is returned from grid3_MAT().");
    }

    /** check and retrieve input */
    /* PARAMS */
    int num_threads = (int) *mxGetPr(prhs[4]);
    int width       = (int) *mxGetPr(prhs[3]); /* grid size */

    /* DATA */
    printf("copy data\n");
        assert(prhs[0] != NULL);     /* check existence */
        assert(mxIsDouble(prhs[0])); /* check for type double */
    int nd = mxGetNumberOfDimensions(prhs[0]); /* get coordinate dimensions */
        assert( nd > 0 ); /* check for valid array size */
    const int *dims = mxGetDimensions(prhs[0]);
        assert( dims[0] == 2 ); /* make sure the fastest varying dim holds a complex vector [Re, Im] */
    unsigned long *dims_l = (unsigned long *) malloc(sizeof(unsigned long) * nd);
    for (i=0;i<nd;i++) dims_l[i] = (unsigned long) dims[i]; 
    dataArray_double *data = new_dataArray_double(nd,dims_l); /* alloc new dataArray_double */
        assert(data != NULL); /* make sure new mem is allocated */
    memcpy( data->data, mxGetPr(prhs[0]), sizeof(double)*(data->num_elem) ); /* copy data */
    free(dims_l);

 
    /* COORDS */
    printf("copy coords\n");
        assert(prhs[1] != NULL);     /* check existence */
        assert(mxIsDouble(prhs[1])); /* check for type double */
    nd = mxGetNumberOfDimensions(prhs[1]); /* get coordinate dimensions */
        assert( nd > 0 ); /* check for valid array size */
    dims = mxGetDimensions(prhs[1]);
        assert( dims[0] == 3 ); /* make sure the fastest varying dim holds an [x,y,z] triplet */
    dims_l = (unsigned long *) malloc(sizeof(unsigned long) * nd);
    for (i=0;i<nd;i++) dims_l[i] = (unsigned long) dims[i]; 
    dataArray_double *coords = new_dataArray_double(nd,dims_l); /* alloc new dataArray_double */
        assert(coords != NULL); /* make sure new mem is allocated */
    memcpy( coords->data, mxGetPr(prhs[1]), sizeof(double)*(coords->num_elem) ); /* copy coords */
    free(dims_l);
 

    /* WEIGHTS */
    printf("copy weights\n");
        assert(prhs[2] != NULL);     /* check existence */
        assert(mxIsDouble(prhs[2])); /* check for type double */
    nd = mxGetNumberOfDimensions(prhs[2]); /* get sample weights */
        assert( nd > 0 ); /* check for valid array size */
    dims = mxGetDimensions(prhs[2]);
    dims_l = (unsigned long *) malloc(sizeof(unsigned long) * nd);
    for (i=0;i<nd;i++) dims_l[i] = (unsigned long) dims[i]; 
    dataArray_double *weights = new_dataArray_double(nd,dims_l); /* alloc new dataArray_double */
        assert(weights != NULL); /* make sure new mem is allocated */
    memcpy( weights->data, mxGetPr(prhs[2]), sizeof(double)*(weights->num_elem) ); /* copy weights */
    free(dims_l);
 
    /* check input data sizes */
    assert( weights->num_elem == data->num_elem/2   );
    assert( weights->num_elem == coords->num_elem/3 );

    /* allocate output array */
    printf("allocate grid\n");
    unsigned long dims_g[4];
    dims_g[0] = 2; /* complex */
    dims_g[1] = width; /* effMtx */
    dims_g[2] = width;
    dims_g[3] = width;
    dataArray_double *gdata = new_dataArray_double(4,dims_g);

    /* allocate kernel table */
    printf("allocate kernel\n");
    unsigned long dim_k[1];
    dim_k[0] = DEFAULT_KERNEL_TABLE_SIZE;
    dataArray_double *kern = (dataArray_double*) new_dataArray_double(1,dim_k);
    for(long i=0;i<kern->num_elem;i++) kern->data[i] = 0.;
    loadGrid3Kernel(kern);

    /* radius FOV product */
    double rfp_d = DEFAULT_RADIUS_FOV_PRODUCT;
    double win_d = DEFAULT_WINDOW_LENGTH;

    if(num_threads == 8)
    {
        /* GRID3_threaded */
        printf("grid3 -threaded\n");

        /* 1) map threads */
        dataArray_double *threadMask = NULL;
        dataArray_double *centerPts  = NULL;
        int partsize = -1;
        partitionGrid_octants(coords,
                              rfp_d,
                              width,
                              &threadMask,
                              &centerPts,
                              &partsize);

        /* 2) grid it */
		grid3_threaded(data, 
                       coords, 
                       weights, 
                       gdata, 
                       kern, 
                       rfp_d, 
                       win_d, 
                       num_threads, 
                       threadMask, 
                       centerPts, 
                       partsize); 

        /* free tmp data */
        free_dataArray_double(threadMask);
        free_dataArray_double(centerPts);
    }
    else
    {
        /* GRID3 non-threaded */
        printf("grid3 -non-threaded\n");

        /* stub out threading vars */
        int nt = 1; 
        int ct = 0; 

        /* grid */
        grid3 ( &nt, &ct, 
                data, 
                coords, 
                weights, 
                &gdata, 
                kern, 
                NULL,
                NULL,
                &width,
                &rfp_d,
                &win_d );
    }

    /* Create an mxArray for the return argument */ 
    printf("copy output grid\n");
    int *odims = (int*) malloc( sizeof(int)*(4) );
    for(i=0;i<4;i++)
        odims[i] = dims_g[i];
    /* matlab doesn't like odims to be static? */
    plhs[0] = mxCreateNumericArray(4, odims, mxDOUBLE_CLASS, mxREAL);  
        assert(plhs[0] != NULL); /* check that mem was allocated */
    memcpy( mxGetPr(plhs[0]), gdata->data, (gdata->num_elem) * sizeof(double));
    free(odims);

    printf("free local memory\n");
    /* free temp data */
    free_dataArray_double(coords);
    free_dataArray_double(data);
    free_dataArray_double(weights);
    free_dataArray_double(gdata);

    /* free kernel table */
    free_dataArray_double(kern);
}


