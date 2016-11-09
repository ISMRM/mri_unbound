/**************************************************************************
 *  GRID3_AVS.C
 *
 *  Author: Nick Zwart
 *  Date: 2011 apr 11
 *  Rev: 2011 apr 11
 *
 *  Summary: AVS module that interpolates 3-D non-cartesian samples onto a 
 *           cartesian grid. 
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  This code 
 *  is for research and academic purposes and is not intended for 
 *  clinical use.
 * 
 **************************************************************************/
#include <avs/avs.h>
#include <avs/port.h>
#include <avs/field.h>
#include <stdio.h>					/* for printf etc. */
#include <assert.h>

extern "C" 
{
  #include "threads.c"
  #include "grid_utils.c"
}

/* copy (memcpy) the data of an AVSfield_double to a dataArray_double 
 */
dataArray_double *AVSfield_2_dataArray(AVSfield_double *in)
{
    assert(in != NULL);

    /* copy input array dimensions, copy veclen as a new dimension */
    unsigned long *dim = (unsigned long*) malloc( sizeof(unsigned long)*((in->ndim)+1) );
    for(int i=0;i<(in->ndim);i++) dim[i+1] = in->dimensions[i];
    dim[0] = in->veclen;

    /* allocate new dataArray_double */
    dataArray_double *out = new_dataArray_double((in->ndim)+1, dim);
    free(dim);

    /* copy the data */
    memcpy(out->data, in->data, sizeof(double)*out->num_elem);

    return(out);
}

/* copy (memcpy) the data of a dataArray_double to an AVSfield_double 
 *  -assumes 2vec input
 */
AVSfield_double *dataArray_2_AVSfield(dataArray_double *in)
{
    assert(in != NULL);
    assert(in->dimensions[0] == 2); /* veclen */

    /* copy input array dimensions, copy veclen as a new dimension */
    int *dim = (int*) malloc( sizeof(int)*((in->nd)) );
    for(int i=0;i<(in->nd)-1;i++) dim[i] = in->dimensions[i+1];

    /* allocate new dataArray_double */
    char field_str[100]; /* generate a data type description for AVS */
    sprintf(field_str,"field %dD %d-vector uniform double", 3, 2);
    AVSfield_double *out = (AVSfield_double *) AVSdata_alloc(field_str, dim);/* let AVS allocate array */
    free(dim);

    /* copy the data */
    memcpy(out->data, in->data, sizeof(double)*(in->num_elem));

    return(out);
}




/*///////////////////////////////////////////////////////////////////////////////////// */
/*///////////////////////////////////////////////////////////////////////////////////// */
/*  BEGIN Module Compute Routine */
int module_compute( AVSfield *data, 
                    AVSfield *coords, 
                    AVSfield *weights, 
                    AVSfield *table, 
                    AVSfield **out, 
                    float *radiusFOVproduct, 
                    int width_in, 
                    float *OverSampFact, 
                    float *windowLength, 
                    float *taper_start, 
                    float *taper_length, 
                    int num_threads,
                    int veclen, /* if data is not present */
                    int compute)
{
    
    /* check to see if the body should run or not */
    if(!compute) 
    {
        /* don't signal downstream modules */
        AVSmark_output_unchanged("spectrum");  
        return(1);
    }

    /* grid size */
    int width = (float)width_in * (*OverSampFact);

    /* allocate kernel table */
    unsigned long dim[1];
    dim[0] = DEFAULT_KERNEL_TABLE_SIZE;
    dataArray_double *kern = (dataArray_double*) new_dataArray_double(1,dim);

    for(long i=0;i<kern->num_elem;i++) kern->data[i]=0.;
    loadGrid3Kernel(kern);

    /* make dataArray_double copies of input data */
    dataArray_double *data_tmp = AVSfield_2_dataArray( (AVSfield_double*)data   );
    dataArray_double *crds_tmp = AVSfield_2_dataArray( (AVSfield_double*)coords );
    dataArray_double *wght_tmp = NULL; 
    if (weights != NULL) wght_tmp = AVSfield_2_dataArray( (AVSfield_double*)weights );
    unsigned long dims[4];
    dims[0] = 2; /* complex */
    dims[1] = width;
    dims[2] = width;
    dims[3] = width;
    dataArray_double *out_tmp = new_dataArray_double(4,dims);
	
    /* choose either octant split gridding or single thread */
    if(num_threads == 8)
    {

        /* GRID3_threaded */

        /* 1) map threads */
        dataArray_double *threadMask = NULL;
        dataArray_double *centerPts  = NULL;
        int partsize = -1;
        double rfp_d = *radiusFOVproduct;
        partitionGrid_octants(crds_tmp,
                              rfp_d,
                              width,
                              &threadMask,
                              &centerPts,
                              &partsize);

        /* 2) grid it */
		grid3_threaded(data_tmp, 
                       crds_tmp, 
                       wght_tmp, 
                       out_tmp, 
                       kern, 
                       *radiusFOVproduct, 
                       *windowLength, 
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
        int nt = 1;
        int ct = 0;
        double rfp_d = *radiusFOVproduct;
        double win_d = *windowLength;
        grid3 ( &nt, &ct, 
                data_tmp, 
                crds_tmp, 
                wght_tmp, 
                &out_tmp, 
                kern, 
                NULL,
                NULL,
                &width,
                &rfp_d,
                &win_d );
    }

    /* allocate output array */
    if(*out != NULL) AVSfield_free((AVSfield*)*out);
    *out = (AVSfield *)dataArray_2_AVSfield(out_tmp);

    /* free kernel table */
    free_dataArray_double(kern);
    free_dataArray_double(out_tmp);

	return 1;
}
/*  END Module Compute Routine */
/*///////////////////////////////////////////////////////////////////////////////////// */
/*///////////////////////////////////////////////////////////////////////////////////// */


/*///////////////////////////////////////////////////////////////////////////////////// */
/*  BEGIN Module Specification */
static int module_desc()
{

	int in_port, out_port, param;

	/* Name and categorize the module. The name is automatically */
	/* generated by the filename */
	char filename[100];
	sprintf(filename, "%s", __FILE__);
	/* MODULE_DATA MODULE_FILTER MODULE_MAPPER MODULE_RENDER */
	AVSset_module_name(strtok(filename,"."), MODULE_MAPPER);

	/* INPUT PORTS ****************** */
	in_port = AVScreate_input_port("data","field 2-vector double", REQUIRED);

	in_port = AVScreate_input_port("coords","field 3-vector double", REQUIRED);

	in_port = AVScreate_input_port("weighting","field 1-vector double", REQUIRED);

	in_port = AVScreate_input_port("kernel table in","field 1D 1-vector", OPTIONAL | INVISIBLE);
	/* END INPUT PORTS ************** */

	/* OUTPUT PORTS ***************** */
	out_port = AVScreate_output_port("gridded","field double");
	/* END OUTPUT PORTS ************* */

	/* PARAMETERS ******************* */
	param = AVSadd_float_parameter("radius fov prod",DEFAULT_RADIUS_FOV_PRODUCT,0.001,FLOAT_UNBOUND);
	AVSconnect_widget(param, "typein_real");

	param = AVSadd_parameter("width","integer",160,0,INT_UNBOUND);
	AVSconnect_widget(param, "typein_integer");

	param = AVSadd_float_parameter ("OverSampFact", 1., 1., FLOAT_UNBOUND);
	AVSconnect_widget (param, "typein_real");

	param = AVSadd_float_parameter("window length",1.0,.0001,FLOAT_UNBOUND);
	AVSconnect_widget(param, "typein_real");

	param = AVSadd_float_parameter("taper start",0,0,100);
	AVSconnect_widget(param, "typein_real");

	param = AVSadd_float_parameter("taper length",0,0,100);
	AVSconnect_widget(param, "typein_real");

	param = AVSadd_parameter("number threads","integer",8,1,8);
	AVSconnect_widget(param, "typein_integer");

	param = AVSadd_parameter("veclen","integer",2,1,2);
	AVSconnect_widget(param, "typein_integer");

	param = AVSadd_parameter("compute","boolean",0,0,1);
	AVSconnect_widget(param, "toggle");

	/* END PARAMETERS **************** */

    /* send function pointer to compute module */
	AVSset_compute_proc((AVS_FNCP)module_compute);

 	/*to keep warning suppressed */
	param = in_port = out_port = 0;

	return(1);
}
/*  END Module Specification */
/*///////////////////////////////////////////////////////////////////////////////////// */

/* AVS module instantiation */
extern "C" 
{
    /* instantiate module */
    void AVSinit_modules()
    {
        AVSmodule_from_desc( module_desc ); 
    }
}


