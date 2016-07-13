#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"
#include <math.h>
// #include <fftw3.h>
#ifdef _OPENMP
   #include <omp.h>
#endif
#define SQR(x) ((x)*(x))
#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))


static PyObject *Cgaussian(PyObject *self, PyObject *args, PyObject *kwargs)
{
	PyArrayObject *m,*r;
	int i,j,k,l;
	double s=-10;

    static char *kwlist[] = {"matrix", "sigma", "result", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OdO", kwlist,
                                     &m, &s, &r))
        return Py_BuildValue("Os", Py_None,"Couldn't parse variable from C function.");

	m=PyArray_GETCONTIGUOUS(m);
	double *matrix = (double *) PyArray_DATA(m);
	r=PyArray_GETCONTIGUOUS(r);
	double *result = (double *) PyArray_DATA(r);
	int shape[3]={PyArray_DIM(m,0),PyArray_DIM(m,1),PyArray_DIM(m,2)};

    double ss=SQR(s);
    int sl=(int)(s*4+0.5);
    if ((sl>shape[0])|(sl>shape[1])|(sl>shape[2])){
    	Py_XDECREF(r);
    	Py_XDECREF(m);
    	result=NULL;matrix=NULL;
    	return Py_BuildValue("Os", Py_None,"Sigma is too big for the matrix.");
    }
	double *temp= (double *)malloc(shape[0]*shape[1]*shape[2]*sizeof(double));

	double *weight=(double *)malloc((sl+1)*sizeof(double)),sum=1.;
	weight[0]=1.;
	for (i=1;i<=sl;i++){
		weight[i]=exp(-0.5*SQR(i)/ss);
		sum+=2*weight[i];
	}
	for (i=0;i<sl+1;i++){
		weight[i]/=sum;
	}

	#define matrix(i,j,k) matrix[((i)*shape[1]+(j))*shape[2]+(k)]
	#define result(i,j,k) result[((i)*shape[1]+(j))*shape[2]+(k)]
	#define temp(i,j,k) temp[((i)*shape[1]+(j))*shape[2]+(k)]
	int index1,index2,index3;
	double *pp=NULL,*rpp=NULL;
	int tempstep=shape[1]*shape[2],ttstep=shape[0]*shape[1]*shape[2];
	int *indhere=(int *)malloc((sl+1)*sizeof(int)),*indhere1=(int *)malloc((sl+1)*sizeof(int));

	// filter last dimension
	for (i=0;i<=sl;i++){
		indhere1[i]=-i+shape[2];
	}
	#pragma omp parallel for schedule(guided) private(j,k,l,index1,index2,index3) firstprivate(sl,weight,pp,rpp)
	for (i=0;i<shape[0];i++){
		index1=i*tempstep;
		for (j=0;j<shape[1];j++){
			index2=index1+j*shape[2];
			for (k=0;k<sl;k++){
				index3=index2+k;
				pp=matrix+index3;
				rpp=result+index3;
				rpp[0]=weight[0]*pp[0];
				for (l=-sl;l<-k;l++)		rpp[0]+=weight[-l]*pp[indhere1[-l]];
				for (l=-k;l<0;l++)			rpp[0]+=weight[-l]*pp[l];
				for (l=1;l<=sl;l++)			rpp[0]+=weight[l]*pp[l];
			}
			for (k=sl;k<shape[2]-sl;k++){
				index3=index2+k;
				pp=matrix+index3;
				rpp=result+index3;
				rpp[0]=weight[0]*pp[0];
				for (l=1;l<=sl;l++)			rpp[0]+=weight[l]*(pp[-l]+pp[l]);
			}
			for (k=shape[2]-sl;k<shape[2];k++){
				index3=index2+k;
				pp=matrix+index3;
				rpp=result+index3;
				rpp[0]=weight[0]*pp[0];
				for (l=-sl;l<0;l++)			rpp[0]+=weight[-l]*pp[l];
				for (l=1;l<shape[2]-k;l++)	rpp[0]+=weight[l]*pp[l];
				for (l=shape[2]-k;l<=sl;l++)	rpp[0]+=weight[l]*pp[-indhere1[l]];
			}
		}
	}

	memcpy(temp,result,shape[0]*shape[1]*shape[2]*sizeof(double));
	// filter second dimension
	for (i=0;i<=sl;i++){
		indhere[i]=i*shape[2];
		indhere1[i]=indhere[i]-tempstep;
	}
	#pragma omp parallel for schedule(guided) private(j,k,l,index1,index2,index3) firstprivate(sl,weight,pp,rpp)
	for (i=0;i<shape[0];i++){
		index1=i*tempstep;
		for (k=0;k<shape[2];k++){
			index2=index1+k;
			for (j=0;j<sl;j++){
				index3=index2+j*shape[2];
				pp=temp+index3;
				rpp=result+index3;
				rpp[0]=weight[0]*pp[0];
				for (l=-sl;l<-j;l++)					rpp[0]+=weight[-l]*pp[-indhere1[-l]];
				for (l=-j;l<0;l++)						rpp[0]+=weight[-l]*pp[-indhere[-l]];
				for (l=1;l<=sl;l++)						rpp[0]+=weight[l]*pp[indhere[l]];
			}
			for (j=sl;j<shape[1]-sl;j++){
				index3=index2+j*shape[2];
				pp=temp+index3;
				rpp=result+index3;
				rpp[0]=weight[0]*pp[0];
				for (l=1;l<=sl;l++)						rpp[0]+=weight[l]*(pp[-indhere[l]]+pp[indhere[l]]);
			}
			for (j=shape[1]-sl;j<shape[1];j++){
				index3=index2+j*shape[2];
				pp=temp+index3;
				rpp=result+index3;
				rpp[0]=weight[0]*pp[0];
				for (l=-sl;l<0;l++)						rpp[0]+=weight[-l]*pp[-indhere[-l]];
				for (l=1;l<shape[1]-j;l++)				rpp[0]+=weight[l]*pp[indhere[l]];
				for (l=shape[1]-j;l<=sl;l++)			rpp[0]+=weight[l]*pp[indhere1[l]];
			}
		}
	}

	// filter the first dimension
	memcpy(temp,result,shape[0]*shape[1]*shape[2]*sizeof(double));
	for (i=0;i<=sl;i++){
		indhere[i]=i*tempstep;
		indhere1[i]=indhere[i]-ttstep;
	}
	#pragma omp parallel for schedule(guided) private(i,k,l,index1,index2,index3) firstprivate(sl,weight,pp,rpp)
	for (j=0;j<shape[1];j++){
		index1=j*shape[2];
		for (k=0;k<shape[2];k++){
			index2=index1+k;
			for (i=0;i<sl;i++){
				index3=index2+i*tempstep;
				pp=temp+index3;
				rpp=result+index3;
				rpp[0]=weight[0]*pp[0];
				for (l=-sl;l<-i;l++)					rpp[0]+=weight[-l]*pp[-indhere1[-l]];
				for (l=-i;l<0;l++)						rpp[0]+=weight[-l]*pp[-indhere[-l]];
				for (l=1;l<=sl;l++)						rpp[0]+=weight[l]*pp[indhere[l]];
			}
			for (i=sl;i<shape[0]-sl;i++){
				index3=index2+i*tempstep;
				pp=temp+index3;
				rpp=result+index3;
				rpp[0]=weight[0]*pp[0];
				for (l=1;l<=sl;l++)						rpp[0]+=weight[l]*(pp[-indhere[l]]+pp[indhere[l]]);
			}
			for (i=shape[0]-sl;i<shape[0];i++){
				index3=index2+i*tempstep;
				pp=temp+index3;
				rpp=result+index3;
				rpp[0]=weight[0]*pp[0];
				for (l=-sl;l<0;l++)						rpp[0]+=weight[-l]*pp[-indhere[-l]];
				for (l=1;l<shape[0]-i;l++)				rpp[0]+=weight[l]*pp[indhere[l]];
				for (l=shape[0]-i;l<=sl;l++)			rpp[0]+=weight[l]*pp[indhere1[l]];
			}
		}
	}
	#undef matrix
	#undef result
	#undef temp

	free(indhere);
	free(indhere1);
	free(temp);
	free(weight);
	Py_XDECREF(r);
	Py_XDECREF(m);
	result=NULL;
	matrix=NULL;
	pp=NULL;
	return Py_BuildValue("O", r);
}

static PyMethodDef Cmrc_analysis_p_methods[] = {
     
    {"Cgaussian",  (PyCFunction)Cgaussian,
     METH_VARARGS | METH_KEYWORDS,
     "Perform Gaussian filter to 3D map.\n"},

    {NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef Cmrc_analysis_pmodule = {
        PyModuleDef_HEAD_INIT,
        "Cmrc_analysis_p",
        "MRC analysis tools with parallel.",
        -1,
        Cmrc_analysis_p_methods,
};
PyMODINIT_FUNC PyInit_Cmrc_analysis_p(void) {
    import_array();
    return PyModule_Create(&Cmrc_analysis_pmodule);
}
#else
PyMODINIT_FUNC initCmrc_analysis_p(void) {

    Py_InitModule3("Cmrc_analysis_p", Cmrc_analysis_p_methods,
        "MRC analysis tools.");

    import_array();
}
#endif