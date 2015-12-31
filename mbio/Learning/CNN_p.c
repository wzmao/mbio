#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"
#include <stdio.h>
// #include <stdlib.h>
// #include <stdarg.h>
// #include <string.h>
// #include <ctype.h>
// #include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

static int *matrixtimes(double *a,double *b,double *c,int n,int m,int p){
    int i,j,k,tempid;

    for (i=0;i<n;i++){
        for (j=0;j<p;j++){
            tempid=i*p+j;
            c[tempid]=b[m*p+j];
            for (k=0;k<m;k++){
                c[tempid]+=a[i*m+k]*b[k*p+j];
            }
        }
    }
    return 0;
}

static PyObject *fit_ANN_BP(PyObject *self, PyObject *args, PyObject *kwargs) {

    PyArrayObject *shapedata,*inputdata,*outputdata,*transdata,**transp;
    int i,j,k,p,q,times=0,transl=0;

    static char *kwlist[] = {"shape", "input", "output", "trans", "times", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOOO|i", kwlist,
                                     &shapedata, &inputdata, &outputdata, &transdata, &times))
        return Py_BuildValue("Os", Py_None,"Couldn't parse variable from C function.");

    shapedata = PyArray_GETCONTIGUOUS(shapedata);
    inputdata = PyArray_GETCONTIGUOUS(inputdata);
    outputdata = PyArray_GETCONTIGUOUS(outputdata);
    transdata = PyArray_GETCONTIGUOUS(transdata);

    int *shape = (int *) PyArray_DATA(shapedata);
    double *input =(double *) PyArray_DATA(inputdata);
    double *output =(double *) PyArray_DATA(outputdata);
    transl=PyArray_DIM(transdata,0);
    double **trans = (double **) malloc(transl*sizeof(double *));
    if (!trans)
        return PyErr_NoMemory();
    transp=(PyArrayObject **)malloc(transl*(sizeof(PyArrayObject *)));
    if (!transp){
        free(trans);
        return PyErr_NoMemory();
    }

    if ((PyArray_NDIM(inputdata)!=2)||((PyArray_NDIM(outputdata)!=2)))
        return Py_BuildValue("Os", Py_None,"Input or Output data are not 2-D data.");
    if (PyArray_DIMS(inputdata)[0]!=PyArray_DIMS(outputdata)[0])
        return Py_BuildValue("Os", Py_None,"Input and Output dim[0] is not the same.");
    if (PyArray_DIMS(inputdata)[1]!=shape[0])
        return Py_BuildValue("Os", Py_None,"Input doesn't fit webshape.");
    if (PyArray_DIMS(outputdata)[1]!=shape[PyArray_DIMS(shapedata)[0]-1])
        return Py_BuildValue("Os", Py_None,"Output doesn't fit webshape.");

    for (i=0;i<transl;i++){
        transp[i]=((PyArrayObject **) PyArray_DATA(transdata))[i];
        transp[i]=PyArray_GETCONTIGUOUS(transp[i]);
        trans[i] = (double *)PyArray_DATA(transp[i]);
    }

    int corenumber=omp_get_max_threads(),mythreadid=-1;
    printf("* Find %d threads avaiable.\n", corenumber);

    double ***temp=NULL;

    temp=(double ***)malloc(corenumber*sizeof(double **));
    if (!temp){
        free(trans);
        free(transp);
        return PyErr_NoMemory();
    }
    for (i=0;i<corenumber;i++){
        temp[i]=(double **)malloc(transl*sizeof(double *));
        if (!temp[i]){
            for (j=0;j<i;j++){
                for (k=0;k<transl;k++)
                    free(temp[j][k]);
                free(temp[j]);
            }
            free(temp);
            free(trans);
            free(transp);
            return PyErr_NoMemory();
        }
        for (j=0;j<transl;j++){
            temp[i][j] = (double *)malloc((shape[j+1])*sizeof(double));
            if (!temp[i][j]){
                for (k=0;k<j;k++)
                    free(temp[i][j]);
                free(temp[i]);
                for (p=0;p<i;p++){
                    for (q=0;q<transl;q++)
                        free(temp[p][q]);
                    free(temp[p]);
                }
                free(temp);
                free(trans);
                free(transp);
                return PyErr_NoMemory();
            }
        }
    }

    
    // #pragma omp parallel for firstprivate(times,i,transl,input,inputdata,trans,shape,temp,j,mythreadid) schedule(dynamic,1)
    for (k=0;k<times;k++){
        mythreadid=omp_get_thread_num();
        double **nowarray=temp[mythreadid];
        for (i=0;i<transl;i++){
            if (i==0)
                matrixtimes(&input[(k%(PyArray_DIM(inputdata,0)))*PyArray_DIM(inputdata,1)],trans[0],nowarray[0],1,shape[0],shape[1]);
            else{
                matrixtimes(nowarray[i-1],trans[i],nowarray[i],1,shape[i],shape[i+1]);
            }
            for (j=0;j<shape[i+1];j++){
                nowarray[i][j]=1./(1.+exp(nowarray[i][j]));
            }
        }
        // #pragma omp barrier
        // printf("%d\n", mythreadid);
        // if (mythreadid==1){
        // }
    }
    for (i=0;i<shape[transl];i++){
        output[i]=temp[0][transl-1][i];
    }

    free(trans);
    free(transp);
    for (i=0;i<corenumber;i++){
        for (j=0;j<transl;j++)
            free(temp[i][j]);
        free(temp[i]);
    }
    free(temp);
    return Py_BuildValue("O", Py_None);
}

static PyMethodDef CNN_p_methods[] = {

    {"fit_ANN_BP",  (PyCFunction)fit_ANN_BP,
     METH_VARARGS | METH_KEYWORDS,
     "Perform BP calculation for ANN\n"},

    {NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef CNN_pmodule = {
    PyModuleDef_HEAD_INIT,
    "CNN_p",
    "Neural Networks tools with parallel.",
    -1,
    CNN_p_methods,
};
PyMODINIT_FUNC PyInit_CNN_p(void) {
  import_array();
  return PyModule_Create(&CNN_pmodule);
}
#else
PyMODINIT_FUNC initCNN_p(void) {

  Py_InitModule3("CNN_p", CNN_p_methods,
    "Neural Networks tools with parallel.");

  import_array();
}
#endif
