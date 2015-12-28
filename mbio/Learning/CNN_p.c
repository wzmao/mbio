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


static PyObject *fit_ANN_BP(PyObject *self, PyObject *args, PyObject *kwargs) {

    PyArrayObject *shapedata,*inputdata,*outputdata,*transdata,**transp;
    int i,j,k,times=0,transl=0;

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
    if (!transp)
        return PyErr_NoMemory();    

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

    free(trans);
    free(transp);
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
