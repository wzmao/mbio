#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"
#include <stdio.h>
#include <stdlib.h>
// #include <stdarg.h>
// #include <string.h>
// #include <ctype.h>
// #include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif


static PyObject *fit_ANN_BP(PyObject *self, PyObject *args, PyObject *kwargs) {

    // char *filename=NULL, SymData[80], *tempchar;
    // PyArrayObject *data;
    // PyObject *header;
    // FILE *m_fp=NULL;
    // int i,j,k,bytesize=4,compress=0;

    // static char *kwlist[] = {"header", "data", "filename", "compress", NULL};

    // if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOsi", kwlist,
    //                                  &header, &data, &filename, &compress))
    //     return Py_BuildValue("Os", Py_None,"Couldn't parse variable from C function.");

    // data = PyArray_GETCONTIGUOUS(data);
    // void *matrix = (void *) PyArray_DATA(data);
    // if (compress){
    //     gzfp = gzopen(filename,"wb");
    //     if(gzfp==NULL)
    //         return Py_BuildValue("Os", Py_None,"Couldn't write file.");
    // }
    // else{
    //     m_fp=fopen(filename,"w");
    //     if(m_fp==NULL)
    //         return Py_BuildValue("Os", Py_None,"Couldn't write file.");
    // }

    // m_header.nx=PyInt_AsLong(PyObject_GetAttrString(header, "nx"));
    // m_header.ny=PyInt_AsLong(PyObject_GetAttrString(header, "ny"));
    // m_header.nz=PyInt_AsLong(PyObject_GetAttrString(header, "nz"));
    // m_header.mode=PyInt_AsLong(PyObject_GetAttrString(header, "mode"));
    // m_header.nxstart=PyInt_AsLong(PyObject_GetAttrString(header, "nxstart"));
    // m_header.nystart=PyInt_AsLong(PyObject_GetAttrString(header, "nystart"));
    // m_header.nzstart=PyInt_AsLong(PyObject_GetAttrString(header, "nzstart"));
    // m_header.mx=PyInt_AsLong(PyObject_GetAttrString(header, "mx"));
    // m_header.my=PyInt_AsLong(PyObject_GetAttrString(header, "my"));
    // m_header.mz=PyInt_AsLong(PyObject_GetAttrString(header, "mz"));
    // m_header.mapc=PyInt_AsLong(PyObject_GetAttrString(header, "mapc"));
    // m_header.mapr=PyInt_AsLong(PyObject_GetAttrString(header, "mapr"));
    // m_header.maps=PyInt_AsLong(PyObject_GetAttrString(header, "maps"));
    // m_header.ispg=PyInt_AsLong(PyObject_GetAttrString(header, "ispg"));
    // m_header.nsymbt=PyInt_AsLong(PyObject_GetAttrString(header, "nsymbt"));
    // m_header.machst=PyInt_AsLong(PyObject_GetAttrString(header, "machst"));
    // m_header.nlabels=PyInt_AsLong(PyObject_GetAttrString(header, "nlabels"));

    // m_header.dmin=PyFloat_AsDouble(PyObject_GetAttrString(header, "dmin"));
    // m_header.dmax=PyFloat_AsDouble(PyObject_GetAttrString(header, "dmax"));
    // m_header.dmean=PyFloat_AsDouble(PyObject_GetAttrString(header, "dmean"));
    // m_header.rms=PyFloat_AsDouble(PyObject_GetAttrString(header, "rms"));

    // tempchar=PyString_AsString(PyObject_GetAttrString(header, "map"));
    // strncpy(m_header.map,tempchar,4);
    // for(i=0;i<4;i++)
    //     if (m_header.map[i]=='\0'){
    //         for(j=i+1;j<4;j++)
    //             m_header.map[j]='\0';
    //         break;
    //     }

    // tempchar=PyString_AsString(PyObject_GetAttrString(header, "extra"));
    // strncpy(m_header.extra,tempchar,100);
    // for(i=0;i<100;i++)
    //     if (m_header.extra[i]=='\0'){
    //         for(j=i+1;j<100;j++)
    //             m_header.extra[j]='\0';
    //         break;
    //     }

    // for (i=0;i<3;i++){
    //   m_header.cella[i]=PyFloat_AsDouble(PyList_GetItem(PyObject_GetAttrString(header,"cella"), i));
    //   m_header.cellb[i]=PyFloat_AsDouble(PyList_GetItem(PyObject_GetAttrString(header,"cellb"), i));
    //   m_header.origin[i]=PyFloat_AsDouble(PyList_GetItem(PyObject_GetAttrString(header,"origin"), i));
    // }

    // for (i=0;i<10;i++){
    //   tempchar=PyString_AsString(PyList_GetItem(PyObject_GetAttrString(header,"label"), i));
    //   strncpy(m_header.label[i],tempchar,80);
    //     for(j=0;j<80;j++)
    //         if (m_header.label[i][j]=='\0'){
    //             for(k=j+1;k<80;k++)
    //                 m_header.label[i][k]='\0';
    //             break;
    //         }
    // }

    // if (m_header.nsymbt==80){
    //     tempchar=PyString_AsString(PyObject_GetAttrString(header, "symdata"));
    //     strncpy(SymData,tempchar,80);
    //     for(i=0;i<80;i++)
    //         if (SymData[i]=='\0'){
    //             for(j=i+1;j<80;j++)
    //                 SymData[j]='\0';
    //             break;
    //         }
    // }
    // else
    //     m_header.nsymbt=0;

    // switch(m_header.mode)
    // {
    //     case 0:
    //         bytesize=1;break;
    //     case 1:
    //         bytesize=2;break;
    //     case 2:
    //         bytesize=4;break;
    //     case 5:
    //         bytesize=1;break;
    //     case 6:
    //         bytesize=2;break;
    // }

    // // Write file.
    // if (compress){
    //     if (gzwrite(gzfp,&m_header,1024)!=1024){
    //         gzclose(gzfp);
    //         return Py_BuildValue("Os", Py_None,"Couldn't write the header.");
    //     }
    //     if (m_header.nsymbt==80){
    //         if (gzwrite(gzfp, SymData, 80)!=80){
    //             gzclose(gzfp);
    //             return Py_BuildValue("Os", Py_None,"Couldn't write Symmetry Data.");
    //         }
    //     }
    //     if (gzwrite(gzfp, matrix, bytesize*m_header.nz*m_header.ny*m_header.nx)!=bytesize*m_header.nz*m_header.ny*m_header.nx){
    //         gzclose(gzfp);
    //         return Py_BuildValue("Os", Py_None,"Couldn't write Matrix.");
    //     }
    //     gzclose(gzfp);
    // }
    // else{
    //     if (fwrite(&m_header,1,1024,m_fp)!=1024){
    //         fclose(m_fp);
    //         return Py_BuildValue("Os", Py_None,"Couldn't write the header.");
    //     }
    //     if (m_header.nsymbt==80){
    //         if (fwrite(SymData, 1, (size_t)80, m_fp)!=80){
    //             fclose(m_fp);
    //             return Py_BuildValue("Os", Py_None,"Couldn't write Symmetry Data.");
    //         }
    //     }
    //     if (fwrite(matrix, bytesize, m_header.nz*m_header.ny*m_header.nx, m_fp)!=m_header.nz*m_header.ny*m_header.nx){
    //         fclose(m_fp);
    //         return Py_BuildValue("Os", Py_None,"Couldn't write Matrix.");
    //     }
    //     fclose(m_fp);
    // }
    return Py_BuildValue("i", 0);
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
