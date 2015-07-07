#include "Python.h"
#include "numpy/arrayobject.h"
#include <stdio.h>

void qs(double a[],int l, int h)
{
/*
It is a quick sorting function
it could qs the array a from l to h.

example:
qs.argtypes=[POINTER(c_double),c_int,c_int]
qs.restype=c_void_p
a=range(100000)
from random import *
shuffle(a)
b=(c_double * len(a))()
for i in range(len(a)):
    b[i]=a[i]
qs(b,0,len(b)-1)
print list(b)==range(100000)
*/

  if (l >= h) 
  {
    return;
  }
  int i, j;
  double key;
  i = l;
  j = h;
  key = a[i];
  while (i < j) 
  {
    while (i < j && a[j] >= key) 
    {
      j--;
    }
    while (i < j && a[i] < key) 
    {
      i++;
    }
    if (i < j) 
    {
      double f = a[j];
      a[j] = a[i];
      a[i] = f;
    }
  }
  // printf("%d %d %d %d\n",l,h,i,j);
  if (l < j ) 
  {
    qs(a, l, j );
  }
  if (j+1 < h) 
  {
    qs(a, j+1 , h);
  }
  return;
}

static PyObject *quicksort(PyObject *self, PyObject *args, PyObject *kwargs){
    PyArrayObject *a;
    static char *kwlist[] = {"a", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O", kwlist,
                                     &a))
        return NULL;
    a = PyArray_GETCONTIGUOUS(a);
    long start=0, end = (a->dimensions[0]) -1;
    double *b = (double *) PyArray_DATA(a);
    qs(b,start,end);
    return Py_BuildValue("O", a);

}

static PyMethodDef c_algorithm_methods[] = {

    {"quicksort",  (PyCFunction)quicksort,
     METH_VARARGS | METH_KEYWORDS,
     "Return an Sorted array."},

    {NULL, NULL, 0, NULL}
};



#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef c_algorithm = {
        PyModuleDef_HEAD_INIT,
        "c_algorithm",
        "Algorithm.",
        -1,
        sort_methods,
};
PyMODINIT_FUNC PyInit_c_algorithm(void) {
    import_array();
    return PyModule_Create(&c_algorithm);
}
#else
PyMODINIT_FUNC initc_algorithm(void) {

    Py_InitModule3("c_algorithm", c_algorithm_methods,
        "Sort tools.");

    import_array();
}
#endif