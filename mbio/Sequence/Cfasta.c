#include "Python.h"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"
#include <string.h>

static PyObject *readFasta(PyObject *self, PyObject *args) {
	
	char *filename;
	if (!PyArg_ParseTuple(args, "s", &filename))
        return NULL;
    PyObject *labels = PyList_New(0), *seqs = PyList_New(0);
    if (!labels)
    	return PyErr_NoMemory();
    FILE *fin = fopen(filename, "r");
    long p,end=0;
    char *line=(char *)malloc(1000*sizeof(char)),*line1;
    p=fscanf(fin, "%s",&line[0]);
    while (1){
    	int l=strlen(line);
    	line1 = (line + 1);
    	PyObject *title=PyString_FromStringAndSize(line1,l-1);
    	PyList_Append(labels,title);
    	p=0;
    	while (1){
    		if (fscanf(fin, "%s", &line[0])!=1){
    			end=1;
    			break;
    		}
    		if (p==0){
    			title=PyString_FromStringAndSize(line,strlen(line));
    			p=1;
    		}
    		else{
    			if (line[0]=='>')
    				break;
    			PyObject *newpart=PyString_FromStringAndSize(line,strlen(line));
    			PyString_Concat(&title,newpart);
    		}
    	}
    	PyList_Append(seqs,title);
        if (end){
            break;
        }
    }
    fclose(fin);
    return Py_BuildValue("OO", labels, seqs);
}

static PyMethodDef Cfasta_methods[] = {

    {"readFasta",  (PyCFunction)readFasta, METH_VARARGS,
     "Read FASTA file."},

    {NULL, NULL, 0, NULL}
};


#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef Cfastamodule = {
        PyModuleDef_HEAD_INIT,
        "Cfasta",
        "FASTA IO.",
        -1,
        Cfasta_methods
};
PyMODINIT_FUNC PyInit_Cfasta(void) {
    import_array();
    return PyModule_Create(&Cfastamodule);
}
#else
PyMODINIT_FUNC initCfasta(void) {

    (void) Py_InitModule3("Cfasta", Cfasta_methods,
                          "FASTA IO.");
    import_array();
}
#endif
