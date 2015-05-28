#include "Python.h"
#include "numpy/arrayobject.h"

typedef struct MRCHeader
{   
    int nx; //number of columns (fastest changing in map)
    int ny;  //number of rows 
    int nz;  //number of sections (slowest changing in map)
    int mode;  //MODE     data type :
              // 0       image : signed 8-bit bytes range -128 to 127
              // 1       image : 16-bit halfwords
              // 2       image : 32-bit reals
              // 3       transform : complex 16-bit integers
              // 4       transform : complex 32-bit reals
              // 5       image : unsigned 8-bit range 0 to 255
              // 6       image : unsigned 16-bit range 0 to 65535
    int nxstart;  //number of first column in map (Default = 0)
    int nystart;  //number of first row in map
    int nzstart;  //number of first section in map
    int mx;  // number of intervals along X   
    int my;  //number of intervals along Y    
    int mz;  //number of intervals along Z   
    float cella[3];  //cell dimensions in angstroms   
    float cellb[3];  //cell angles in degrees       
    int mapc;  //axis corresp to cols (1,2,3 for X,Y,Z)    
    int mapr;  //axis corresp to rows (1,2,3 for X,Y,Z)    
    int maps;  // axis corresp to sections (1,2,3 for X,Y,Z)   
    float dmin;  //minimum density value    
    float dmax;  //maximum density value    
    float dmean;  //mean density value    
    int ispg;  //space group number 0 or 1 (default=0)    
    int nsymbt;  //number of bytes used for symmetry data (0 or 80)    
    char extra[100];  //extra space used for anything   - 0 by default
    float origin[3];  //origin in X,Y,Z used for transforms        
    char map[4];  //character string 'MAP ' to identify file type
    int machst;  //machine stamp    
    float rms;  //rms deviation of map from mean density    
    int nlabels;  //number of labels being used    
    char label[10][80];  //ten 80-character text labels
                          //Symmetry records follow - if any - stored as text 
                          //as in International Tables, operators separated 
                          //by * and grouped into 'lines' of 80 characters 
                          //(ie. symmetry operators do not cross the ends of 
                          //the 80-character 'lines' and the 'lines' do not 
                          //terminate in a *). 
                          //Data records follow.
} MRCHeader;

static PyObject *readHeader(PyObject *self, PyObject *args, PyObject *kwargs) {

    char *filename, SymData[80];
    PyObject *header;
    FILE *m_fp;
    MRCHeader m_header;
    int i;

    static char *kwlist[] = {"filename", "header", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "sO", kwlist,
                                     &filename, &header))
        return Py_BuildValue("Os", Py_None,"Couldn't parse variable from C function.");

    m_fp=fopen(filename,"r");

    if(m_fp==NULL)
        return Py_BuildValue("Os", Py_None,"Couldn't read file.");

    if(fread(&m_header,1,1024,m_fp)<1024)
        return Py_BuildValue("Os", Py_None,"File header is not complete.");
    
    PyObject_SetAttrString(header, "nx", PyInt_FromLong(m_header.nx));
    PyObject_SetAttrString(header, "ny", PyInt_FromLong(m_header.ny));
    PyObject_SetAttrString(header, "nz", PyInt_FromLong(m_header.nz));
    PyObject_SetAttrString(header, "mode", PyInt_FromLong(m_header.mode));
    PyObject_SetAttrString(header, "nxstart", PyInt_FromLong(m_header.nxstart));
    PyObject_SetAttrString(header, "nystart", PyInt_FromLong(m_header.nystart));
    PyObject_SetAttrString(header, "nzstart", PyInt_FromLong(m_header.nzstart));
    PyObject_SetAttrString(header, "mx", PyInt_FromLong(m_header.mx));
    PyObject_SetAttrString(header, "my", PyInt_FromLong(m_header.my));
    PyObject_SetAttrString(header, "mz", PyInt_FromLong(m_header.mz));
    PyObject_SetAttrString(header, "mapc", PyInt_FromLong(m_header.mapc));
    PyObject_SetAttrString(header, "mapr", PyInt_FromLong(m_header.mapr));
    PyObject_SetAttrString(header, "maps", PyInt_FromLong(m_header.maps));
    PyObject_SetAttrString(header, "ispg", PyInt_FromLong(m_header.ispg));
    PyObject_SetAttrString(header, "nsymbt", PyInt_FromLong(m_header.nsymbt));
    PyObject_SetAttrString(header, "machst", PyInt_FromLong(m_header.machst));
    PyObject_SetAttrString(header, "nlabels", PyInt_FromLong(m_header.nlabels));

    PyObject_SetAttrString(header, "dmin", PyFloat_FromDouble(m_header.dmin));
    PyObject_SetAttrString(header, "dmax", PyFloat_FromDouble(m_header.dmax));
    PyObject_SetAttrString(header, "dmean", PyFloat_FromDouble(m_header.dmean));
    PyObject_SetAttrString(header, "rms", PyFloat_FromDouble(m_header.rms));

    PyObject_SetAttrString(header, "map", PyString_FromStringAndSize(m_header.map,4));
    PyObject_SetAttrString(header, "extra", PyString_FromString(m_header.extra));

    for (i=0;i<3;i++){
      PyList_SetItem(PyObject_GetAttrString(header,"cella"), i, PyFloat_FromDouble(m_header.cella[i]));
      PyList_SetItem(PyObject_GetAttrString(header,"cellb"), i, PyFloat_FromDouble(m_header.cellb[i]));
      PyList_SetItem(PyObject_GetAttrString(header,"origin"), i, PyFloat_FromDouble(m_header.origin[i]));
    }

    for (i=0;i<10;i++){
      PyList_SetItem(PyObject_GetAttrString(header,"label"), i, PyString_FromString(m_header.label[i]));
    }

    if (m_header.nsymbt==0){
        PyObject_SetAttrString(header, "symdata", Py_None);
    }
    else if (m_header.nsymbt==80){
        if(fseek(m_fp, 1024, SEEK_SET)!=0)
            return Py_BuildValue("Os", Py_None,"Symmetry data couldn't be located.");
        if (fread(SymData, 80, sizeof(char), m_fp)!=0)
            return Py_BuildValue("Os", Py_None,"Couldn't parse symmetry data from file."); 
        PyObject_SetAttrString(header, "symdata", PyString_FromString(SymData));
    }
    else{
        return Py_BuildValue("Os", Py_None,"Symmetry data size is not 0 or 80.");
    }
    fclose(m_fp);
    // printf("%ld\n", PyInt_AsLong(PyObject_GetAttrString(header,"maps")));
    return Py_BuildValue("O", header);
}

static PyObject *readData(PyObject *self, PyObject *args, PyObject *kwargs) {

    char *filename=NULL;
    PyArrayObject *data;
    FILE *m_fp=NULL;
    int nsymbt=0, datamode=-1,size=0,bytesize=4;

    static char *kwlist[] = {"filename", "nsymbt", "datamode", "data", "size", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "siiOi", kwlist,
                                     &filename, &nsymbt, &datamode, &data, &size))
        return Py_BuildValue("Os", Py_None,"Couldn't parse variable from C function.");

    data = PyArray_GETCONTIGUOUS(data);
    float *c = (float *) PyArray_DATA(data);
    m_fp=fopen(filename,"r");
    if(m_fp==NULL)
        return Py_BuildValue("Os", Py_None,"Couldn't read file.");
    if(fseek(m_fp, (size_t)(1024+nsymbt), SEEK_SET)!=0)
        return Py_BuildValue("Os", Py_None,"Matrix data couldn't be located.");
    switch(datamode)
    {
        case 0:
            bytesize=1;
        case 1:
            bytesize=2;
        case 2:
            bytesize=4;
        case 5:
            bytesize=1;
        case 6:
            bytesize=2;
    }
    if(fread(c, bytesize,size, m_fp)!=size)
        return Py_BuildValue("Os", Py_None,"Parsing data Error.");
    fclose(m_fp);
    return Py_BuildValue("O", data);
}

static PyMethodDef Cmrc_methods[] = {

    {"readHeader",  (PyCFunction)readHeader,
     METH_VARARGS | METH_KEYWORDS,
     "Read the MRC file header into a header variable.\n"},

    {"readData",  (PyCFunction)readData,
     METH_VARARGS | METH_KEYWORDS,
     "Read the MRC Data into a Numpy variable.\n"},

    {NULL, NULL, 0, NULL}
};



#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef Cmrcmodule = {
        PyModuleDef_HEAD_INIT,
        "Cmrc",
        "MRC file tools.",
        -1,
        Cmrc_methods,
};
PyMODINIT_FUNC PyInit_Cmrc(void) {
    import_array();
    return PyModule_Create(&Cmrcmodule);
}
#else
PyMODINIT_FUNC initCmrc(void) {

    Py_InitModule3("Cmrc", Cmrc_methods,
        "MSA mrc tools.");

    import_array();
}
#endif
