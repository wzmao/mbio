#include "Python.h"
// #include "stdio.h"
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

    char *filename, *mode;
    PyObject header;
    FILE *m_fp;
    MRCHeader m_header;

    static char *kwlist[] = {"filename", "mode", "header", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ssO", kwlist,
                                     &filename, &mode, &header))
        return NULL;


    m_fp=fopen(filename,mode);

    if(m_fp==NULL)
        return Py_BuildValue("OS", Py_None,"Couldn't read file.");
    rewind(m_fp);

    if(fread(&m_header,1,1024,m_fp)<1024)
        return Py_BuildValue("OS", Py_None,"File header is not complete.");
    /*unfinished*/
    return Py_BuildValue("O", header);
}


static PyMethodDef Cmrc_methods[] = {

    {"readHeader",  (PyCFunction)readHeader,
     METH_VARARGS | METH_KEYWORDS,
     "Read the MRC file header into a header variable.\n"},

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
