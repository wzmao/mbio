#include "Python.h"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"
#include <math.h>
#include "fftw3.h"
#include <omp.h>


typedef struct MRCHeader
{	
	int nx;	//number of columns (fastest changing in map)
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


typedef struct CMPLXd
{
	double x;
	double y;
} CMPLXd;

typedef struct CMPLXf
{
	float x;
	float y;
} CMPLXf;

static CMPLXf GetStructureFactor(int *atomid,float *occ,float *bf,float *cor,int atomnumber, float h, float k, float l)
{
	int i;
	float PI=3.14159265359;
	float a1[]={0.04924,0.06056,0.18100,0.18017,0.18848,0.23238,0.13897,
	0.24115,0.23422,0.23449,0.54110,0.52369,0.56914,0.53363,0.50346,
	0.49157,0.48380,0.47103,0.83435,0.81886,0.82184,0.83718,0.85465,
	0.88038,0.88646,0.90628,0.90998,0.93017,0.95193,0.95958,1.09595,
	1.07488,1.04029,1.02041,1.00375,0.98899,1.42780,1.38199,1.34006,
	1.30145,1.25788,1.23849,1.24122,1.21842,1.20244,1.05766,1.21087,
	1.22503,1.35495,1.33656,1.30727,1.29090,1.30488,1.42532,2.21883,
	2.19676,2.13774,2.14007,2.19312,2.18399,2.17908,2.18332,2.17938,
	2.11909,2.19660,2.16544,2.04451,2.15671,2.15114,2.14582,2.07541,
	2.03115,1.98850,1.97392,1.92672,1.90732,1.87864,1.84670,1.82271,
	1.83029,1.93415,1.94124,1.92664,1.88230,1.86076,1.86088,2.60915,
	2.61407,2.56421,2.50061,2.55227,2.56086,2.53867,2.57218,2.57125,
	2.53626,2.54045,2.25859};
	float b1[]={1.12357,0.69147,1.37982,1.11557,0.97203,1.01435,0.52475,
	0.80966,0.70369,0.63245,1.19898,1.08412,1.10085,0.97815,0.87932,
	0.82092,0.77733,0.72503,1.26295,1.19024,1.14499,1.11535,1.08878,
	1.07250,1.03534,1.01383,0.97791,0.96164,0.94385,0.91866,0.99943,
	0.95257,0.89757,0.85673,0.82180,0.78947,1.13570,1.07975,1.03163,
	0.98560,0.93706,0.90855,0.90093,0.87173,0.85021,0.72821,0.83473,
	0.83438,0.91033,0.88439,0.85098,0.82652,0.82164,0.87668,1.28716,
	1.25255,1.20039,1.17941,1.18570,1.16032,1.13818,1.12260,1.10027,
	1.05620,1.07789,1.04740,0.97532,1.01427,0.99878,0.98379,0.93924,
	0.90829,0.87840,0.86344,0.83350,0.81707,0.79535,0.77396,0.75692,
	0.75599,0.79990,0.79971,0.78860,0.76281,0.74794,0.74375,1.06688,
	1.05942,1.02943,0.99429,1.00772,1.00272,0.98578,0.99197,0.98329,
	0.96055,0.95368,0.83727};
	float a2[]={0.13162,0.15040,0.46885,0.56677,0.64565,0.81994,0.48164,
	0.66042,0.62045,0.56495,0.90328,0.89281,1.06080,1.08891,1.09013,
	1.13064,1.16254,1.18116,2.57342,2.38629,2.29740,2.22146,2.14893,
	2.09996,1.99407,1.92631,1.85730,1.78977,1.74900,1.67373,1.72261,
	1.67658,1.61472,1.56993,1.52874,1.53397,3.42834,3.33289,3.27025,
	3.22818,3.20895,3.20473,3.20720,3.19676,3.16863,2.69010,3.16064,
	3.12183,3.24016,3.11583,2.98093,2.85075,2.76224,2.77836,4.60571,
	4.67806,4.59844,4.55426,4.52007,4.42520,4.35491,4.28398,4.22665,
	4.09793,4.11165,3.98593,3.74949,3.83978,3.76811,3.70152,3.58591,
	3.51892,3.48239,3.48896,3.44733,3.45952,3.46471,3.46733,3.47527,
	3.53807,3.85934,3.85731,3.82543,3.71207,3.63868,3.58711,5.08860,
	5.18681,5.10621,4.98823,5.16312,5.20664,5.15760,5.22851,5.21316,
	5.10320,5.07297,4.49760};
	float b2[]={4.81011,3.08412,10.0558,7.28333,5.81321,5.78949,2.68683,
	3.71987,3.09136,2.64199,6.49087,5.89505,6.73577,5.88225,5.13366,
	4.73512,4.38180,4.00402,8.16010,7.17228,6.56826,6.13917,5.78576,
	5.53692,5.18729,4.97780,4.71826,4.56669,4.46409,4.30323,5.03468,
	4.78250,4.43901,4.20859,4.01238,3.88699,8.55214,7.78631,7.15377,
	6.62142,6.13119,5.80547,5.60056,5.28704,5.01511,4.00891,4.69018,
	4.55327,4.91209,4.61908,4.30739,4.05350,3.93610,4.20012,8.64777,
	8.37051,7.83720,7.64689,7.71404,7.46517,7.27662,7.13181,6.98080,
	6.54065,6.79783,6.50638,5.79447,6.24679,6.12492,6.01377,5.61207,
	5.35170,5.12588,5.02772,4.79093,4.68197,4.53871,4.38830,4.26626,
	4.26088,4.66429,4.60863,4.47990,4.23226,4.07523,3.98911,7.01375,
	6.96304,6.62342,6.24105,6.35932,6.28157,6.06470,6.06710,5.94195,
	5.68828,5.57817,4.52510};
	float a3[]={0.22393,0.15627,1.47726,1.44076,1.33187,1.12081,1.06113,
	0.79431,0.71679,0.63292,1.60421,2.10328,2.67657,2.70558,2.58035,
	2.42939,2.26346,2.09559,2.42561,3.38497,3.22354,3.04722,2.89030,
	2.20830,2.54856,2.41077,2.28316,2.16054,1.64503,1.95071,2.61438,
	2.87946,2.98528,3.02495,2.99457,2.97758,2.98652,4.00265,4.12802,
	4.06357,3.55467,3.44083,3.64441,3.12230,2.96386,2.64994,2.66163,
	2.86597,3.44836,3.80499,4.04458,4.19425,4.38534,4.55840,4.43838,
	5.20415,5.46274,5.28185,4.71689,4.54596,4.42783,4.31155,4.23139,
	4.45265,4.09941,3.90800,4.00012,3.74202,3.65004,3.58929,3.87601,
	4.01070,4.07937,4.08018,4.03028,3.97995,3.91500,3.61207,3.51273,
	3.60589,3.92418,4.20251,4.52248,4.73640,4.97434,5.18974,5.62536,
	6.02657,6.43140,6.70978,6.10216,5.93139,5.71736,5.13516,4.98168,
	5.21947,5.06255,4.74670};
	float b3[]={14.8439,9.62423,45.2163,26.6073,20.0682,19.2115,9.74277,
	11.1230,9.44919,7.79134,40.7770,29.8227,29.0652,23.3404,18.8552,
	16.1274,14.0357,12.2894,45.7999,39.3149,34.5258,31.4668,29.2795,
	26.8071,25.4956,24.2683,22.9255,21.9676,21.0681,20.4870,24.7417,
	22.0912,19.0809,16.8080,14.8316,13.5568,41.4786,39.6810,34.5591,
	30.6316,26.3035,24.5162,24.3207,21.7096,20.3657,13.5298,19.1689,
	19.4157,23.6490,22.0613,20.0491,18.1887,17.1101,17.1970,37.6962,
	39.1838,35.2049,34.6029,36.6823,35.5218,35.0384,34.6951,34.4883,
	30.8713,34.6966,32.6372,26.7471,31.8179,31.3490,31.0735,27.8278,
	25.4086,23.4938,22.2866,20.4165,19.3994,18.3615,16.7894,15.9633,
	16.2594,19.7903,19.8331,19.4982,18.1512,17.3069,16.6754,31.4890,
	32.0739,30.0266,27.5340,28.2608,27.9933,26.8172,27.2516,26.8301,
	25.4559,25.0739,18.4511};
	float a4[]={0.11705,0.04248,1.13754,0.84079,0.60166,0.30422,0.50718,
	0.25232,0.19442,0.18389,1.66459,1.62238,1.50965,1.42748,1.24278,
	1.03713,0.87424,0.75940,3.03079,3.20375,2.84374,2.54718,2.28594,
	1.65155,1.94847,1.79055,1.67125,1.55296,1.11668,1.34175,1.52306,
	1.59500,1.53121,1.44175,1.38444,1.24952,3.72877,4.18667,3.72874,
	3.36196,2.44883,2.17559,2.55125,1.80913,1.69441,0.99422,1.42689,
	1.80448,2.15730,2.37062,2.41286,2.43357,2.22441,1.79305,4.92676,
	5.87447,5.28736,5.08329,5.23029,5.12534,4.95194,4.78637,4.60577,
	4.27187,4.22118,4.24760,4.26511,4.00651,3.90289,3.78259,3.61559,
	3.29314,2.98409,2.68047,2.54236,2.32049,2.14741,1.57926,1.45598,
	1.67871,2.05207,2.25982,2.48360,2.69494,2.66395,2.51993,4.93145,
	6.27116,5.92985,5.46739,5.29750,4.97860,4.89271,4.79656,4.61062,
	4.39750,4.27205,5.03828};
	float b4[]={39.0726,23.5306,128.547,73.1093,57.3252,54.1132,30.7532,
	30.3521,26.2880,21.7321,129.803,88.1764,90.5010,71.1450,54.9253,
	45.3868,38.5357,33.0660,166.443,121.103,108.232,100.008,94.2114,
	98.1169,83.2112,79.4409,75.6734,72.4782,79.6493,67.4006,84.9517,
	71.1665,57.6297,48.6300,41.5922,37.1324,173.152,130.817,112.588,
	101.129,96.6631,92.2924,83.6231,84.2304,80.6371,44.4972,78.9906,
	69.6184,85.1136,73.8053,62.3919,53.7125,48.8038,46.9535,189.044,
	145.486,126.133,123.673,134.769,131.059,128.845,126.872,125.712,
	109.416,124.661,118.443,96.1710,114.901,112.622,111.822,96.1033,
	86.5407,80.1163,75.7871,69.5880,66.0742,62.9838,60.4880,57.7965,
	56.3570,74.7244,69.2077,64.7147,57.2473,52.3519,48.6098,171.743,
	134.285,114.588,100.577,107.961,106.931,103.069,112.848,110.605,
	97.5280,95.3855,76.9542};
	float c[]={0.00683,0.00819,0.02026,0.02330,0.02600,0.03179,0.02179,
	0.03513,0.03534,0.03567,0.06161,0.06365,0.07005,0.06988,0.06975,
	0.07118,0.07318,0.07339,0.11302,0.11416,0.11645,0.11901,0.12164,
	0.12471,0.12617,0.12842,0.12965,0.13262,0.13474,0.13659,0.14866,
	0.14870,0.14722,0.14656,0.14669,0.14660,0.19685,0.19625,0.19642,
	0.19601,0.19452,0.19520,0.19966,0.19977,0.20136,0.18384,0.20701,
	0.21128,0.22698,0.22688,0.22472,0.22396,0.22646,0.23859,0.29953,
	0.29951,0.29740,0.29857,0.30398,0.30473,0.30577,0.30792,0.30810,
	0.30556,0.31428,0.31329,0.30399,0.31613,0.31797,0.31964,0.32463,
	0.31215,0.30946,0.31075,0.30761,0.30798,0.30603,0.30430,0.30394,
	0.30863,0.32681,0.33191,0.33331,0.32945,0.32869,0.33168,0.42517,
	0.42739,0.42435,0.41959,0.42682,0.42920,0.42879,0.43444,0.43596,
	0.43314,0.43449,0.40399};

	CMPLXd Fhkl;
	CMPLXf tmp;

	float r2=h*h+k*k+l*l;

	Fhkl.x=0.0;
	Fhkl.y=0.0;

	int curAtomID=-1;
	float AtomDiffFactor=0.0;
	float amp,phase;

	for(i=0;i<atomnumber;i++)
	{
		if(atomid[i]!=curAtomID) 
		{
			curAtomID=atomid[i];
			AtomDiffFactor=(a1[curAtomID]*exp(-(b1[curAtomID])*r2/4.0)+a2[curAtomID]*exp(-(b2[curAtomID])*r2/4.0)+a3[curAtomID]*exp(-(b3[curAtomID])*r2/4.0)+a4[curAtomID]*exp(-(b4[curAtomID])*r2/4.0)+c[curAtomID]);
		}
		amp=exp(-r2*bf[i]/2.0)*occ[i]*AtomDiffFactor;
		phase=-2.0*PI*(h*cor[i*3]+k*cor[i*3+1]+l*cor[i*3+2]);
		Fhkl.x+=cos(phase)*amp;
		Fhkl.y+=sin(phase)*amp;
	}
	tmp.x=Fhkl.x;
	tmp.y=Fhkl.y;
	return tmp;
}

static void ifft3d(float* buf, int nsam, float* out)
{
	fftwf_plan plan_fft=fftwf_plan_dft_c2r_3d(nsam,nsam,nsam,(fftwf_complex *) buf,out,FFTW_ESTIMATE);  
	fftwf_execute(plan_fft);
	fftwf_destroy_plan(plan_fft);
}

static PyObject *Cpdb2mrc(PyObject *self, PyObject *args, PyObject *kwargs)
{
	PyArrayObject *atomiddata, *occdata,*bfdata, *cordata, *map;
	int nsam;
	float psize,res;
    static char *kwlist[] = {"atomid", "occ", "bf", "cor", "map", "nsam", "psize", "res", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOOOOiff", kwlist,
                                     &atomiddata, &occdata, &bfdata, &cordata, &map, &nsam, &psize, &res))
        return Py_BuildValue("Os", Py_None,"Couldn't parse variable from C function.");

	atomiddata=PyArray_GETCONTIGUOUS(atomiddata);
	int *atomid = (int *) PyArray_DATA(atomiddata);
	occdata=PyArray_GETCONTIGUOUS(occdata);
	float *occ = (float *) PyArray_DATA(occdata);
	bfdata=PyArray_GETCONTIGUOUS(bfdata);
	float *bf = (float *) PyArray_DATA(bfdata);
	cordata=PyArray_GETCONTIGUOUS(cordata);
	float *cor = (float *) PyArray_DATA(cordata);
	map=PyArray_GETCONTIGUOUS(map);
	float *data = (float *) PyArray_DATA(map);


	int atomnumber = PyArray_DIMS(atomiddata)[0];

	int i=0,j=0,k=0,ii=0,jj=0,kk=0,id=0,id0=0,id1=0;
	int isodd=nsam%2;
	int hnsaml=(nsam/2+1);
	int nsaml=hnsaml*2;
	int hnsam=nsam/2;
	float Fpsize=1.0/(psize*nsam);
	float rmax=psize*nsam/res;

	if(rmax>hnsam) 
	{
		rmax=hnsam;
		printf("* Warning: resolution is larger than Nyquist.(From C code.)\n");
	}

	float *pfft=malloc(sizeof(float)*nsam*nsam*nsaml);
	CMPLXf *pfftc=(CMPLXf *)pfft;

	int step=(int)(nsam/10.0+0.5),count=0,p=0;

	#pragma omp parallel for schedule(dynamic) firstprivate(isodd,hnsaml,hnsam,ii,nsam,id0,j,jj,id1,k,kk,id,rmax,Fpsize,atomid,occ,bf,cor,atomnumber) shared(count,p)
	for(i=-hnsam;i<hnsam+isodd;i++)
	{
		ii=(i+nsam)%nsam;
		id0=ii*nsam;
		if (count%step==0){
			#pragma omp critical
			if(p==0){
				printf(" Progress: %%%.0f\n",count*100.0/nsam);
				p=1;
			}
		}
		for(j=-hnsam;j<hnsam+isodd;j++)
		{
			jj=(j+nsam)%nsam;
			id1=(jj+id0)*hnsaml;
			for(k=0;k<hnsaml;k++)
			{
				kk=k;
				id=kk+id1;
				if(sqrt(i*i+j*j+k*k)>rmax) 
				{
					pfftc[id].x=0.0;
					pfftc[id].y=0.0;
				}
				else pfftc[id]=GetStructureFactor(atomid,occ,bf,cor,atomnumber, k*Fpsize, j*Fpsize, i*Fpsize);
			}
		}
		#pragma omp critical
		{
			count++;
			p=0;
		}
	}
	printf(" Progress: %%100\n");
	
	ifft3d(pfft,nsam,data);

	float scale=nsam*nsam*nsam; 
	for(i=0;i<nsam;i++)
	{
		id=i*nsam;
		for(j=0;j<nsam;j++)
		{
			id1=(j+id)*nsam;
			for(k=0;k<nsam;k++)
			{
				data[id1+k]=data[id1+k]/scale;
			}
		}
	}

	free(pfft);

	return Py_None;
}

static PyMethodDef Cmrcmodel_p_methods[] = {

    {"Cpdb2mrc",  (PyCFunction)Cpdb2mrc,
     METH_VARARGS | METH_KEYWORDS,
     "Read the MRC file header into a header variable.\n"},

    {NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef Cmrcmodel_pmodule = {
        PyModuleDef_HEAD_INIT,
        "Cmrcmodel_p",
        "MRC model tools.",
        -1,
        Cmrcmodel_p_methods,
};
PyMODINIT_FUNC PyInit_Cmrcmodel_p(void) {
    import_array();
    return PyModule_Create(&Cmrcmodel_pmodule);
}
#else
PyMODINIT_FUNC initCmrcmodel_p(void) {

    Py_InitModule3("Cmrcmodel_p", Cmrcmodel_p_methods,
        "MRC model tools.");

    import_array();
}
#endif