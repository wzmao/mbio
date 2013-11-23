#include "Python.h"
#include "numpy/arrayobject.h"
#include <stdlib.h>
#include <time.h>

#define NUMCHARS 27
const int twenty[20] = {1, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13,
                        14, 16, 17, 18, 19, 20, 22, 23, 25};
const int unambiguous[23] = {0, 1, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14,
                             15, 16, 17, 18, 19, 20, 21, 22, 23, 25};

static void zeroJoint(double **joint) {
    /* Fill NUMCHARSxNUMCHARS joint array with zeros. */
    int k, l;
    double *jrow;
    for (k = 0; k < NUMCHARS; k++) {
        jrow = joint[k];
        for (l = 0; l < NUMCHARS; l++)
            jrow[l] = 0;
    }
}

static void sortJoint(double **joint) {

    /* Sort probability of ambiguous amino acids. */

    int k, l, t;
    double *jrow, jp, *krow;

    /* X */
    jrow = joint[24];
    /* XX */
    jp = jrow[24];
    if (jp > 0) {
        jp = jp / 400;
        for (k = 0; k < 20; k++) {
            t = twenty[k];
            krow = joint[t];
            for (l = 0; l < 20; l++)
                joint[t][twenty[l]] += jp;
        }
        jrow[24] = 0;
    }
    /* XB */
    jp = jrow[2];
    if (jp > 0) {
        jp = jp / 40;
        for (k = 0; k < 20; k++) {
            krow = joint[twenty[k]];
            krow[4] += jp;
            krow[14] += jp;
        }
        jrow[2] = 0;
    }
    /* XJ */
    jp = jrow[10];
    if (jp > 0) {
        jp = jp / 40;
        for (k = 0; k < 20; k++) {
            krow = joint[twenty[k]];
            krow[9] += jp;
            krow[12] += jp;
        }
        jrow[10] = 0;
    }
    /* XZ */
    jp = jrow[26];
    if (jp > 0) {
        jp = jp / 40;
        for (k = 0; k < 20; k++) {
            krow = joint[twenty[k]];
            krow[5] += jp;
            krow[17] += jp;
        }
        jrow[26] = 0;
    }

    /* B */
    jrow = joint[2];
    /* BB */
    jp = jrow[2];
    if (jp > 0) {
        jp = jp / 4;
        joint[4][4] += jp;
        joint[4][14] += jp;
        joint[14][4] += jp;
        joint[14][14] += jp;
        jrow[2] = 0;
    }
    /* BX */
    jp = jrow[24];
    if (jp > 0) {
        jp = jp / 40;
        for (k = 0; k < 20; k++) {
            t = twenty[k];
            joint[4][t] += jp;
            joint[14][t] += jp;
        }
        jrow[24] = 0;
    }
    /* BJ */
    jp = jrow[10];
    if (jp > 0) {
        jp = jp / 4;
        joint[4][9] += jp;
        joint[4][12] += jp;
        joint[14][9] += jp;
        joint[14][12] += jp;
        jrow[10] = 0;
    }
    /* BZ */
    jp = jrow[26];
    if (jp > 0) {
        jp = jp / 4;
        joint[4][5] += jp;
        joint[4][17] += jp;
        joint[14][5] += jp;
        joint[14][17] += jp;
        jrow[26] = 0;
    }

    /* Z */
    jrow = joint[26];
    /* ZZ */
    jp = jrow[26];
    if (jp > 0) {
        jp = jp / 4;
        joint[5][5] += jp;
        joint[5][17] += jp;
        joint[17][5] += jp;
        joint[17][17] += jp;
        jrow[26] = 0;
    }
    /* ZX */
    jp = jrow[24];
    if (jp > 0) {
        jp = jp / 40;
        for (k = 0; k < 20; k++) {
            t = twenty[k];
            joint[5][t] += jp;
            joint[17][t] += jp;
        }
        jrow[24] = 0;
    }
    /* ZJ */
    jp = jrow[10];
    if (jp > 0) {
        jp = jp / 4;
        joint[5][9] += jp;
        joint[5][12] += jp;
        joint[17][9] += jp;
        joint[17][12] += jp;
        jrow[10] = 0;
    }
    /* ZB */
    jp = jrow[2];
    if (jp > 0) {
        jp = jp / 4;
        joint[5][4] += jp;
        joint[5][14] += jp;
        joint[17][4] += jp;
        joint[17][14] += jp;
        jrow[2] = 0;
    }

    /* J */
    jrow = joint[10];
    /* JJ */
    jp = jrow[10];
    if (jp > 0) {
        jp = jp / 4;
        joint[9][9] += jp;
        joint[9][12] += jp;
        joint[12][9] += jp;
        joint[12][12] += jp;
        joint[10][10] = 0;
    }
    /* JX */
    jp = jrow[24];
    if (jp > 0) {
        jp = jp / 40;
        for (k = 0; k < 20; k++) {
            t = twenty[k];
            joint[9][t] += jp;
            joint[12][t] += jp;
        }
        jrow[24] = 0;
    }
    /* JB */
    jp = jrow[2];
    if (jp > 0) {
        jp = jp / 4;
        joint[9][4] += jp;
        joint[9][14] += jp;
        joint[12][4] += jp;
        joint[12][14] += jp;
        jrow[2] = 0;
    }
    /* BZ */
    jp = jrow[26];
    if (jp > 0) {
        jp = jp / 4;
        joint[9][5] += jp;
        joint[9][17] += jp;
        joint[12][5] += jp;
        joint[12][17] += jp;
        jrow[26] = 0;
    }

    /*for (k = 0; k < NUMCHARS; k++) {*/
    for (t = 0; t < 23; t++) {
        k = unambiguous[t];
        jrow = joint[k];

        /* B */
        jp = jrow[2];
        if (jp > 0) {
            jp = jp / 2;
            jrow[4] += jp;
            jrow[14] += jp;
            jrow[2] = 0;
        }
        jp = joint[2][k];
        if (jp > 0) {
            jp  = jp / 2;
            joint[4][k] += jp;
            joint[14][k] += jp;
            joint[2][k] = 0;
        }

        /* J */
        jp = jrow[10];
        if (jp > 0) {
            jp  = jp / 2;
            jrow[9] += jp;
            jrow[12] += jp;
            jrow[10] = 0;
        }
        jp = joint[10][k];
        if (jp > 0) {
            jp = jp / 2;
            joint[9][k] += jp;
            joint[12][k] += jp;
            joint[10][k] = 0;
        }

        /* Z */
        jp = jrow[26];
        if (jp > 0) {
            jp = jp / 2;
            jrow[5] += jp;
            jrow[17] += jp;
            jrow[26] = 0;
        }
        jp = joint[26][k];
        if (jp > 0) {
            jp = jp / 2;
            joint[5][k] += jp;
            joint[17][k] += jp;
            joint[26][k] = 0;
        }

        /* X */
        jp = jrow[24];
        if (jp > 0) {
            jp = jp / 20.;
            for (l = 0; l < 20; l++)
                jrow[twenty[l]] += jp;
            jrow[24] = 0;
        }
        jp = joint[24][k];
        if (jp > 0) {
            jp = jp / 20.;
            for (l = 0; l < 20; l++)
                joint[twenty[l]][k] += jp;
            joint[24][k] = 0;
        }
    }

}
static double calcMI(double **joint, double **probs, long i, long j) {
    /* Calculate mutual information for a pair of columns in MSA. */
    int k, l;
    double *jrow, *iprb = probs[i], *jprb = probs[j], jp, mi = 0., inside;
    for (k = 0; k < NUMCHARS; k++) {
        jrow = joint[k];
        for (l = 0; l < NUMCHARS; l++) {
            jp = jrow[l];
            if (jp > 0) {
                inside = jp / iprb[k] / jprb[l];
                if (inside != 1)
                    mi += jp * log(inside);
            }
        }
    }
    return mi;
}

static PyObject *shufflemi(PyObject *self, PyObject *args, PyObject *kwargs) {

    PyArrayObject *msa, *pinfo;
    long ambiguity = 1, times = 10000;

    static char *kwlist[] = {"msa", "p",
                             "ambiguity", "times", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|ll", kwlist,
                                     &msa, &pinfo,
                                     &ambiguity, &times))
        return NULL;

    /* make sure to have a contiguous and well-behaved array */
    msa = PyArray_GETCONTIGUOUS(msa);

    /* check dimensions */
    long number = msa->dimensions[0], length = msa->dimensions[1];

    /* get pointers to data */
    char *seq = (char *) PyArray_DATA(msa); /*size: number x length */
    double *p = (double *) PyArray_DATA(pinfo);


    long i, j;
    /* allocate memory */
    double *mut = malloc(length*length*sizeof(double));
    if (!mut)
        return PyErr_NoMemory();
    unsigned char *temp = malloc(number * sizeof(unsigned char));
    if (!temp)
        return PyErr_NoMemory();
    unsigned char *iseq = malloc(number * sizeof(unsigned char));
    if (!iseq){
        free(mut);
        free(temp);
        return PyErr_NoMemory();
    }

    /* hold transpose of the sorted character array */
    unsigned char **trans = malloc(length * sizeof(unsigned char *));
    if (!trans) {
        free(mut);
        free(temp);
        free(iseq);
        return PyErr_NoMemory();
    }

    /* allocate rows that will store columns of MSA */
    trans[0] = iseq;
    for (i = 1; i < length; i++) {
        trans[i] = malloc(number * sizeof(unsigned char));
        if (!trans[i]) {
            for (j = 1; j < i; j++)
                free(trans[j]);
            free(trans);
            free(iseq);
            free(temp);
            free(mut);
            return PyErr_NoMemory();
        }
    }
    unsigned char *jseq = iseq; /* so that we don't get uninitialized warning*/

    /* length*27, a row for each column in the MSA */
    double **probs = malloc(length * sizeof(double *)), *prow;
    if (!probs) {
        for (j = 1; j < length; j++)
            free(trans[j]);
        free(trans);
        free(iseq);
        free(mut);
        free(temp);
        return PyErr_NoMemory();
    }

    /* 27x27, alphabet characters and a gap*/
    double **joint = malloc(NUMCHARS * sizeof(double *)), *jrow;
    if (!joint) {
        for (j = 1; j < length; j++)
            free(trans[j]);
        free(trans);
        free(iseq);
        free(probs);
        free(mut);
        free(temp);
        return PyErr_NoMemory();
    }

    for (i = 0; i < length; i++) {
        prow = malloc(NUMCHARS * sizeof(double));
        if (!prow) {
            for (j = 0; j < i; j++)
                free(probs[j]);
            free(probs);
            free(joint);
            for (j = 1; j < length; j++)
                free(trans[j]);
            free(trans);
            free(iseq);
            free(mut);
            free(temp);
            return PyErr_NoMemory();
        }
        probs[i] = prow;
        for (j = 0; j < NUMCHARS; j++)
            prow[j] = 0;
    }

    for (i = 0; i < NUMCHARS; i++)  {
        joint[i] = malloc(NUMCHARS * sizeof(double));
        if (!joint[i]) {
            for (j = 0; j < i; j++)
                free(joint[j]);
            free(joint);
            for (j = 0; j < length; j++)
                free(probs[j]);
            free(probs);
            for (j = 1; j < length; j++)
                free(trans[j]);
            free(trans);
            free(iseq);
            free(mut);
            free(temp);
            return PyErr_NoMemory();
        }
    }

    unsigned char a;
    long k, l, offset, e;
    double p_incr = 1. / number, prb = 0;
    prow = probs[0];

    /* START mutinfo calculation */
    /* calculate first row of MI matrix and all column probabilities */
    for (i=0;i<length*length;i++){
        p[i]=0;
        mut[i]=0.;
    }
    for (i = 0; i < length; i++) {
        jrow = probs[i];
        jseq = trans[i];
        for (j = 0; j < number; j++) {
            offset = j * length;
            a = (unsigned char) seq[offset + i];
            if (a > 90)
                a -= 96;
            else
                a -= 64;
            if (a < 1 || a > 26)
                a = 0; /* gap character */
            jseq[j] = a;
            jrow[a] += p_incr;
        }

        if (ambiguity) {
            prb = jrow[2];
            if (prb > 0) { /* B -> D, N  */
                prb = prb / 2.;
                jrow[4] += prb;
                jrow[14] += prb;
                jrow[2] = 0.;
            }
            prb = jrow[10];
            if (prb > 0) { /* J -> I, L  */
                prb = prb / 2.;
                jrow[9] += prb;
                jrow[12] += prb;
                jrow[10] = 0.;
            }
            prb = jrow[26];
            if (prb > 0) { /* Z -> E, Q  */
                prb = prb / 2.;
                jrow[5] += prb;
                jrow[17] += prb;
                jrow[26] = 0.;
            }
            if (jrow[24] > 0) { /* X -> 20 AA */
                prb = jrow[24] / 20.;
                for (l = 0; l < 20; l++)
                    jrow[twenty[l]] += prb;
                jrow[24] = 0.;
            }
        }
    }

    /* calculate rest of MI matrix */
    long ioffset;
    for (i = 0; i < length; i++) {
        ioffset = i * length;
        mut[ioffset+i] = 0.;
        iseq = trans[i];
        for (j = i + 1; j < length; j++) {
            zeroJoint(joint);
            jseq = trans[j];
            for (k = 0; k < number; k++)
                joint[iseq[k]][jseq[k]] += p_incr;
            if (ambiguity)
                sortJoint(joint);
            mut[ioffset + j] = mut[i + length * j] =
                calcMI(joint, probs, i, j);
        }
    }

    double newmi=0.0,tiny=1e-15;
    long id;
    unsigned char sw;
    srand((unsigned)time(NULL));
    for (i=0;i<length;i++){
        iseq=trans[i];
        ioffset=i*length;
        for (j=0;j<number;j++)
            temp[j]=iseq[j];
        for (e=0;e<times;e++){
            /*shuffle*/
            for (j=0;j<number;j++){
                id=(rand()%(number-j))+j;
                sw=temp[id];
                temp[id]=temp[j];
                temp[j]=sw;
            }
            for (j=0;j<length;j++){
                if (j==i)
                    continue;
                /*build joint*/
                zeroJoint(joint);
                jseq=trans[j];
                for (k = 0; k < number; k++)
                    joint[temp[k]][jseq[k]] += p_incr;
                if (ambiguity)
                    sortJoint(joint);
                newmi = calcMI(joint, probs, i, j);
                if (newmi >= mut[ioffset +j]-tiny){
                    p[ioffset +j]++;
                    p[j*length+i]++;
                }
            }
        }
    }
    for (i=0;i<length*length;i++)
        p[i]/=times*2;

    /* free memory */
    for (i = 0; i < length; i++){
        free(probs[i]);
    }
    free(probs);
    for (i = 0; i < NUMCHARS; i++){
        free(joint[i]);
    }
    free(joint);
    for (j = 1; j < length; j++)
        free(trans[j]);
    free(trans);
    free(mut);
    free(temp);
    return Py_BuildValue("O", pinfo);
}

static double calcOMES(double **joint, double **probs, long i, long j, int n) {

    /* Calculate OMES for a pair of columns in MSA. */
    int k, l;
    double *jrow, *iprb = probs[i], *jprb = probs[j], jp, omes = 0, inside;
    for (k = 0; k < NUMCHARS; k++) {
        jrow = joint[k];
        for (l = 0; l < NUMCHARS; l++) {
            jp = jrow[l];
            inside = iprb[k] * jprb[l];
            if (inside != 0)
                omes += n * (jp - inside) * (jp - inside) / inside;
        }
    }
    return omes;
}

static PyObject *shuffleomes(PyObject *self, PyObject *args, PyObject *kwargs) {

    PyArrayObject *msa, *pinfo;
    long ambiguity = 1, times = 10000;

    static char *kwlist[] = {"msa", "p",
                             "ambiguity", "times", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|ll", kwlist,
                                     &msa, &pinfo,
                                     &ambiguity, &times))
        return NULL;

    /* make sure to have a contiguous and well-behaved array */
    msa = PyArray_GETCONTIGUOUS(msa);

    /* check dimensions */
    long number = msa->dimensions[0], length = msa->dimensions[1];

    /* get pointers to data */
    char *seq = (char *) PyArray_DATA(msa); /*size: number x length */
    double *p = (double *) PyArray_DATA(pinfo);


    long i, j;
    /* allocate memory */
    double *omes = malloc(length*length*sizeof(double));
    if (!omes)
        return PyErr_NoMemory();
    unsigned char *temp = malloc(number * sizeof(unsigned char));
    if (!temp)
        return PyErr_NoMemory();
    unsigned char *iseq = malloc(number * sizeof(unsigned char));
    if (!iseq){
        free(omes);
        free(temp);
        return PyErr_NoMemory();
    }

    /* hold transpose of the sorted character array */
    unsigned char **trans = malloc(length * sizeof(unsigned char *));
    if (!trans) {
        free(omes);
        free(temp);
        free(iseq);
        return PyErr_NoMemory();
    }

    /* allocate rows that will store columns of MSA */
    trans[0] = iseq;
    for (i = 1; i < length; i++) {
        trans[i] = malloc(number * sizeof(unsigned char));
        if (!trans[i]) {
            for (j = 1; j < i; j++)
                free(trans[j]);
            free(trans);
            free(iseq);
            free(temp);
            free(omes);
            return PyErr_NoMemory();
        }
    }
    unsigned char *jseq = iseq; /* so that we don't get uninitialized warning*/

    /* length*27, a row for each column in the MSA */
    double **probs = malloc(length * sizeof(double *)), *prow;
    if (!probs) {
        for (j = 1; j < length; j++)
            free(trans[j]);
        free(trans);
        free(iseq);
        free(omes);
        free(temp);
        return PyErr_NoMemory();
    }

    /* 27x27, alphabet characters and a gap*/
    double **joint = malloc(NUMCHARS * sizeof(double *)), *jrow;
    if (!joint) {
        for (j = 1; j < length; j++)
            free(trans[j]);
        free(trans);
        free(iseq);
        free(probs);
        free(omes);
        free(temp);
        return PyErr_NoMemory();
    }

    for (i = 0; i < length; i++) {
        prow = malloc(NUMCHARS * sizeof(double));
        if (!prow) {
            for (j = 0; j < i; j++)
                free(probs[j]);
            free(probs);
            free(joint);
            for (j = 1; j < length; j++)
                free(trans[j]);
            free(trans);
            free(iseq);
            free(omes);
            free(temp);
            return PyErr_NoMemory();
        }
        probs[i] = prow;
        for (j = 0; j < NUMCHARS; j++)
            prow[j] = 0;
    }

    for (i = 0; i < NUMCHARS; i++)  {
        joint[i] = malloc(NUMCHARS * sizeof(double));
        if (!joint[i]) {
            for (j = 0; j < i; j++)
                free(joint[j]);
            free(joint);
            for (j = 0; j < length; j++)
                free(probs[j]);
            free(probs);
            for (j = 1; j < length; j++)
                free(trans[j]);
            free(trans);
            free(iseq);
            free(omes);
            free(temp);
            return PyErr_NoMemory();
        }
    }

    unsigned char a;
    long k, l, offset, e;
    double p_incr = 1. / number, prb = 0;
    prow = probs[0];

    /* START omesinfo calculation */
    /* calculate first row of OMES matrix and all column probabilities */
    for (i=0;i<length*length;i++){
        p[i]=0;
        omes[i]=0.;
    }
    for (i = 0; i < length; i++) {
        jrow = probs[i];
        jseq = trans[i];
        for (j = 0; j < number; j++) {
            offset = j * length;
            a = (unsigned char) seq[offset + i];
            if (a > 90)
                a -= 96;
            else
                a -= 64;
            if (a < 1 || a > 26)
                a = 0; /* gap character */
            jseq[j] = a;
            jrow[a] += p_incr;
        }

        if (ambiguity) {
            prb = jrow[2];
            if (prb > 0) { /* B -> D, N  */
                prb = prb / 2.;
                jrow[4] += prb;
                jrow[14] += prb;
                jrow[2] = 0.;
            }
            prb = jrow[10];
            if (prb > 0) { /* J -> I, L  */
                prb = prb / 2.;
                jrow[9] += prb;
                jrow[12] += prb;
                jrow[10] = 0.;
            }
            prb = jrow[26];
            if (prb > 0) { /* Z -> E, Q  */
                prb = prb / 2.;
                jrow[5] += prb;
                jrow[17] += prb;
                jrow[26] = 0.;
            }
            if (jrow[24] > 0) { /* X -> 20 AA */
                prb = jrow[24] / 20.;
                for (l = 0; l < 20; l++)
                    jrow[twenty[l]] += prb;
                jrow[24] = 0.;
            }
        }
    }

    /* calculate rest of OMES matrix */
    long ioffset;
    for (i = 0; i < length; i++) {
        ioffset = i * length;
        omes[ioffset+i] = 0.;
        iseq = trans[i];
        for (j = i + 1; j < length; j++) {
            zeroJoint(joint);
            jseq = trans[j];
            for (k = 0; k < number; k++)
                joint[iseq[k]][jseq[k]] += p_incr;
            if (ambiguity)
                sortJoint(joint);
            omes[ioffset + j] = omes[i + length * j] =
                calcOMES(joint, probs, i, j, number);
        }
    }

    double newomes=0.0,tiny=1e-15;
    long id;
    unsigned char sw;
    srand((unsigned)time(NULL));
    for (i=0;i<length;i++){
        iseq=trans[i];
        ioffset=i*length;
        for (j=0;j<number;j++)
            temp[j]=iseq[j];
        for (e=0;e<times;e++){
            /*shuffle*/
            for (j=0;j<number;j++){
                id=(rand()%(number-j))+j;
                sw=temp[id];
                temp[id]=temp[j];
                temp[j]=sw;
            }
            for (j=0;j<length;j++){
                if (j==i)
                    continue;
                /*build joint*/
                zeroJoint(joint);
                jseq=trans[j];
                for (k = 0; k < number; k++)
                    joint[temp[k]][jseq[k]] += p_incr;
                if (ambiguity)
                    sortJoint(joint);
                newomes = calcOMES(joint, probs, i, j, number);
                if (newomes >= omes[ioffset +j]-tiny){
                    p[ioffset +j]++;
                    p[j*length+i]++;
                }
            }
        }
    }
    for (i=0;i<length*length;i++)
        p[i]/=times*2;

    /* free memory */
    for (i = 0; i < length; i++){
        free(probs[i]);
    }
    free(probs);
    for (i = 0; i < NUMCHARS; i++){
        free(joint[i]);
    }
    free(joint);
    for (j = 1; j < length; j++)
        free(trans[j]);
    free(trans);
    free(omes);
    free(temp);
    return Py_BuildValue("O", pinfo);
}

static PyObject *shufflemip(PyObject *self, PyObject *args, PyObject *kwargs) {

    PyArrayObject *msa, *pinfo;
    long ambiguity = 1, times = 10000;

    static char *kwlist[] = {"msa", "p",
                             "ambiguity", "times", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|ll", kwlist,
                                     &msa, &pinfo,
                                     &ambiguity, &times))
        return NULL;

    /* make sure to have a contiguous and well-behaved array */
    msa = PyArray_GETCONTIGUOUS(msa);

    /* check dimensions */
    long number = msa->dimensions[0], length = msa->dimensions[1];

    /* get pointers to data */
    char *seq = (char *) PyArray_DATA(msa); /*size: number x length */
    double *p = (double *) PyArray_DATA(pinfo);


    long i, j;
    /* allocate memory */
    double *mi = malloc(length*length*sizeof(double));
    if (!mi)
        return PyErr_NoMemory();
    double *mip = malloc(length*length*sizeof(double));
    if (!mip){
        free(mi);
    }
    unsigned char *temp = malloc(number * sizeof(unsigned char));
    if (!temp)
        return PyErr_NoMemory();
    unsigned char *iseq = malloc(number * sizeof(unsigned char));
    if (!iseq){
        free(mi);
        free(mip);
        free(temp);
        return PyErr_NoMemory();
    }

    /* hold transpose of the sorted character array */
    unsigned char **trans = malloc(length * sizeof(unsigned char *));
    if (!trans) {
        free(mi);
        free(mip);
        free(temp);
        free(iseq);
        return PyErr_NoMemory();
    }

    /* allocate rows that will store columns of MSA */
    trans[0] = iseq;
    for (i = 1; i < length; i++) {
        trans[i] = malloc(number * sizeof(unsigned char));
        if (!trans[i]) {
            for (j = 1; j < i; j++)
                free(trans[j]);
            free(trans);
            free(iseq);
            free(temp);
            free(mi);
            free(mip);
            return PyErr_NoMemory();
        }
    }
    unsigned char *jseq = iseq; /* so that we don't get uninitialized warning*/

    /* length*27, a row for each column in the MSA */
    double **probs = malloc(length * sizeof(double *)), *prow;
    if (!probs) {
        for (j = 1; j < length; j++)
            free(trans[j]);
        free(trans);
        free(iseq);
        free(mi);
        free(temp);
        free(mip);
        return PyErr_NoMemory();
    }

    /* 27x27, alphabet characters and a gap*/
    double **joint = malloc(NUMCHARS * sizeof(double *)), *jrow;
    if (!joint) {
        for (j = 1; j < length; j++)
            free(trans[j]);
        free(trans);
        free(iseq);
        free(probs);
        free(mi);
        free(temp);
        free(mip);
        return PyErr_NoMemory();
    }

    for (i = 0; i < length; i++) {
        prow = malloc(NUMCHARS * sizeof(double));
        if (!prow) {
            for (j = 0; j < i; j++)
                free(probs[j]);
            free(probs);
            free(joint);
            for (j = 1; j < length; j++)
                free(trans[j]);
            free(trans);
            free(iseq);
            free(mi);
            free(temp);
            free(mip);
            return PyErr_NoMemory();
        }
        probs[i] = prow;
        for (j = 0; j < NUMCHARS; j++)
            prow[j] = 0;
    }

    for (i = 0; i < NUMCHARS; i++)  {
        joint[i] = malloc(NUMCHARS * sizeof(double));
        if (!joint[i]) {
            for (j = 0; j < i; j++)
                free(joint[j]);
            free(joint);
            for (j = 0; j < length; j++)
                free(probs[j]);
            free(probs);
            for (j = 1; j < length; j++)
                free(trans[j]);
            free(trans);
            free(iseq);
            free(mi);
            free(temp);
            free(mip);
            return PyErr_NoMemory();
        }
    }

    unsigned char a;
    long k, l, offset, e;
    double p_incr = 1. / number, prb = 0;
    prow = probs[0];

    /* START miinfo calculation */
    /* calculate first row of MI matrix and all column probabilities */
    for (i=0;i<length*length;i++){
        p[i]=0;
        mi[i]=0.;
        mip[i]=0.;
    }
    for (i = 0; i < length; i++) {
        jrow = probs[i];
        jseq = trans[i];
        for (j = 0; j < number; j++) {
            offset = j * length;
            a = (unsigned char) seq[offset + i];
            if (a > 90)
                a -= 96;
            else
                a -= 64;
            if (a < 1 || a > 26)
                a = 0; /* gap character */
            jseq[j] = a;
            jrow[a] += p_incr;
        }

        if (ambiguity) {
            prb = jrow[2];
            if (prb > 0) { /* B -> D, N  */
                prb = prb / 2.;
                jrow[4] += prb;
                jrow[14] += prb;
                jrow[2] = 0.;
            }
            prb = jrow[10];
            if (prb > 0) { /* J -> I, L  */
                prb = prb / 2.;
                jrow[9] += prb;
                jrow[12] += prb;
                jrow[10] = 0.;
            }
            prb = jrow[26];
            if (prb > 0) { /* Z -> E, Q  */
                prb = prb / 2.;
                jrow[5] += prb;
                jrow[17] += prb;
                jrow[26] = 0.;
            }
            if (jrow[24] > 0) { /* X -> 20 AA */
                prb = jrow[24] / 20.;
                for (l = 0; l < 20; l++)
                    jrow[twenty[l]] += prb;
                jrow[24] = 0.;
            }
        }
    }

    /* calculate rest of MI matrix */
    long ioffset;
    for (i = 0; i < length; i++) {
        ioffset = i * length;
        mi[ioffset+i] = 0.;
        iseq = trans[i];
        for (j = i + 1; j < length; j++) {
            zeroJoint(joint);
            jseq = trans[j];
            for (k = 0; k < number; k++)
                joint[iseq[k]][jseq[k]] += p_incr;
            if (ambiguity)
                sortJoint(joint);
            mi[ioffset + j] = mi[i + length * j] =
                calcMI(joint, probs, i, j);
        }
    }
    double sum[length],allsum=0.,newsum[length],newallsum=0.;
    for (i=0;i<length;i++){
        sum[i]=0.;
        newsum[i]=0.;
    }
    for (i=0;i<length;i++){
        for (j=i+1;j<length;j++){
            allsum+=2*mi[i*length+j];
            sum[i]+=mi[i*length+j];
            sum[j]+=mi[i*length+j];
        }
    }
    double factor=1.*length/(length-1);
    for(i=0;i<length;i++)
        for(j=0;j<length;j++){
            mip[i*length+j]=mip[j*length+i]=mi[i*length+j]-sum[i]*sum[j]/allsum*factor;
        }
    double newmi[length],tiny=1e-15;
    long id;
    unsigned char sw;
    srand((unsigned)time(NULL));
    for (i=0;i<length;i++){
        iseq=trans[i];
        ioffset=i*length;
        newmi[i]=0.;
        for (j=0;j<number;j++)
            temp[j]=iseq[j];
        for (e=0;e<times;e++){
            /*shuffle*/
            for (j=0;j<number;j++){
                id=(rand()%(number-j))+j;
                sw=temp[id];
                temp[id]=temp[j];
                temp[j]=sw;
            }
            for (j=0;j<length;j++)
                newsum[j]=sum[j];
            newallsum=allsum;
            for (j=0;j<length;j++){
                if (j==i)
                    continue;
                /*build joint*/
                zeroJoint(joint);
                jseq=trans[j];
                for (k = 0; k < number; k++)
                    joint[temp[k]][jseq[k]] += p_incr;
                if (ambiguity)
                    sortJoint(joint);
                newmi[j] = calcMI(joint, probs, i, j);
                newsum[j]-=mi[ioffset+j]+newmi[j];
                newsum[i]-=mi[ioffset+j]+newmi[j];
                newallsum-=2*(mi[ioffset+j]+newmi[j]);
            }
            for (j=0;j<length;j++){
                if (j==i)
                    continue;
                if ((newmi[j]-newsum[i]*newsum[j]/newallsum*factor)>=mip[ioffset+j]-tiny){
                    p[ioffset +j]++;
                    p[j*length+i]++;
                }
            }
        }
    }
    for (i=0;i<length*length;i++)
        p[i]/=times*2;

    /* free memory */
    for (i = 0; i < length; i++){
        free(probs[i]);
    }
    free(probs);
    for (i = 0; i < NUMCHARS; i++){
        free(joint[i]);
    }
    free(joint);
    for (j = 1; j < length; j++)
        free(trans[j]);
    free(trans);
    free(mi);
    free(temp);
    free(mip);
    return Py_BuildValue("O", pinfo);
}

static PyMethodDef c_shuffle_methods[] = {

    {"shufflemi",  (PyCFunction)shufflemi,
     METH_VARARGS | METH_KEYWORDS,
     "Return shuffled result for MI."},

    {"shuffleomes",  (PyCFunction)shuffleomes,
     METH_VARARGS | METH_KEYWORDS,
     "Return shuffled result for OMES."},

    {"shufflemip",  (PyCFunction)shufflemip,
     METH_VARARGS | METH_KEYWORDS,
     "Return shuffled result for MIp."},

    {NULL, NULL, 0, NULL}
};


#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef c_shuffle = {
        PyModuleDef_HEAD_INIT,
        "c_shuffle",
        "MSA shuffle tools.",
        -1,
        c_shuffle_methods,
};
PyMODINIT_FUNC PyInit_c_shuffle(void) {
    import_array();
    return PyModule_Create(&c_shuffle);
}
#else
PyMODINIT_FUNC initc_shuffle(void) {

    Py_InitModule3("c_shuffle", c_shuffle_methods,
        "MSA shuffle tools.");

    import_array();
}
#endif
