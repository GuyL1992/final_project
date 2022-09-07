#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "spkmeans.h"



static PyObject* getDataPoints_capi(PyObject *self, PyObject *args);
static void convert_cMatrix_to_pyMAtrix(PyObject* pyMatrix, double** cMatrix, int n,int d);
static void convert_data_to_c(PyObject* data_points_py, double** data_points, int n, int d);
static void convert_point_to_c(PyObject* point_py, double* point, int d);
static PyObject* wam_capi(PyObject *self, PyObject *args);

static PyObject* getDataPoints_capi(PyObject *self, PyObject *args){   
   
    PyObject* source_dataPy;
    PyObject* PyNewDataPoints;
    int n,d,k;
    double** sourceDataForC, **newDataPoints;

    if(!PyArg_ParseTuple(args,"iiiO",&n,&d,&k,&source_dataPy)){
        return NULL;
    }
    PyNewDataPoints = PyList_New(n);
    sourceDataForC = allocationMatrix(n,d);
    convert_data_to_c(source_dataPy,sourceDataForC,n,d);
    newDataPoints = getDataPoints(sourceDataForC,n,d,k);
    convert_cMatrix_to_pyMAtrix(PyNewDataPoints,newDataPoints,n,3); // adress of k!!!


    // convert_data_to_c(data_points_py,data_points,n,d);
    // newDataPoints = getNewDataPointsDimK(data_points,n,d,&k);
    // for (i=0;i<n;i++){
    //     currentLine = PyList_New(k);
    //     for (j=0;j<k;j++){
    //         item = PyFloat_FromDouble(newDataPoints[i][j]);
    //         PyList_SetItem(currentLine,j,item);
    //     }
    //     PyList_SetItem(PyNewDataPoints,i,currentLine);
    //     free(newDataPoints[i]);
    // }
    // for (i=0;i<n;i++){
    //     free(data_points[i]);
    // }
    // free(newDataPoints);
    // free(data_points);
    return PyNewDataPoints;

}

static PyObject* kmeans_pp_capi(PyObject *self, PyObject *args)
{
    PyObject* data_points_py;
    PyObject* centroids_py;
    int n;
    int d;
    int k;
    int i;
    int max_iter;
    double** data_points;
    double** centroids;

    if(!PyArg_ParseTuple(args,"iiiiOO",&n,&d,&k,&max_iter,&data_points_py,&centroids_py)){
        return NULL;
    }

    data_points = allocationMatrix(n,k);
    centroids = allocationMatrix(k,k);
    convert_data_to_c(data_points_py,data_points,n,d);
    convert_data_to_c(centroids_py,centroids,k,k);
    calculateKCentroids(n,d,k,max_iter,data_points,centroids);
    
    for (i=0;i<k;i++){
        free(centroids[i]);
    }
    for (i=0;i<n;i++){
        free(data_points[i]);
    }
    free(centroids);
    free(data_points);
    Py_RETURN_NONE;

}

static PyObject* wam_capi(PyObject *self, PyObject *args){
    PyObject* source_dataPy;
    PyObject* PyNewDataPoints;
    int n,d;
    double** sourceDataForC;

    if(!PyArg_ParseTuple(args,"iiO",&n,&d,&source_dataPy)){
        return NULL;
    }
    PyNewDataPoints = PyList_New(n);
    sourceDataForC = allocationMatrix(n,d);
    convert_data_to_c(source_dataPy,sourceDataForC,n,d);
    wamProcess(sourceDataForC, n,d);
    Py_RETURN_NONE;
}

static PyObject* ddg_capi(PyObject *self, PyObject *args){
    PyObject* source_dataPy;
    PyObject* PyNewDataPoints;
    int n,d;
    double** sourceDataForC;

    if(!PyArg_ParseTuple(args,"iiO",&n,&d,&source_dataPy)){
        return NULL;
    }
    PyNewDataPoints = PyList_New(n);
    sourceDataForC = allocationMatrix(n,d);
    convert_data_to_c(source_dataPy,sourceDataForC,n,d);
    ddgProcess(sourceDataForC, n,d);
    freeMatrix(sourceDataForC,n);
    Py_RETURN_NONE;
}

static PyObject* jacobi_capi(PyObject *self, PyObject *args){
    PyObject* source_dataPy;
    PyObject* PyNewDataPoints;
    int n,d;
    double** sourceDataForC;

    if(!PyArg_ParseTuple(args,"iiO",&n,&d,&source_dataPy)){
        return NULL;
    }
    PyNewDataPoints = PyList_New(n);
    sourceDataForC = allocationMatrix(n,d);
    convert_data_to_c(source_dataPy,sourceDataForC,n,d);
    jacobiProcess(sourceDataForC, n,d);
    freeMatrix(sourceDataForC,n);
    Py_RETURN_NONE;
}

static PyObject* lnorm_capi(PyObject *self, PyObject *args){
    PyObject* source_dataPy;
    PyObject* PyNewDataPoints;
    int n,d;
    double** sourceDataForC;

    if(!PyArg_ParseTuple(args,"iiO",&n,&d,&source_dataPy)){
        return NULL;
    }
    PyNewDataPoints = PyList_New(n);
    sourceDataForC = allocationMatrix(n,d);
    convert_data_to_c(source_dataPy,sourceDataForC,n,d);
    lnormProcess(sourceDataForC, n,d);
    freeMatrix(sourceDataForC,n);
    Py_RETURN_NONE;
}

static void convert_data_to_c(PyObject* data_points_py, double** data_points, int n, int d){
    int i;
    int j;
    double current_item;
    PyObject* point_py;
    double* point;
    for (i=0;i<n;i++){
        point_py = PyList_GET_ITEM(data_points_py,i);
        assert(point!=NULL && "An Error Has Occured");
        point = (double*)calloc(d,sizeof(double));
        data_points[i] = (double *)calloc(d,sizeof(double));
        assert(data_points[i]!=NULL && "An Error Has Occured");
        convert_point_to_c(point_py, point, d);
        for(j=0;j<d;j++){
            current_item = point[j];
            data_points[i][j] = current_item;
        }
        free(point);  
    }
}

static void convert_point_to_c(PyObject* point_py, double* point, int d){
    int i;
    double temp;
    PyObject* item;
    for (i=0;i<d;i++){
        item = PyList_GetItem(point_py,i);
        temp = PyFloat_AsDouble(item);
        point[i] = temp;
        
    }
}

static void convert_cMatrix_to_pyMAtrix(PyObject* pyMatrix, double** cMatrix, int n,int d){

    int i;
    int j;
    PyObject* line;

    for(i = 0; i < n; i++){
        line = PyList_New(d);
        for (j = 0; j < d; j++){
            PyList_SetItem(line,j,PyFloat_FromDouble(cMatrix[i][j]));
        }
        PyList_SetItem(pyMatrix,i,line);
    }

    
}

static PyMethodDef _capiMethods[] = {
    {"kmeans_pp", (PyCFunction) kmeans_pp_capi, METH_VARARGS, PyDoc_STR("Get datapoints and centroids and calc kmeans")},
    {"get_new_data", (PyCFunction) getDataPoints_capi, METH_VARARGS, PyDoc_STR("A function to run Jacobi algorithm to find Eignvalues and Eignvectors")},
    {"jacobi", (PyCFunction) jacobi_capi, METH_VARARGS, PyDoc_STR("A function to run Jacobi algorithm to find Eignvalues and Eignvectors")},
    {"wam", (PyCFunction) wam_capi, METH_VARARGS, PyDoc_STR("A function to find WAM for set of data points")},
    {"ddg", (PyCFunction) ddg_capi, METH_VARARGS, PyDoc_STR("A function to find DDG for set of data points")},
    {"lnorm", (PyCFunction) lnorm_capi, METH_VARARGS, PyDoc_STR("A function to find LNORM for set of data points")},
    {NULL,NULL,0,NULL}
};

static struct PyModuleDef _moduledef = {
    PyModuleDef_HEAD_INIT,
    "spkmeans",
    NULL,
    -1,
    _capiMethods
};

PyMODINIT_FUNC
PyInit_spkmeans(void)
{
    PyObject *m;
    m = PyModule_Create(&_moduledef);
    if(!m){
        return NULL;
    }
    return m;
}