/*=========================================================================

  Program:   Utilities
  Module:    netcdfio.cpp

  This module handles the reading and eventual writing of NetCDF files
  for the Simulation Paraview Suite

  Copyright (c) 2018 Arnau Miro, Simulation
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef NETCDFIO_H
#define NETCDFIO_H

#include <string>
#include <vector>

#include "V3.h"
#include "field.h"

#define NETCDF_ERR 0
#define NETCDF_OK 1
#define MISSING_VALUE 1.e20

namespace OGS::NetCDF {

/* Function to read a single variable */
template <typename T>
int NetCDFGetVar(int fid, const char *varname, unsigned int n, T *out,
                 bool m2zero);


/* Read NetCDF Routines */
template <typename T>
int readNetCDF(const char *fname, const char *varname, unsigned int n, T *out,
               bool m2zero);

template <typename T>
int readNetCDF(const char *fname, const char *varname, std::vector<T> &out,
               bool m2zero = true) {
  return readNetCDF(fname, varname, (int)(out.size()), &out[0], m2zero);
}

template <typename T>
int readNetCDF(const char *fname, const char *varname, OGS::field::Field<T> &f,
               bool m2zero = true) {
  return readNetCDF(fname, varname, f.get_n(), f.data(), m2zero);
}

template <typename T>
int readNetCDF(const char *fname, std::vector<std::string> var_names,
               OGS::field::Field<T> &f);


/* Write NetCDF Routines */
int writeNetCDF(const char *fname, const char *varname, int dims[], double *lon,
                double *lat, double *depth, double *data);

int writeNetCDF(const char *fname, const char *varname, int dims[], float *lon,
                float *lat, float *depth, float *data);

int writeNetCDF(const char *fname, std::string *varname, unsigned int nvars,
                int dims[], double *lon, double *lat, double *depth,
                double **data);

int writeNetCDF(const char *fname, std::string *varname, unsigned int nvars,
                int dims[], float *lon, float *lat, float *depth, float **data);

int writeNetCDF(const char *fname, const char *varname, int dims[], double *lon,
                double *lat, double *depth, OGS::field::Field<double> &f);

int writeNetCDF(const char *fname, const char *varname, int dims[], float *lon,
                float *lat, float *depth, OGS::field::Field<float> &f);

int writeNetCDFProfile(const char *fname, const char *varname, int dims,
                       double *lon, double *lat, double *depth, double *data);

int writeNetCDFProfile(const char *fname, const char *varname, int dims,
                       float *lon, float *lat, float *depth, float *data);

int writeNetCDFProfile(const char *fname, std::string *varname, int nvars,
                       int dims, double *lon, double *lat, double *depth,
                       double **data);

int writeNetCDFProfile(const char *fname, std::string *varname, int nvars,
                       int dims, float *lon, float *lat, float *depth,
                       float **data);

int writeNetCDFProfile(const char *fname, const char *varname, int dims,
                       double *lon, double *lat, double *depth,
                       OGS::field::Field<double> &f);

int writeNetCDFProfile(const char *fname, const char *varname, int dims,
                       float *lon, float *lat, float *depth,
                       OGS::field::Field<float> &f);

}  // namespace OGS::NetCDF
#endif
