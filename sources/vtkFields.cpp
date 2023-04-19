/*=========================================================================

  Program:   VTK Fields
  Module:    vtkFields.cpp

  Wrapper to create many kinds of VTK fields
  that might be needed across plugins.

  Copyright (c) 2018 Arnau Miro, Simulation
  All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "OGS/vtkFields.h"

namespace VTK {

/* CREATEVTKSTRF

        Creates a vtk string array given the name
        and data to be converted.

*/
vtkStringArray *createVTKstrf(const char *name, const unsigned int n) {
  // Create string array
  vtkStringArray *vtkArray = vtkStringArray::New();
  vtkArray->SetName(name);
  vtkArray->SetNumberOfTuples(n);

  return vtkArray;
}

vtkStringArray *createVTKstrf(const char *name, const unsigned int n,
                              const char *data) {
  vtkStringArray *vtkArray = createVTKstrf(name, n);

  // Set the value
  for (unsigned int ii = 0; ii < n; ii++) vtkArray->SetValue(ii, data);

  return vtkArray;
}

/* CREATERECTILINEARGRID

        Creates a VTK rectilinear grid given the mesh dimensions and the
        conversion to meters.

*/
void createRectilinearGrid(int nx, int ny, int nz, double *x, double *y,
                           double *z, double scalf, vtkRectilinearGrid *rgrid) {
  // Set dimension arrays
  vtkDoubleArray *vtkx;
  vtkx = createVTKscaf<vtkDoubleArray, double>("x coord", nx, x);
  vtkDoubleArray *vtky;
  vtky = createVTKscaf<vtkDoubleArray, double>("y coord", ny, y);
  vtkDoubleArray *vtkz;
  vtkz = createVTKscaf<vtkDoubleArray, double>("z coord", nz, z);

  // Fix scaling in z
  for (int ii = 0; ii < nz; ii += 1) vtkz->SetTuple1(ii, -scalf * z[ii]);

  // Set rectilinear grid
  rgrid->SetDimensions(nx, ny, nz);
  rgrid->SetXCoordinates(vtkx);
  vtkx->Delete();
  rgrid->SetYCoordinates(vtky);
  vtky->Delete();
  rgrid->SetZCoordinates(vtkz);
  vtkz->Delete();
}

}  // namespace VTK
