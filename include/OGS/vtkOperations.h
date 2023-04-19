/*=========================================================================

  Program:   VTK Operations
  Module:    vtkOperations.hpp

  Wrapper to perform operations with VTK
  that might be needed accross plugins.

  Copyright (c) 2018 Arnau Miro, Simulation
  All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef VTKOPERATIONS_H
#define VTKOPERATIONS_H

#include <vector>

#include "V3.h"
#include "field.h"
#include "vtkDataSet.h"

/* GETVTKCELLPOINTS

        Recovers the cell point coordinates given a VTKDataSet.

*/
OGS::V3::V3v getVTKCellPoints(vtkDataSet *mesh, double fact);

/* GETVTKCELLCENTERS

        Recovers the cell center coordinates given a VTKDataSet.

*/
OGS::V3::V3v getVTKCellCenters(vtkDataSet *mesh, double fact);

/* GETCONNECTEDPOINTS

  Recovers all the points connected to a point in the mesh.
*/
std::vector<int> getConnectedPoints(vtkDataSet *mesh, int pointId);

/* GETCONNECTEDCELLS

  Recovers all the cells connected to a cell in the mesh.
*/
std::vector<int> getConnectedCells(vtkDataSet *mesh, int cellId);

#endif