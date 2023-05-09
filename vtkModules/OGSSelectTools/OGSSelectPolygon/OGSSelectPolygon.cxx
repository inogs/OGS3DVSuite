/*=========================================================================

  Program:   OGSSelectTools
  Module:    OGSSelectPolygon.cxx

  Copyright (c) 2019 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "OGSSelectPolygon.h"

#include <algorithm>
#include <cstdint>
#include <vector>

#include "vtkCellData.h"
#include "vtkDataArraySelection.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkTypeUInt8Array.h"
#include "vtkUnstructuredGrid.h"

#ifdef PARAVIEW_USE_MPI
#include "vtkMultiProcessController.h"
vtkCxxSetObjectMacro(OGSSelectPolygon, Controller, vtkMultiProcessController);
#endif

vtkStandardNewMacro(OGSSelectPolygon);

//----------------------------------------------------------------------------
#include "OGS/V3.h"
#include "OGS/field.h"
#include "OGS/macros.h"
#include "OGS/projection.h"
#include "OGS/vtkFields.h"
#include "OGS/vtkOperations.h"

//----------------------------------------------------------------------------
void strsplit(const std::string &str, std::vector<std::string> &cont,
              char delim) {
  std::size_t current, previous = 0;
  current = str.find(delim);
  while (current != std::string::npos) {
    cont.push_back(str.substr(previous, current - previous));
    previous = current + 1;
    current = str.find(delim, previous);
  }
  cont.push_back(str.substr(previous, current - previous));
}

//----------------------------------------------------------------------------
OGSSelectPolygon::OGSSelectPolygon() {
  this->nProcs = 0;
  this->Invert = false;
  this->dfact = 1000.;
  this->projName = std::string("Mercator");

#ifdef PARAVIEW_USE_MPI
  this->Controller = NULL;
  this->SetController(vtkMultiProcessController::GetGlobalController());
#endif
}

//----------------------------------------------------------------------------
OGSSelectPolygon::~OGSSelectPolygon() {
  this->poly.clear();

#ifdef PARAVIEW_USE_MPI
  this->SetController(NULL);
#endif
}

//----------------------------------------------------------------------------
int OGSSelectPolygon::RequestData(vtkInformation *vtkNotUsed(request),
                                  vtkInformationVector **inputVector,
                                  vtkInformationVector *outputVector) {
  // Get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

// Stop all threads except from the master to execute
#ifdef PARAVIEW_USE_MPI
  if (this->procId > 0) return 1;
#endif

  // Get the input and output
  vtkDataSet *input =
      vtkDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(
      outInfo->Get(vtkDataObject::DATA_OBJECT()));

  this->UpdateProgress(0.0);

  // Obtain information on the projection (Metadata array)
  vtkStringArray *vtkmetadata = vtkStringArray::SafeDownCast(
      input->GetFieldData()->GetAbstractArray("Metadata"));
  this->dfact = (vtkmetadata != nullptr) ? std::stod(vtkmetadata->GetValue(2))
                                         : this->dfact;
  this->projName = (vtkmetadata != nullptr) ? vtkmetadata->GetValue(7)
                                            : std::string("Mercator");
  std::transform(this->projName.begin(), this->projName.end(),
                 this->projName.begin(), ::tolower);

  // Understand whether we are under cell or point data and compute the points
  // of the mesh
  int n_cell_vars = input->GetCellData()->GetNumberOfArrays();
  int n_point_vars = input->GetPointData()->GetNumberOfArrays();

  bool iscelld = (n_cell_vars > n_point_vars);

  // Obtain the cell centers or points
  OGS::V3::V3v xyz = (iscelld) ? getVTKCellCenters(input, this->dfact)
                               : getVTKCellPoints(input, this->dfact);

  // Generate a new field that will be used as cutting mask
  OGS::field::Field<FLDMASK> cutmask(xyz.len(), (FLDMASK)(1));

  // At this point "this->poly" should contain the polygon in degrees or be
  // empty
  OGS::Geom::Polygon<double> projected_poly;
  auto projected_poly_points = std::vector<OGS::Geom::Point<double>>();

  if (!this->poly.is_empty()) {
    // Loop the polygon and project it
    OGS::PROJ::Projection projection("degrees", this->projName);
    for (int ip = 0; ip < this->poly.get_npoints(); ++ip) {
      const auto p = this->poly.get_point(ip);
      projected_poly_points.push_back(projection.transform_point(p[0], p[1]));
    }
    projected_poly.set_points(std::move(projected_poly_points));

    // Compute the polygon bounding box for a faster performance
    this->poly.get_bbox();
    this->UpdateProgress(0.1);

// Loop and update cutting mask (Mesh loop, can be parallelized)
#pragma omp parallel shared(cutmask, xyz)
    {
      for (int ii = 0 + OMP_THREAD_NUM; ii < cutmask.get_n();
           ii += OMP_NUM_THREADS) {
        // Check if the point of the mesh is inside the polygon and set the
        // cutmask accordingly
        cutmask[ii][0] = this->poly.is_inside(OGS::Geom::Point<double>(
                             xyz[ii][0], xyz[ii][1], 0.))
                             ? 1
                             : 0;
      }
    }

    // Force ThresholdBetween to obtain values that are greater than 0
    if (this->Invert) {
      this->SetThresholdFunction(vtkThreshold::THRESHOLD_LOWER);
      this->SetLowerThreshold(0.5);
    } else {
      this->SetThresholdFunction(vtkThreshold::THRESHOLD_UPPER);
      this->SetUpperThreshold(0.5);
    }
  } else {
    this->SetThresholdFunction(vtkThreshold::THRESHOLD_UPPER);
    this->SetUpperThreshold(1e12);
    vtkWarningMacro("Problems parsing points! Check the input points.");
  }
  this->UpdateProgress(0.4);

  // Convert field to vtkArray and add it to input
  VTKMASK *vtkcutmask;
  vtkcutmask = VTK::createVTKfromField<VTKMASK, FLDMASK>("CutMask", cutmask);

  if (iscelld) {
    input->GetCellData()->AddArray(vtkcutmask);
    // Force to use the CutMask to produce the Threshold
    this->Superclass::SetInputArrayToProcess(
        0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "CutMask");
  } else {
    input->GetPointData()->AddArray(vtkcutmask);
    // Force to use the CutMask to produce the Threshold
    this->Superclass::SetInputArrayToProcess(
        0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "CutMask");
  }

  this->UpdateProgress(0.6);

  // Run the actual threshold filter
  this->Superclass::RequestData(nullptr, inputVector, outputVector);

  this->UpdateProgress(0.8);

  // Cleanup the output by deleting the CutMask and the basins mask
  if (iscelld)
    output->GetCellData()->RemoveArray("CutMask");
  else
    output->GetPointData()->RemoveArray("CutMask");

  vtkcutmask->Delete();

  // Return
  this->UpdateProgress(1.0);
  return 1;
}

//----------------------------------------------------------------------------
void OGSSelectPolygon::GetPolygon(const char *arg) {
  /*
          This function parses the text box and obtains the polygon to be
          used for cutting the mesh.
  */
  this->poly.clear();  // Make sure the polygon is empty before we start

  if (arg) {
    // Split the string by endline
    std::vector<std::string> aux;
    strsplit(std::string(arg), aux, '\n');

    // Loop the user inputed points
    bool error = false;
    std::vector<OGS::Geom::Point<double>> points;
    for (const std::string &str : aux) {
      // Split again by spaces
      std::vector<std::string> aux2;
      strsplit(str, aux2, ' ');
      // Check correctness
      if (str.empty() || aux2.size() != 2) {
        error = true;
        break;
      }
      // Convert to double
      double lon = std::stod(aux2[1]), lat = std::stod(aux2[0]);
      // Store point
      points.emplace_back(lon, lat, 0.);
    }
    // Only update polygon if the size is positive and there has been no error
    if (!points.empty() && !error) {
      points.push_back(
          points[0]);  // We need to close the polygon for a correct definition

      // Define the polygon
      this->poly.set_points(points);
    }
  }
  this->Superclass::Modified();
  this->Modified();
}