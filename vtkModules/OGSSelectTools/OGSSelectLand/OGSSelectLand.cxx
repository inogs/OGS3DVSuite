/*=========================================================================

  Program:   OGSSelectTools
  Module:    OGSSelectLand.cxx

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "OGSSelectLand.h"

#include "vtkCellData.h"
#include "vtkDataArraySelection.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkUnstructuredGrid.h"

#ifdef PARAVIEW_USE_MPI
#include "vtkMultiProcessController.h"
vtkCxxSetObjectMacro(OGSSelectLand, Controller, vtkMultiProcessController);
#endif

vtkStandardNewMacro(OGSSelectLand);

//----------------------------------------------------------------------------

#include "OGS/field.h"
#include "OGS/macros.h"
#include "OGS/vtkFields.h"

//----------------------------------------------------------------------------
OGSSelectLand::OGSSelectLand() {
  this->mask_field = static_cast<char const *>("land mask");
  this->nProcs = 0;
  this->procId = 0;

#ifdef PARAVIEW_USE_MPI
  this->Controller = NULL;
  this->SetController(vtkMultiProcessController::GetGlobalController());
#endif
}

//----------------------------------------------------------------------------
OGSSelectLand::~OGSSelectLand() {
  this->mask_field = nullptr;

#ifdef PARAVIEW_USE_MPI
  this->SetController(NULL);
#endif
}

//----------------------------------------------------------------------------
void OGSSelectLand::PrintSelf(ostream &os, vtkIndent indent) {
  this->Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------
int OGSSelectLand::RequestData(vtkInformation *vtkNotUsed(request),
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

  // Decide whether we have cell or point data
  int n_cell_vars = input->GetCellData()->GetNumberOfArrays();
  int n_point_vars = input->GetPointData()->GetNumberOfArrays();

  bool iscelld = (n_cell_vars > n_point_vars);

  if (iscelld) {
    // Force to use the CutMask to produce the Threshold
    this->Superclass::SetInputArrayToProcess(
        0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, this->mask_field);
  } else {
    // Force to use the CutMask to produce the Threshold
    this->Superclass::SetInputArrayToProcess(
        0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, this->mask_field);
  }

  this->UpdateProgress(0.4);

  this->SetThresholdFunction(vtkThreshold::THRESHOLD_UPPER);
  this->SetUpperThreshold(0.5);

  this->UpdateProgress(0.6);

  // Run the actual threshold filter
  this->Superclass::RequestData(nullptr, inputVector, outputVector);

  this->UpdateProgress(0.8);

  // Cleanup the output by deleting the CutMask and the basins mask
  if (iscelld) {
    output->GetCellData()->RemoveArray(this->mask_field);
  } else {
    output->GetPointData()->RemoveArray(this->mask_field);
  }

  // Return
  this->UpdateProgress(1.0);
  return 1;
}
