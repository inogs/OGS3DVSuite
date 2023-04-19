/*=========================================================================

  Program:   OGSSelectTools
  Module:    OGSSelectCoast.cxx

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "OGSSelectCoast.h"

#include <cstdint>

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
vtkCxxSetObjectMacro(OGSSelectCoast, Controller, vtkMultiProcessController);
#endif

vtkStandardNewMacro(OGSSelectCoast)
//----------------------------------------------------------------------------
#include "OGS/field.h"
#include "OGS/macros.h"
#include "OGS/vtkFields.h"

    //----------------------------------------------------------------------------
    void addCoasts(vtkDataArraySelection *CoastsDataArraySelection) {
  CoastsDataArraySelection->AddArray("Continental shelf");
  CoastsDataArraySelection->AddArray("Open Sea");
}

//----------------------------------------------------------------------------
OGSSelectCoast::OGSSelectCoast() {
  // Add the sub basins into the array
  this->CoastsDataArraySelection = vtkDataArraySelection::New();
  addCoasts(this->CoastsDataArraySelection);

  this->mask_field = static_cast<char const *>("coast mask");
  this->nProcs = 0;
  this->procId = 0;

#ifdef PARAVIEW_USE_MPI
  this->Controller = NULL;
  this->SetController(vtkMultiProcessController::GetGlobalController());
#endif
}

//----------------------------------------------------------------------------
OGSSelectCoast::~OGSSelectCoast() {
  this->CoastsDataArraySelection->Delete();
  this->mask_field = nullptr;

#ifdef PARAVIEW_USE_MPI
  this->SetController(NULL);
#endif
}

//----------------------------------------------------------------------------
void OGSSelectCoast::PrintSelf(ostream &os, vtkIndent indent) {
  this->Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------
int OGSSelectCoast::RequestData(vtkInformation *vtkNotUsed(request),
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

  // Recover basins mask as a field
  VTKMASK *vtkmask = nullptr;
  if (iscelld)
    vtkmask =
        VTKMASK::SafeDownCast(input->GetCellData()->GetArray(this->mask_field));
  else
    vtkmask = VTKMASK::SafeDownCast(
        input->GetPointData()->GetArray(this->mask_field));

  if (vtkmask == nullptr) {
    vtkErrorMacro("Cannot load mask field " << this->mask_field << "!");
    return 0;
  }

  // Recover basins mask as a field
  OGS::field::Field<FLDMASK> mask =
      VTK::createFieldfromVTK<VTKMASK, FLDMASK>(vtkmask);

  // Generate a new field (initialized at zero) that will be used as cutting
  // mask
  OGS::field::Field<FLDMASK> cutmask(mask.get_n(), 1);

  this->UpdateProgress(0.2);

// Loop and update cutting mask (Mesh loop, can be parallelized)
#pragma omp parallel shared(mask, cutmask)
  {
    for (int ii = 0 + OMP_THREAD_NUM; ii < mask.get_n();
         ii += OMP_NUM_THREADS) {
      cutmask[ii][0] = 0;
      // Loop on the basins array selection
      for (int bid = 0; bid < this->GetNumberOfCoastsArrays(); ++bid)
        // Test to zero is allowed when working with integers
        if (this->GetCoastsArrayStatus(this->GetCoastsArrayName(bid)) &&
            (mask[ii][0] - (bid + 1)) == 0) {
          cutmask[ii][0] = 1;
          continue;
        }
    }
  }

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

  this->UpdateProgress(0.4);

  // Force ThresholdBetween to obtain values that are greater than 0
  this->SetThresholdFunction(vtkThreshold::THRESHOLD_UPPER);
  this->SetUpperThreshold(0.5);

  this->UpdateProgress(0.6);

  // Run the actual threshold filter
  this->Superclass::RequestData(nullptr, inputVector, outputVector);

  this->UpdateProgress(0.8);

  // Cleanup the output by deleting the CutMask and the basins mask
  if (iscelld) {
    output->GetCellData()->RemoveArray("CutMask");
    output->GetCellData()->RemoveArray(this->mask_field);
  } else {
    output->GetPointData()->RemoveArray("CutMask");
    output->GetPointData()->RemoveArray(this->mask_field);
  }
  vtkcutmask->Delete();

  // Return
  this->UpdateProgress(1.0);
  return 1;
}


//----------------------------------------------------------------------------
void OGSSelectCoast::DisableAllCoastsArrays() {
  this->CoastsDataArraySelection->DisableAllArrays();
}

void OGSSelectCoast::EnableAllCoastsArrays() {
  this->CoastsDataArraySelection->EnableAllArrays();
}

int OGSSelectCoast::GetNumberOfCoastsArrays() {
  return this->CoastsDataArraySelection->GetNumberOfArrays();
}

const char *OGSSelectCoast::GetCoastsArrayName(int index) {
  if (index >= (int)this->GetNumberOfCoastsArrays() || index < 0)
    return nullptr;
  else
    return this->CoastsDataArraySelection->GetArrayName(index);
}

int OGSSelectCoast::GetCoastsArrayIndex(const char *name) {
  return this->CoastsDataArraySelection->GetArrayIndex(name);
}

int OGSSelectCoast::GetCoastsArrayStatus(const char *name) {
  return this->CoastsDataArraySelection->ArrayIsEnabled(name);
}

void OGSSelectCoast::SetCoastsArrayStatus(const char *name, int status) {
  if (status)
    this->CoastsDataArraySelection->EnableArray(name);
  else
    this->CoastsDataArraySelection->DisableArray(name);

  this->Modified();
}
