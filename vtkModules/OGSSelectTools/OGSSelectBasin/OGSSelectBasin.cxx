/*=========================================================================

  Program:   OGSSelectTools
  Module:    OGSSelectBasin.cxx

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "OGSSelectBasin.h"

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
vtkCxxSetObjectMacro(OGSSelectBasin, Controller, vtkMultiProcessController);
#endif

vtkStandardNewMacro(OGSSelectBasin);

//----------------------------------------------------------------------------
#include "OGS/field.h"
#include "OGS/macros.h"
#include "OGS/vtkFields.h"

//----------------------------------------------------------------------------
void addSubBasins(vtkDataArraySelection *BasinsDataArraySelection) {
  BasinsDataArraySelection->AddArray("Alboran Sea");
  BasinsDataArraySelection->AddArray("South Western Mediterranean (west)");
  BasinsDataArraySelection->AddArray("South Western Mediterranean (east)");
  BasinsDataArraySelection->AddArray("North Western Mediterranean");
  BasinsDataArraySelection->AddArray("Northern Tyrrhenian");
  BasinsDataArraySelection->AddArray("Southern Tyrrhenian");
  BasinsDataArraySelection->AddArray("Northern Adriatic");
  BasinsDataArraySelection->AddArray("Southern Adriatic");
  BasinsDataArraySelection->AddArray("Aegean Sea");
  BasinsDataArraySelection->AddArray("Western Ionian");
  BasinsDataArraySelection->AddArray("Eastern Ionian");
  BasinsDataArraySelection->AddArray("Northern Ionian");
  BasinsDataArraySelection->AddArray("Western Levantine");
  BasinsDataArraySelection->AddArray("Northern Levantine");
  BasinsDataArraySelection->AddArray("Southern Levantine");
  BasinsDataArraySelection->AddArray("Eastern Levantine");
}

//----------------------------------------------------------------------------
OGSSelectBasin::OGSSelectBasin() {
  // Add the sub basins into the array
  this->BasinsDataArraySelection = vtkDataArraySelection::New();
  addSubBasins(this->BasinsDataArraySelection);

  this->mask_field = "basins mask";
  this->nProcs = 0;
  this->procId = 0;
  this->basins_field = 0;

#ifdef PARAVIEW_USE_MPI
  this->Controller = NULL;
  this->SetController(vtkMultiProcessController::GetGlobalController());
#endif
}

//----------------------------------------------------------------------------
OGSSelectBasin::~OGSSelectBasin() {
  this->BasinsDataArraySelection->Delete();
  this->mask_field = nullptr;

#ifdef PARAVIEW_USE_MPI
  this->SetController(NULL);
#endif
}

//----------------------------------------------------------------------------
void OGSSelectBasin::PrintSelf(ostream &os, vtkIndent indent) {
  this->Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------
int OGSSelectBasin::RequestData(vtkInformation *vtkNotUsed(request),
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

  OGS::field::Field<FLDMASK> mask =
      VTK::createFieldfromVTK<VTKMASK, FLDMASK>(vtkmask);

  // Generate a new field (initialized at zero) that will be used as cutting
  // mask
  OGS::field::Field<FLDMASK> cutmask(mask.get_n(), 1), bfield(mask.get_n(), 1);

  this->UpdateProgress(0.2);

// Loop and update cutting mask (Mesh loop, can be parallelized)
#pragma omp parallel shared(mask, cutmask, bfield)
  {
    for (int ii = 0 + OMP_THREAD_NUM; ii < mask.get_n();
         ii += OMP_NUM_THREADS) {
      cutmask[ii][0] = 0;
      // Loop on the basins array selection
      for (int bid = 0; bid < this->GetNumberOfBasinsArrays(); ++bid)
        if (this->GetBasinsArrayStatus(this->GetBasinsArrayName(bid)) &&
            mask[ii][bid]) {
          cutmask[ii][0] = 1;
          bfield[ii][0] = (FLDMASK)(bid);
          continue;
        }
    }
  }

  // Convert field to vtkArray and add it to input
  VTKMASK *vtkcutmask, *vtkbfield;
  vtkcutmask = VTK::createVTKfromField<VTKMASK, FLDMASK>("CutMask", cutmask);
  if (this->basins_field)
    vtkbfield =
        VTK::createVTKfromField<VTKMASK, FLDMASK>("basins field", bfield);

  if (iscelld) {
    input->GetCellData()->AddArray(vtkcutmask);
    if (this->basins_field) input->GetCellData()->AddArray(vtkbfield);
    // Force to use the CutMask to produce the Threshold
    this->Superclass::SetInputArrayToProcess(
        0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "CutMask");
  } else {
    input->GetPointData()->AddArray(vtkcutmask);
    if (this->basins_field) input->GetPointData()->AddArray(vtkbfield);
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
void OGSSelectBasin::DisableAllBasinsArrays() {
  this->BasinsDataArraySelection->DisableAllArrays();
}

void OGSSelectBasin::EnableAllBasinsArrays() {
  this->BasinsDataArraySelection->EnableAllArrays();
}

int OGSSelectBasin::GetNumberOfBasinsArrays() {
  return this->BasinsDataArraySelection->GetNumberOfArrays();
}

const char *OGSSelectBasin::GetBasinsArrayName(int index) {
  if (index >= (int)this->GetNumberOfBasinsArrays() || index < 0)
    return nullptr;
  else
    return this->BasinsDataArraySelection->GetArrayName(index);
}

int OGSSelectBasin::GetBasinsArrayIndex(const char *name) {
  return this->BasinsDataArraySelection->GetArrayIndex(name);
}

int OGSSelectBasin::GetBasinsArrayStatus(const char *name) {
  return this->BasinsDataArraySelection->ArrayIsEnabled(name);
}

void OGSSelectBasin::SetBasinsArrayStatus(const char *name, int status) {
  if (status)
    this->BasinsDataArraySelection->EnableArray(name);
  else
    this->BasinsDataArraySelection->DisableArray(name);

  this->Modified();
}
