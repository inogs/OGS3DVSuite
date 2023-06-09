/*=========================================================================

  Program:   OGSDepthProfile
  Module:    OGSDepthProfile.cxx

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "OGSDepthProfile.h"

#include <string>

#include "vtkCell.h"
#include "vtkCellData.h"
#include "vtkCellLocator.h"
#include "vtkDataArray.h"
#include "vtkDataSet.h"
#include "vtkFloatArray.h"
#include "vtkGenericCell.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPointSet.h"
#include "vtkStaticCellLocator.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#define CELL_TOLERANCE_FACTOR_SQR 1e-6

vtkStandardNewMacro(OGSDepthProfile);
vtkCxxSetObjectMacro(OGSDepthProfile, CellLocatorPrototype,
                     vtkAbstractCellLocator);

#ifdef PARAVIEW_USE_MPI
#include "vtkMultiProcessController.h"
vtkCxxSetObjectMacro(OGSDepthProfile, Controller, vtkMultiProcessController);
#endif

//----------------------------------------------------------------------------
#include "OGS/V3.h"
#include "OGS/field.h"
#include "OGS/macros.h"
#include "OGS/vtkFields.h"

//----------------------------------------------------------------------------
OGSDepthProfile::OGSDepthProfile() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);

  this->PointList = nullptr;
  this->CellList = nullptr;
  this->nProcs = 0;
  this->procId = 0;

  this->CellLocatorPrototype = nullptr;

#ifdef PARAVIEW_USE_MPI
  this->Controller = NULL;
  this->SetController(vtkMultiProcessController::GetGlobalController());
#endif
}

//----------------------------------------------------------------------------
OGSDepthProfile::~OGSDepthProfile() {
  delete this->PointList;
  delete this->CellList;

  this->OGSDepthProfile::SetCellLocatorPrototype(nullptr);

#ifdef PARAVIEW_USE_MPI
  this->SetController(NULL);
#endif
}

//----------------------------------------------------------------------------
int OGSDepthProfile::RequestData(vtkInformation *vtkNotUsed(request),
                                 vtkInformationVector **inputVector,
                                 vtkInformationVector *outputVector) {
// Stop all threads except from the master to execute
#ifdef PARAVIEW_USE_MPI
  if (this->procId > 0) return 1;
#endif

  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *sourceInfo = inputVector[1]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input and output
  // input contains the interpolating line information (number of points, etc)
  vtkDataSet *input =
      vtkDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  // Source contains the data where to interpolate from
  vtkDataSet *source =
      vtkDataSet::SafeDownCast(sourceInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkDataSet *output =
      vtkDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  // First, copy the input to the output as a starting point
  output->CopyStructure(input);
  output->GetPointData()->SetCopyAttribute(vtkDataSetAttributes::SCALARS, 2,
                                           vtkDataSetAttributes::INTERPOLATE);
  output->GetFieldData()->PassData(source->GetFieldData());

  // If there is data to interpolate, begin the interpolation
  if (source) {
    this->Initialize(input, source, output);
    this->Interpolate(input, source, output);
  }

  return 1;
}

//----------------------------------------------------------------------------
int OGSDepthProfile::RequestInformation(vtkInformation *vtkNotUsed(request),
                                        vtkInformationVector **inputVector,
                                        vtkInformationVector *outputVector) {
/* SET UP THE PARALLEL CONTROLLER

        The MPI threads come initialized by the ParaView server. Here
        we set up the environment for this filter.

*/
#ifdef PARAVIEW_USE_MPI
  if (this->Controller->GetNumberOfProcesses() > 1) {
    this->nProcs = this->Controller->GetNumberOfProcesses();
    this->procId = this->Controller->GetLocalProcessId();
  }

  // Stop all threads except from the master to execute
  if (this->procId > 0) return 1;
#endif

  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *sourceInfo = inputVector[1]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  outInfo->CopyEntry(sourceInfo,
                     vtkStreamingDemandDrivenPipeline::TIME_STEPS());
  outInfo->CopyEntry(sourceInfo,
                     vtkStreamingDemandDrivenPipeline::TIME_RANGE());

  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),
               inInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT()),
               6);

  return 1;
}

//----------------------------------------------------------------------------
int OGSDepthProfile::RequestUpdateExtent(vtkInformation *vtkNotUsed(request),
                                         vtkInformationVector **inputVector,
                                         vtkInformationVector *outputVector) {
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *sourceInfo = inputVector[1]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  int usePiece = 0;

  // What ever happened to CopyUpdateExtent in vtkDataObject?
  // Copying both piece and extent could be bad.  Setting the piece
  // of a structured data set will affect the extent.
  vtkDataObject *output = outInfo->Get(vtkDataObject::DATA_OBJECT());
  if (output && (!strcmp(output->GetClassName(), "vtkUnstructuredGrid") ||
                 !strcmp(output->GetClassName(), "vtkPolyData")))
    usePiece = 1;

  inInfo->Set(vtkStreamingDemandDrivenPipeline::EXACT_EXTENT(), 1);

  sourceInfo->Remove(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT());
  if (sourceInfo->Has(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT()))
    sourceInfo->Set(
        vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(),
        sourceInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT()), 6);

  // SpatialMatch does not exist in this implementation
  sourceInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(), 0);
  sourceInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(),
                  1);
  sourceInfo->Set(
      vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(), 0);

  if (usePiece) {
    inInfo->Set(
        vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(),
        outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER()));
    inInfo->Set(
        vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(),
        outInfo->Get(
            vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES()));
    inInfo->Set(
        vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(),
        outInfo->Get(
            vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS()));
  } else {
    inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(),
                outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT()),
                6);
  }

  return 1;
}

//----------------------------------------------------------------------------
void OGSDepthProfile::Initialize(vtkDataSet *input, vtkDataSet *source,
                                 vtkDataSet *output) {
  // In this function, we will allocate space in output for the arrays to be
  // interpolated We will recover them and create VTKARRAYS

  int npoints = input->GetNumberOfPoints();
  vtkDataArray *vtkDArray;

  // Copy Cell Arrays
  for (int varId = 0; varId < source->GetCellData()->GetNumberOfArrays();
       ++varId) {
    // Recover the array and the array name
    vtkDArray = source->GetCellData()->GetArray(varId);
    std::string arrName = vtkDArray->GetName();

    // Do not work with the basins, coasts mask, e1, e2 or e3
    if (std::string("basins mask") == arrName) continue;
    if (std::string("coast mask") == arrName) continue;
    if (std::string("land mask") == arrName) continue;
    if (std::string("Okubo-Weiss mask") == arrName) continue;
    if (std::string("Q-criterion mask") == arrName) continue;
    if (std::string("e1") == arrName) continue;
    if (std::string("e2") == arrName) continue;
    if (std::string("e3") == arrName) continue;

    OGS::field::Field<FLDARRAY> array(npoints,
                                      vtkDArray->GetNumberOfComponents(), 0.);
    auto *vtkArray =
        VTK::createVTKfromField<VTKARRAY, FLDARRAY>(arrName.c_str(), array);
    output->GetPointData()->AddArray(vtkArray);
  }

  // Copy Point Arrays
  for (int varId = 0; varId < source->GetPointData()->GetNumberOfArrays();
       ++varId) {
    // Recover the array and the array name
    vtkDArray = source->GetPointData()->GetArray(varId);
    std::string arrName = vtkDArray->GetName();

    // Do not work with the basins, coasts mask, e1, e2 or e3
    if (std::string("basins mask") == arrName) continue;
    if (std::string("coast mask") == arrName) continue;
    if (std::string("land mask") == arrName) continue;
    if (std::string("Okubo-Weiss mask") == arrName) continue;
    if (std::string("Q-criterion mask") == arrName) continue;
    if (std::string("e1") == arrName) continue;
    if (std::string("e2") == arrName) continue;
    if (std::string("e3") == arrName) continue;

    OGS::field::Field<FLDARRAY> array(npoints,
                                      vtkDArray->GetNumberOfComponents(), 0.);
    auto *vtkArray =
        VTK::createVTKfromField<VTKARRAY, FLDARRAY>(arrName.c_str(), array);
    output->GetPointData()->AddArray(vtkArray);
  }
}

void OGSDepthProfile::Interpolate(vtkDataSet *input, vtkDataSet *source,
                                  vtkDataSet *output) {
  // Update progress
  this->UpdateProgress(0.);

  // #pragma omp parallel shared(input,source,output)
  //{
  //  Preallocate weights
  double *weights;
  weights = new double[source->GetMaxCellSize()];

  // Create the cell locator object
  vtkCellLocator *cellLocator = vtkCellLocator::New();
  cellLocator->SetDataSet(source);
  cellLocator->BuildLocator();

  // Recover the number of points in the input
  int npoints = input->GetNumberOfPoints();

  // Loop the number of points
  vtkNew<vtkGenericCell> gcell;
  VTKARRAY *outArray, *srcArray;

  // for (int pId = omp_get_thread_num(); pId < npoints;
  // pId+=omp_get_num_threads()) {
  for (int pId = 0; pId < npoints; ++pId) {
    // Get the xyz coordinate of the point in the input dataset
    // then, find the cell id that contains xyz
    vtkIdType cellId = 0;
    int subId = 0;
    double dist2;
    OGS::V3::V3 xyz, pcoords, closestPoint;

    input->GetPoint(pId, &xyz[0]);
    cellLocator->FindClosestPoint(&xyz[0], &closestPoint[0], gcell.GetPointer(),
                                  cellId, subId, dist2);

    // Evaluate interpolation weights
    if (cellId >= 0) {
      // Compute a tolerance proportional to the cell length.
      gcell->EvaluatePosition(&xyz[0], &closestPoint[0], subId, &pcoords[0],
                              dist2, weights);
      // Abort if the distance is too big
      if (dist2 > (gcell->GetLength2() * CELL_TOLERANCE_FACTOR_SQR)) continue;

      // Interpolate the point
      for (int varId = 0; varId < output->GetPointData()->GetNumberOfArrays();
           ++varId) {
        // Recover data arrays
        outArray =
            VTKARRAY::SafeDownCast(output->GetPointData()->GetArray(varId));
        srcArray = VTKARRAY::SafeDownCast(
            source->GetCellData()->GetArray(outArray->GetName()));

        if (srcArray) {
          // We have a CellData array
          output->GetPointData()->CopyTuple(srcArray, outArray, cellId, pId);
        } else {
          // We have a PointData array
          srcArray = VTKARRAY::SafeDownCast(
              source->GetPointData()->GetArray(outArray->GetName()));
          outArray->InterpolateTuple(pId, gcell->PointIds, srcArray, weights);
        }
      }
    }
    this->UpdateProgress(0. + 1. / (npoints - 0.) * pId);
  }
  delete[] weights;
  cellLocator->Delete();
  //}

  // Update progress
  this->UpdateProgress(1.);
}

//----------------------------------------------------------------------------
void OGSDepthProfile::SetSourceConnection(vtkAlgorithmOutput *algOutput) {
  this->SetInputConnection(1, algOutput);
}

void OGSDepthProfile::SetSourceData(vtkDataObject *input) {
  this->SetInputData(1, input);
}

vtkDataObject *OGSDepthProfile::GetSource() {
  if (this->GetNumberOfInputConnections(1) < 1) return nullptr;

  return this->GetExecutive()->GetInputData(1, 0);
}

//----------------------------------------------------------------------------
