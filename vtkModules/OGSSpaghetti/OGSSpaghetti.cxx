/*=========================================================================

  Program:   OGSSpaghetti
  Module:    vtkOGSSpaghetti.cxx

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "OGSSpaghetti.h"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <ctime>

#include "vtkCell.h"
#include "vtkCellData.h"
#include "vtkCellLocator.h"
#include "vtkDataArray.h"
#include "vtkDataSet.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkGenericCell.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPointSet.h"
#include "vtkStaticCellLocator.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStringArray.h"
#include "vtkTable.h"
#include "vtkTypeUInt8Array.h"

vtkStandardNewMacro(OGSSpaghetti);
vtkCxxSetObjectMacro(OGSSpaghetti, CellLocatorPrototype,
                     vtkAbstractCellLocator);

#ifdef PARAVIEW_USE_MPI
#include "vtkMultiProcessController.h"
vtkCxxSetObjectMacro(vtkOGSSpaghetti, Controller, vtkMultiProcessController);
#endif

#define CELL_TOLERANCE_FACTOR_SQR 1e-6

//----------------------------------------------------------------------------
#include "OGS/OGS3DVSuite.h"
#include "OGS/fieldOperations.h"
#include "OGS/macros.h"
#include "OGS/netcdfio.h"
#include "OGS/utilities.h"
#include "OGS/vtkFields.h"
#include "OGS/vtkOperations.h"

//----------------------------------------------------------------------------
void RecoverMasterFileName(std::string &fname, vtkDataSet *input) {
  // Recover the master file name from the metadata array
  // Return whether we need to stop executing or not

  vtkStringArray *vtkmetadata = vtkStringArray::SafeDownCast(
      input->GetFieldData()->GetAbstractArray("Metadata"));

  // If successful, recover the file name
  if (vtkmetadata) fname = vtkmetadata->GetValue(4);
}

//----------------------------------------------------------------------------
int ComputeCellId(vtkDataSet *input, vtkDataSet *source) {
  int cId = -1;

  // We need the weights as they are part of a function
  double *weights;
  weights = new double[source->GetMaxCellSize()];

  // Create the cell locator object
  vtkCellLocator *cellLocator = vtkCellLocator::New();
  cellLocator->SetDataSet(source);
  cellLocator->BuildLocator();

  // Get the xyz coordinate of the point in the input dataset
  // then, find the cell id that contains xyz
  vtkIdType cellId = 0;
  int subId = 0;
  double dist2;
  OGS::V3::V3 xyz, pcoords, closestPoint;
  vtkNew<vtkGenericCell> gcell;

  input->GetPoint(0, &xyz[0]);
  cellLocator->FindClosestPoint(&xyz[0], &closestPoint[0], gcell.GetPointer(),
                                cellId, subId, dist2);

  // Evaluate interpolation weights
  if (cellId >= 0) {
    // Compute a tolerance proportional to the cell length.
    gcell->EvaluatePosition(&xyz[0], &closestPoint[0], subId, &pcoords[0],
                            dist2, weights);
    // Abort if the distance is too big
    if (dist2 > (gcell->GetLength2() * CELL_TOLERANCE_FACTOR_SQR)) return cId;
    // Store the cell id
    cId = (int)(cellId);
  }

  return cId;
}

//----------------------------------------------------------------------------
void ComputeCellIds(std::vector<int> &cellIds, vtkDataSet *input,
                    vtkDataSet *source) {
  // We need the weights as they are part of a function
  double *weights;
  weights = new double[source->GetMaxCellSize()];

  // Create the cell locator object
  vtkCellLocator *cellLocator = vtkCellLocator::New();
  cellLocator->SetDataSet(source);
  cellLocator->BuildLocator();

  // Recover the number of points in the input
  int npoints = input->GetNumberOfPoints();
  cellIds.resize(npoints, -1);

  // Loop the number of points
  vtkNew<vtkGenericCell> gcell;
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
      // Store the cell id
      cellIds[pId] = (int)(cellId);
    }
  }
  delete[] weights;
  cellLocator->Delete();
}

//----------------------------------------------------------------------------
void BuildTimeList(OGS::Time::TimeList &TL, vtkInformation *Info) {
  // Build a TimeList object using the pipeline temporal data.
  // This TimeList will be later used for computing the averages
  // given a TimeRequestor. The metadata array might not be available
  // from the beginning.

  int ntsteps = Info->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
  double *tsteps;
  tsteps = Info->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());

  // Create a TimeObjectList where to store all the instants
  OGS::Time::TimeObjectList TOL(ntsteps);

  // Iterate the number of steps and set the values of the list
  for (int ii = 0; ii < ntsteps; ++ii) {
    // Convert to struct tm
    time_t time = std::chrono::system_clock::to_time_t(
        std::chrono::system_clock::time_point(
            std::chrono::duration_cast<std::chrono::seconds>(
                std::chrono::duration<double>(tsteps[ii]))));
    struct tm tm = *localtime(&time);

    // Set up the TimeObjectList
    TOL[ii] = OGS::Time::TimeObject(tm);
  }

  // Sort and print the list (for debugging purposes)
  TOL.sort();
  //	printf("TOL: %s\n",TOL.as_string("%Y%m%d").c_str());

  // Now create the TimeList from the TimeObjectList
  TL = OGS::Time::TimeList(TOL);
  //	printf("Defined list of %d elements: %s ...
  //%s\n",TL.len(),TL[0].as_string("%Y-%m-%d
  //%H:%M:%S").c_str(),TL[-1].as_string("%Y-%m-%d %H:%M:%S").c_str());
}

//----------------------------------------------------------------------------
OGSSpaghetti::OGSSpaghetti() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);

  this->field = nullptr;
  this->FolderName = NULL;
  this->bmask_field = NULL;
  this->cmask_field = NULL;
  this->TL_computed = false;
  this->use_average = false;
  this->use_files = false;
  this->CurrentTimeIndex = 0;
  this->sId = 0;
  this->per_coast = 0;
  this->procId = 0;
  this->nProcs = 0;
  this->isReqInfo = false;
  this->epsi = 1.e-3;
  this->dfact = 1000.;

  this->CellLocatorPrototype = nullptr;

#ifdef PARAVIEW_USE_MPI
  this->Controller = nullptr;
  this->SetController(vtkMultiProcessController::GetGlobalController());
#endif
}

//----------------------------------------------------------------------------
OGSSpaghetti::~OGSSpaghetti() {
  this->Setfield(NULL);
  this->SetFolderName(NULL);
  this->Setbmask_field(NULL);
  this->Setcmask_field(NULL);

  this->OGSSpaghetti::SetCellLocatorPrototype(nullptr);

#ifdef PARAVIEW_USE_MPI
  this->SetController(nullptr);
#endif
}

//----------------------------------------------------------------------------
int OGSSpaghetti::FillInputPortInformation(int vtkNotUsed(port),
                                           vtkInformation *info) {
  info->Remove(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE());
  info->Append(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
  return 1;
}

//----------------------------------------------------------------------------
int OGSSpaghetti::RequestInformation(vtkInformation *vtkNotUsed(request),
                                     vtkInformationVector **inputVector,
                                     vtkInformationVector *outputVector) {
  // Grab the input
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkInformation *srcInfo = inputVector[1]->GetInformationObject(0);

// Handle the parallel controller
#ifdef PARAVIEW_USE_MPI
  if (this->Controller->GetNumberOfProcesses() > 1) {
    this->procId = this->Controller->GetLocalProcessId();
    this->nProcs = this->Controller->GetNumberOfProcesses();
  }
#endif

  // Compute the TimeList if it hasn't been computed
  if (!this->TL_computed) {
    BuildTimeList(this->TL, srcInfo);
    this->TL_computed = true;
  }

  // The output data of this filter has no time associated with it. It is the
  // result of computations that happen over all time.
  outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
  outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_RANGE());

  this->isReqInfo = true;
  return 1;
}

//------------------------------------------------------------------------------
int OGSSpaghetti::RequestUpdateExtent(
    vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector,
    vtkInformationVector *vtkNotUsed(outputVector)) {
  vtkInformation *inInfo = inputVector[1]->GetInformationObject(0);

  // The RequestData method will tell the pipeline executive to iterate the
  // upstream pipeline to get each time step in order.  The executive in turn
  // will call this method to get the extent request for each iteration (in this
  // case the time step).
  double *inTimes = inInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());

  if (inTimes) {
    int inst = (this->instants.size() == 0)
                   ? this->CurrentTimeIndex
                   : this->instants[this->CurrentTimeIndex];
    //		printf("Time: <%d> %s\n",inst,this->TL[inst].as_string("%Y-%m-%d
    //%H:%M:%S").c_str());
    inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(),
                inTimes[inst]);
  }

  return 1;
}

//----------------------------------------------------------------------------
int OGSSpaghetti::RequestData(vtkInformation *request,
                              vtkInformationVector **inputVector,
                              vtkInformationVector *outputVector) {
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *srcInfo = inputVector[1]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // input contains the interpolating line information (number of points, etc)
  vtkDataSet *input =
      vtkDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  // Source contains the mesh where to interpolate from
  vtkDataSet *source =
      vtkDataSet::SafeDownCast(srcInfo->Get(vtkDataObject::DATA_OBJECT()));
  // Output is a vtkTable with the interpolated data per each timestep
  vtkTable *output =
      vtkTable::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  if (this->CurrentTimeIndex == 0) {
    // Check input field
    if (std::string(this->field) == std::string("coast mask")) {
      vtkErrorMacro("Wrong input field! Aborting...");
      return 0;
    }
    if (std::string(this->field) == std::string("basins mask")) {
      vtkErrorMacro("Wrong input field! Aborting...");
      return 0;
    }
    if (std::string(this->field) == std::string("land mask")) {
      vtkErrorMacro("Wrong input field! Aborting...");
      return 0;
    }
    if (std::string(this->field) == std::string("e1")) {
      vtkErrorMacro("Wrong input field! Aborting...");
      return 0;
    }
    if (std::string(this->field) == std::string("e2")) {
      vtkErrorMacro("Wrong input field! Aborting...");
      return 0;
    }
    if (std::string(this->field) == std::string("e3")) {
      vtkErrorMacro("Wrong input field! Aborting...");
      return 0;
    }

    // Master computes the instants to loop
    if (this->procId == 0) {
      // Compute the TimeList if it hasn't been computed
      if (!this->TL_computed) {
        BuildTimeList(this->TL, srcInfo);
        this->TL_computed = true;
      }
      //			printf("Defined list of %d elements: %s ...
      //%s\n",this->TL.len(),
      // this->TL[0].as_string("%Y-%m-%d %H:%M:%S").c_str(),
      // this->TL[-1].as_string("%Y-%m-%d %H:%M:%S").c_str());

      // At this point check that the TimeInterval length is not zero
      // if so, crash and stop, there is no sense to compute an average.
      //			printf("Time interval is %s of length %ld\n",
      //				this->TI.as_string("%Y%m%d-%H:%M:%S").c_str(),this->TI.length());
      if (this->TI.length() == 0) {
        vtkErrorMacro("TimeInterval is <"
                      << this->TI.as_string("%Y%m%d-%H:%M:%S")
                      << ">! Cannot continue...");
        return 0;
      }

      // We can now create a GenericRequestor with this TimeInterval.
      // We will use it to obtain all the instants and average weights.
      OGS::Time::Requestor req(this->TI, "%Y%m%d-%H:%M:%S");
      this->instants.clear();
      this->weights.clear();  // make sure vectors are empty

      if (this->TL.select(&req, this->instants, this->weights) == TIME_ERR) {
        vtkErrorMacro(
            "Error selecting instants with Requestor! Cannot continue...");
        return 0;
      }
      //			for (int ii=0; ii<instants.size();++ii)
      //				printf("ind: %d, weight: %f,
      //%s\n",this->instants[ii],this->weights[ii],
      //					this->TL[instants[ii]].as_string("%Y-%m-%d
      //%H:%M:%S").c_str());
    }

// Broadcast instants and weights to all the processors if needed
#ifdef PARAVIEW_USE_MPI
    if (this->nProcs > 1) {
      // Array length
      int arr_len = (int)(this->instants.size());
      this->Controller->Broadcast(&arr_len, 1, 0);
      // Allocate
      if (this->procId > 0) {
        this->instants.resize(arr_len, 0);
        this->weights.resize(arr_len, 0.);
      }
      // Broadcast array
      this->Controller->Broadcast(this->instants.data(), arr_len, 0);
    }
#endif
  }

  // Decide which algorithm to use
  int retval = 0;
  if (this->use_average) {
    int ntsteps =
        srcInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
    retval = this->AveragesIterationAlgorithm(ntsteps, request, input, source,
                                              output);
  } else {
    if (this->use_files)
      retval = this->FileIterationAlgorithm(request, input, source, output);
    else
      retval = this->PipelineIterationAlgorithm(request, input, source, output);
  }

  this->UpdateProgress(1.);
  return retval;
}

//----------------------------------------------------------------------------
int OGSSpaghetti::PipelineIterationAlgorithm(vtkInformation *request,
                                             vtkDataSet *input,
                                             vtkDataSet *source,
                                             vtkTable *output) {
// Stop all threads except from the master to execute
#ifdef PARAVIEW_USE_MPI
  if (this->procId > 0) return 1;
#endif

  // Compute the cell list corresponding to the
  // input interpolating point or line.
  int cList = ComputeCellId(input, source);

  this->UpdateProgress(.25);

  // Current time instant
  int itime = this->instants[this->CurrentTimeIndex];
  //	double weight = this->weights[this->CurrentTimeIndex];
  //	printf("Time: <%d,%d>
  //%s\n",this->CurrentTimeIndex,itime,this->TL[itime].as_string("%Y-%m-%d
  //%H:%M:%S").c_str());

  // Recover input array values
  VTKARRAY *vtkInArray;
  vtkInArray =
      VTKARRAY::SafeDownCast(source->GetCellData()->GetArray(this->field));
  if (vtkInArray == nullptr)
    vtkInArray =
        VTKARRAY::SafeDownCast(source->GetPointData()->GetArray(this->field));

  if (vtkInArray == nullptr) {
    vtkErrorMacro("Cannot load variable <" << this->field << ">!");
    return 0;
  }

  OGS::field::Field<FLDARRAY> inArray;
  inArray = VTK::createFieldfromVTK<VTKARRAY, FLDARRAY>(vtkInArray);

  // Create output array
  // We store each instant as a column on the output
  OGS::field::Field<FLDARRAY> outArray(1, 1, 0.);

  // Set the output value
  if (inArray.get_m() > 1) {
    double mod = 0.f;
    for (int jj = 0; jj < inArray.get_m(); ++jj)
      mod += (cList > 0) ? inArray[cList][0] * inArray[cList][0] : 0.f;
    outArray[0][0] = std::sqrt(mod);
  } else {
    outArray[0][0] = (cList > 0) ? inArray[cList][0] : 0.;
  }
  this->UpdateProgress(.5);

  // Write column to the output table
  VTKARRAY *vtkColumn;
  vtkColumn = VTK::createVTKfromField<VTKARRAY, FLDARRAY>(
      this->TL[itime].as_string("%Y%m%d-%H:%M:%S").c_str(), outArray);
  output->AddColumn(vtkColumn);
  vtkColumn->Delete();

  this->UpdateProgress(.75);

  // Increment the time-step
  this->CurrentTimeIndex++;

  // Continue executing or finish
  if (this->CurrentTimeIndex < this->instants.size()) {
    // There is still more to do.
    request->Set(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 1);
  } else {
    // We are done. Finish up.
    request->Remove(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING());
    this->CurrentTimeIndex = 0;
  }

  this->UpdateProgress(1.);
  return 1;
}

//----------------------------------------------------------------------------
int OGSSpaghetti::FileIterationAlgorithm(vtkInformation *request,
                                         vtkDataSet *input, vtkDataSet *source,
                                         vtkTable *output) {
  // Check inputs
  if (!source->IsA("vtkRectilinearGrid")) {
    vtkErrorMacro(
        "A rectilinear grid is needed for this method! Cannot continue...");
    return 0;
  }

  if (source->GetPointData()->GetNumberOfArrays() > 0) {
    vtkErrorMacro("Data contains point arrays! Cannot continue...");
    return 0;
  }

  /* INITIALIZATION PHASE

          Master computes the cell list corresponding to the spaghetti point.
          Recover (and broadcast) the file name of the OGS master file.

          Initialize arrays.

  */

  std::string FileName, projName;
  int cList, npoints;

#ifdef PARAVIEW_USE_MPI

  // Thread 0 has all the information in input and output, therefore is the
  // only computing the array which will later be broadcasted to all ranks
  if (this->procId == 0) {
    npoints = source->GetNumberOfCells();
    // Compute the cell list corresponding to the Hovmoeller interpolating line.
    cList = ComputeCellId(input, source);
    // Recover the master file name from source
    RecoverMasterFileName(FileName, source);
  }

  // Broadcast the information to all other threads if applicable
  if (this->nProcs > 1) {
    // Broadcast master file name
    char buff[128] = "";
    sprintf(buff, "%s", FileName.c_str());
    this->Controller->Broadcast(buff, FileName.length(), 0);
    FileName = std::string(buff);
    // Broadcast cList
    this->Controller->Broadcast(&cList, 1, 0);
    // Broadcast number of points
    this->Controller->Broadcast(&npoints, 1, 0);
  }

#else

  // This is the normal non-parallel algorithm
  npoints = static_cast<int>(source->GetNumberOfCells());
  cList = ComputeCellId(input, source);
  RecoverMasterFileName(FileName, source);

#endif

  // Read the OGS file to be able to generate the paths to the variable files
  OGS::Simulation *current_simulation = nullptr;
  try {
    current_simulation = new OGS::Simulation(FileName);
  } catch (...) {
    vtkErrorMacro("Unable to read file " + FileName);
    return 0;
  }

  // Read mesh. ProjID is not important since only the dimensions
  // and the e1, e2 and e3 arrays are necessary.
  try {
    current_simulation->load_mesh(0);
  } catch (...) {
    vtkErrorMacro("Problems reading the mesh!\nAborting.");
    delete current_simulation;
    return 0;
  }

  this->UpdateProgress(0.1);

  /* TRAVERSE PHASE

          Open and read the NetCDF files containing the variables that we
     previously defined. Store the data in a master field for posterior
     reduction.

  */
  // Parallel partition, define starting and ending ranges for each worker
  int ii_start = 0, ii_end = this->instants.size();

#ifdef PARAVIEW_USE_MPI
  if (this->nProcs > 1) {
    // Everyone should have the range of start and end times
    // We must split equally among the threads. We must also control
    // that the number of threads is less or equal than the intervals
    // requested.
    int range = ii_end - ii_start;
    if (this->nProcs < range) {
      // We split normally among processes assuming no remainder
      int rangePerProcess = std::floor(range / this->nProcs);
      ii_start = ii_start + this->procId * rangePerProcess;
      ii_end = ii_start + rangePerProcess;
      // Handle the remainder
      int remainder = range - rangePerProcess * this->nProcs;
      if (remainder > this->procId) {
        ii_start += this->procId;
        ii_end += this->procId + 1;
      } else {
        ii_start += remainder;
        ii_end += remainder;
      }
    } else {
      // Each process will forcefully conduct one instant.
      ii_start = (this->procId < ii_end) ? this->procId : ii_end;
      ii_end = (this->procId < ii_end) ? this->procId + 1 : ii_end;
    }
  }
#endif

  // Loop the instants, prepare arrays and variables for the computation of the
  // mean
  OGS::field::Field<FLDARRAY> arrayTemp,
      spag_array(this->instants.size(), 1, 0.);

  std::vector<std::string> vel_vars;
  const auto var_indices = current_simulation->get_var_indices(this->field);

  if (!var_indices) {
    vtkErrorMacro("Can not find variable named " + std::string(this->field));
    delete current_simulation;
    return 0;
  }

  for (int inst = ii_start; inst < ii_end; ++inst) {
    int itime = this->instants[inst];

    // Deal with the velocity
    if (std::string(this->field) == std::string("Velocity")) {
      OGS::field::Field<FLDARRAY> tmp(npoints, 3, 0.);

      OGS::utilities::strsplit(current_simulation->var_vname(
                                   var_indices->first, var_indices->second),
                               ",", vel_vars);

      const auto var_path = current_simulation->var_path(
          var_indices->first, var_indices->second, itime);
      if (OGS::NetCDF::readNetCDF(var_path.c_str(), vel_vars, tmp) !=
          NETCDF_OK) {
        vtkErrorMacro("Cannot read variable <" << this->field
                                               << "> in NetCDF! Aborting!");
        delete current_simulation;
        return 0;
      }

      const auto current_mesh = current_simulation->loaded_mesh();
      // Projecting the velocity field to the UVW grid
      arrayTemp = OGS::field::UVW2T(
          tmp, current_mesh->e1(), current_mesh->e2(), current_mesh->e3(),
          current_mesh->nlon_cells(), current_mesh->nlat_cells(),
          current_mesh->nlev_cells());
    } else {
      const auto var_vname = current_simulation->var_vname(var_indices->first,
                                                           var_indices->second);
      arrayTemp.set_dim(npoints, 1);
      const auto var_path = current_simulation->var_path(
          var_indices->first, var_indices->second, itime);
      if (OGS::NetCDF::readNetCDF(var_path.c_str(), var_vname.c_str(),
                                  arrayTemp) != NETCDF_OK) {
        vtkErrorMacro("Cannot read variable <" << this->field << " ("
                                               << var_vname
                                               << ")> in NetCDF! Aborting!");
        delete current_simulation;
        return 0;
      }
    }

    delete current_simulation;

    if (arrayTemp.get_m() > 1) {
      double mod = 0.;
      for (int jj = 0; jj < arrayTemp.get_m(); ++jj)
        mod += (cList > 0) ? arrayTemp[cList][0] * arrayTemp[cList][0] : 0.;
      spag_array[inst][0] = std::sqrt(mod);
    } else {
      spag_array[inst][0] = (cList > 0) ? arrayTemp[cList][0] : 0.;
    }

    arrayTemp.clear();
    this->UpdateProgress(0.1 + 0.7 / (ii_end - ii_start) * (inst - ii_start));
  }

/* REDUCTION PHASE

        If run in parallel with more than one rank,
        the reduction of the fields is done here.

*/
#ifdef PARAVIEW_USE_MPI

  // Only reduce if we have more than 1 process
  if (this->nProcs > 1) {
    field::Field<FLDARRAY> spag_array_tmp(this->instants.size(), 1, 0.);

    // Reduce the Spaghetti array
    this->Controller->Reduce(spag_array.data(), spag_array_tmp.data(),
                             spag_array.get_sz(),
                             vtkCommunicator::StandardOperations::SUM_OP, 0);

    spag_array.set_val(spag_array_tmp.data());
    spag_array_tmp.clear();
  }

#endif

  this->UpdateProgress(0.9);

  /* FINALIZE

          Finalization, the master process stores the arrays
          inside the output.

  */
  // Build output table
  if (this->procId == 0) {
    // Auxiliar field for columns
    OGS::field::Field<FLDARRAY> column(1, 1, 0.);
    VTKARRAY *vtkColumn;

    // Write the time instants
    for (int ii = 0; ii < this->instants.size(); ++ii) {
      int itime = this->instants[ii];
      column[0][0] = spag_array[ii][0];
      vtkColumn = VTK::createVTKfromField<VTKARRAY, FLDARRAY>(
          this->TL[itime].as_string("%Y%m%d-%H:%M:%S").c_str(), column);
      output->AddColumn(vtkColumn);
      vtkColumn->Delete();
    }
  }

  this->UpdateProgress(1.);
  return 1;
}

//----------------------------------------------------------------------------
int OGSSpaghetti::AveragesIterationAlgorithm(int ntsteps,
                                             vtkInformation *request,
                                             vtkDataSet *input,
                                             vtkDataSet *source,
                                             vtkTable *output) {
// There is no need to parallelize this algorithm
#ifdef PARAVIEW_USE_MPI
  if (this->procId > 0) return 1;
#endif

  // Check inputs
  if (!source->IsA("vtkRectilinearGrid")) {
    vtkErrorMacro(
        "A rectilinear grid is needed for this method! Cannot continue...");
    return 0;
  }

  if (source->GetPointData()->GetNumberOfArrays() > 0) {
    vtkErrorMacro("Data contains point arrays! Cannot continue...");
    return 0;
  }

  // Recover Metadata array (depth factor)
  vtkStringArray *vtkmetadata = vtkStringArray::SafeDownCast(
      source->GetFieldData()->GetAbstractArray("Metadata"));

  // This section is only executed once, to populate the xyz and
  // cId2zId arrays. Successive iterations should not execute.
  // This section is included here since RequestInformation gives
  // troubles when restarting.
  if (this->xyz.isempty() || this->isReqInfo) {
    this->isReqInfo = false;

    this->dfact =
        (vtkmetadata != nullptr) ? std::stod(vtkmetadata->GetValue(2)) : 1000.;

    if (vtkmetadata == nullptr)
      vtkWarningMacro(
          "Field array Metadata not found! Depth factor set to 1000. "
          "automatically.");

    // Recover cell or point coordinates
    this->xyz = getVTKCellCenters(source, this->dfact);

    // Up to this point we have the cell centers or point coordinates correctly
    // stored under "xyz". Now we shall find the number of unique z coordinates
    // or, depending on the user input, the coordinates of each depth level, as
    // well as its mesh connectivity (cId2zId).
    this->cId2zId = OGS::field::countDepthLevels(this->xyz, this->zcoords,
                                                 this->epsi, false);
  }
  this->UpdateProgress(.1);

  // Recover the basins and coasts mask
  // and add them to the output
  int nbasins = 21, ncoasts = 3, nStat = 9;
  VTKMASK *vtkMask;
  OGS::field::Field<FLDMASK> bmask, cmask;
  // Basins mask
  vtkMask =
      VTKMASK::SafeDownCast(source->GetCellData()->GetArray(this->bmask_field));
  bmask = VTK::createFieldfromVTK<VTKMASK, FLDMASK>(vtkMask);
  // Coast mask
  vtkMask =
      VTKMASK::SafeDownCast(source->GetCellData()->GetArray(this->cmask_field));
  cmask = VTK::createFieldfromVTK<VTKMASK, FLDMASK>(vtkMask);

  // Load the mesh weights
  VTKARRAY *vtke3 = NULL;
  vtke3 = VTKARRAY::SafeDownCast(source->GetCellData()->GetArray("e3"));

  if (vtke3 == NULL) {
    vtkErrorMacro("Mesh weights (e1, e2 and e3) need to be loaded to proceed!");
    return 0;
  }

  // Convert to field arrays
  OGS::field::Field<FLDARRAY> e3;
  e3 = VTK::createFieldfromVTK<VTKARRAY, FLDARRAY>(vtke3);

  // Load variable stat profile
  OGS::field::Field<FLDARRAY> statArray(
      ntsteps * nbasins * ncoasts * this->zcoords.size() * nStat, 1, 0.);
  std::string filename = std::string(this->FolderName) + std::string("/") +
                         std::string(this->field) + std::string(".nc");

  if (OGS::NetCDF::readNetCDF(filename.c_str(), this->field, statArray) !=
      NETCDF_OK) {
    // If file cannot be read or variable does not exist
    vtkErrorMacro("File <" << filename.c_str() << "> or variable <"
                           << this->field << "> cannot be read!");
    return 0;
  }
  // Retrieve cellIds from interpolating line
  std::vector<int> cellIds;
  ComputeCellIds(cellIds, input, source);

  // Compute the mesh dependent position that does not depend on time
  // we only need to do this once
  std::vector<int> pos1;
  pos1.resize(cellIds.size(), 0);
  double weight_sum = 0.;
  for (int cc = 0; cc < cellIds.size(); ++cc) {
    // Recover cell id
    int cellId = cellIds[cc];
    // Depth index
    int zId = this->cId2zId[cellId][0];
    // In which basin are we? (we need to loop the basins and find which is
    // true)
    int bId = -1;
    bool isbasin = false;
    for (bId = 0; bId < bmask.get_m(); ++bId) {
      if (bmask[cellId][bId]) {
        isbasin = true;
        break;
      }
    }
    // In which coast are we?
    int cId = this->per_coast ? cmask[cellId][0] - 1 : 2;
    // Compute the position
    int p = ncoasts * this->zcoords.size() * nStat * bId +
            this->zcoords.size() * nStat * cId + nStat * zId + this->sId;
    p = (bId > 0 || cId > 0) ? p : -1;
    // Store the position and the sum of weights
    pos1[cc] = p;
    weight_sum += e3[cellId][0];
  }

  this->UpdateProgress(.25);

  // For each time instant, loop on the cell list and load the data into
  // the table
  OGS::field::Field<FLDARRAY> column(1, 1);
  VTKARRAY *vtkColumn;
  for (int ii = 0; ii < this->instants.size(); ii += 1) {
    int itime = this->instants[ii];
    column[0][0] = 0.;
    for (int cc = 0; cc < cellIds.size(); ++cc) {
      // Recover ids
      int p = pos1[cc];
      int cellId = cellIds[cc];
      // Compute position
      int pos = nbasins * ncoasts * this->zcoords.size() * nStat * itime + p;
      // Retrieve value from array
      column[0][0] += (p > 0) ? e3[cellId][0] * statArray[pos][0] : 0.;
    }
    column[0][0] /= weight_sum;
    vtkColumn = VTK::createVTKfromField<VTKARRAY, FLDARRAY>(
        this->TL[itime].as_string("%Y%m%d-%H:%M:%S").c_str(), column);
    output->AddColumn(vtkColumn);
    vtkColumn->Delete();

    this->UpdateProgress(0.25 + 0.75 / (this->instants.size() - 0) * (ii - 0));
  }

  this->UpdateProgress(1.);
  return 1;
}

//----------------------------------------------------------------------------
void OGSSpaghetti::SetStartTI(const char *tstep) {
  // Set the start point of the TimeInterval

  OGS::Time::TimeObject TO(tstep, "%Y%m%d-%H:%M:%S");
  this->TI.set_start_time(TO);
  this->Modified();
}

void OGSSpaghetti::SetEndTI(const char *tstep) {
  // Set the end point of the TimeInterval

  OGS::Time::TimeObject TO(tstep, "%Y%m%d-%H:%M:%S");
  this->TI.set_end_time(TO);
  this->Modified();
}

//----------------------------------------------------------------------------
vtkStringArray *OGSSpaghetti::GetTimeValues() {
  // Recover time values as a string from the TimeList

  vtkStringArray *TimeValues;
  TimeValues = VTK::createVTKstrf("TimeValues", this->TL.len(), NULL);

  // Loop the TimeList to recover the instants
  for (int ii = 0; ii < this->TL.len(); ++ii)
    TimeValues->SetValue(ii, this->TL[ii].as_string("%Y%m%d-%H:%M:%S"));

  return TimeValues;
}

//----------------------------------------------------------------------------
void OGSSpaghetti::SetSourceConnection(vtkAlgorithmOutput *algOutput) {
  this->SetInputConnection(1, algOutput);
}

void OGSSpaghetti::SetSourceData(vtkDataObject *input) {
  this->SetInputData(1, input);
}

vtkDataObject *OGSSpaghetti::GetSource() {
  if (this->GetNumberOfInputConnections(1) < 1) return nullptr;

  return this->GetExecutive()->GetInputData(1, 0);
}