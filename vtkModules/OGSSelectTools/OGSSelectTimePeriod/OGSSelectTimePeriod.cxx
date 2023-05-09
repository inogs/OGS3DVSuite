/*=========================================================================

  Program:   OGSSelectTimePeriod
  Module:    vtkOGSSelectTimePeriod.cxx

  Copyright (c) 2020 Arnau Miro, OGS
  All rights reserved.

         This software is distributed WITHOUT ANY WARRANTY; without even
         the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
         PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "OGSSelectTimePeriod.h"

#include <chrono>
#include <ctime>
#include <string>
#include <vector>

#include "vtkCellData.h"
#include "vtkDataSet.h"
#include "vtkFieldData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStringArray.h"

vtkStandardNewMacro(OGSSelectTimePeriod)

#ifdef PARAVIEW_USE_MPI
#include "vtkMultiProcessController.h"
    vtkCxxSetObjectMacro(vtkOGSSelectTimePeriod, Controller,
                         vtkMultiProcessController);
#endif

//----------------------------------------------------------------------------
#include "OGS/TimeList.h"
#include "OGS/TimeObject.h"
#include "OGS/macros.h"
#include "OGS/vtkFields.h"

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
OGSSelectTimePeriod::OGSSelectTimePeriod() {
  this->TL_computed = false;
  this->CurrentTimeIndex = 0;

  this->procId = 0;
  this->nProcs = 0;

#ifdef PARAVIEW_USE_MPI
  this->Controller = nullptr;
  this->SetController(vtkMultiProcessController::GetGlobalController());
#endif
}


//----------------------------------------------------------------------------
int OGSSelectTimePeriod::RequestInformation(
    vtkInformation *vtkNotUsed(request),
    vtkInformationVector **vtkNotUsed(inputVector),
    vtkInformationVector *outputVector) {
  // Grab the input
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

// Handle the parallel controller
#ifdef PARAVIEW_USE_MPI
  if (this->Controller->GetNumberOfProcesses() > 1) {
    this->procId = this->Controller->GetLocalProcessId();
    this->nProcs = this->Controller->GetNumberOfProcesses();
  }
#endif

  // Compute the TimeList if it hasn't been computed
  if (!this->TL_computed) {
    BuildTimeList(this->TL, outInfo);
    this->TL_computed = true;
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

  // Set the output timesteps
  outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
  outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_RANGE());

  auto *timeSteps = new double[this->instants.size()];
  for (int ii = 0; ii < this->instants.size(); ii++) {
    struct tm tm = this->TL[this->instants[ii]].get_tm();
    timeSteps[ii] = difftime(mktime(&tm), 0);
  }
  outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), &timeSteps[0],
               static_cast<int>(this->instants.size()));

  // Set up the time range
  double timeRange[2] = {timeSteps[0], timeSteps[this->instants.size() - 1]};
  outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), timeRange, 2);

  delete[] timeSteps;
  this->CurrentTimeIndex == 1;

  return 1;
}

//----------------------------------------------------------------------------
int OGSSelectTimePeriod::RequestData(vtkInformation *request,
                                     vtkInformationVector **inputVector,
                                     vtkInformationVector *outputVector) {
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input and output
  vtkDataSet *input =
      vtkDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkDataSet *output =
      vtkDataSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  if (this->CurrentTimeIndex == 0) {
    // Master computes the instants to loop
    if (this->procId == 0) {
      // Compute the TimeList if it hasn't been computed
      if (!this->TL_computed) {
        BuildTimeList(this->TL, outInfo);
        this->TL_computed = true;
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

      // Set the output timesteps
      outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
      outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_RANGE());

      auto *timeSteps = new double[this->instants.size()];
      for (int ii = 0; ii < this->instants.size(); ii++) {
        struct tm tm = this->TL[this->instants[ii]].get_tm();
        timeSteps[ii] = difftime(mktime(&tm), 0);
      }
      outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),
                   &timeSteps[0], static_cast<int>(this->instants.size()));

      // Set up the time range
      double timeRange[2] = {timeSteps[0],
                             timeSteps[this->instants.size() - 1]};
      outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), timeRange,
                   2);

      this->CurrentTimeIndex == 1;
      delete[] timeSteps;
    }
  }

  // Output the requested instant
  output->ShallowCopy(input);

  this->UpdateProgress(1.);
  return 1;
}

//----------------------------------------------------------------------------
void OGSSelectTimePeriod::SetStartTI(const char *tstep) {
  // Set the start point of the TimeInterval

  OGS::Time::TimeObject TO(tstep, "%Y%m%d-%H:%M:%S");
  this->TI.set_start_time(TO);
  this->Modified();
}

void OGSSelectTimePeriod::SetEndTI(const char *tstep) {
  // Set the end point of the TimeInterval

  OGS::Time::TimeObject TO(tstep, "%Y%m%d-%H:%M:%S");
  this->TI.set_end_time(TO);
  this->Modified();
}

//----------------------------------------------------------------------------
vtkStringArray *OGSSelectTimePeriod::GetTimeValues() {
  // Recover time values as a string from the TimeList

  vtkStringArray *TimeValues;
  TimeValues = VTK::createVTKstrf("TimeValues", this->TL.len());

  // Loop the TimeList to recover the instants
  for (int ii = 0; ii < this->TL.len(); ++ii)
    TimeValues->SetValue(ii, this->TL[ii].as_string("%Y%m%d-%H:%M:%S"));

  return TimeValues;
}
