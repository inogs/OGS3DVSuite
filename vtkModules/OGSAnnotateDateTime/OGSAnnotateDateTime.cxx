/*=========================================================================

  Program:   OGSAnnotateDateTime
  Module:    OGSAnnotateDateTime.cxx

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "OGSAnnotateDateTime.h"

#include <chrono>
#include <ctime>

#include "vtkDataSet.h"
#include "vtkFieldData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStringArray.h"

vtkStandardNewMacro(OGSAnnotateDateTime)

    //----------------------------------------------------------------------------
    OGSAnnotateDateTime::OGSAnnotateDateTime() {
  this->TimeFormat = nullptr;
  this->useMetadata = 0;
  this->procId = 0;
}

//----------------------------------------------------------------------------
OGSAnnotateDateTime::~OGSAnnotateDateTime() { this->TimeFormat = nullptr; }

//----------------------------------------------------------------------------
int OGSAnnotateDateTime::RequestData(vtkInformation* request,
                                     vtkInformationVector** inputVector,
                                     vtkInformationVector* outputVector) {
  // Recover the Date string
  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);

  // Decide how to write the time, whether to use metadata or to use a
  // conversion from the timestep to a string format
  struct tm current_tm = {0};
  std::ostringstream buf;

  // Ensure that we are working in UTC (but probably this does not make any
  // difference)
  auto current_locale_time = setlocale(LC_TIME, nullptr);
  setlocale(LC_TIME, "C");

  if (this->useMetadata) {
    vtkDataSet* input =
        vtkDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

    // Recover metadata
    vtkStringArray* vtkmetadata = vtkStringArray::SafeDownCast(
        input->GetFieldData()->GetAbstractArray("Metadata"));

    // If successful, modify the entry by parsing the string
    // using the time library
    if (vtkmetadata) {
      strptime((vtkmetadata->GetValue(0)).c_str(), "%Y%m%d-%H:%M:%S",
               &current_tm);
      buf << std::put_time(&current_tm, this->TimeFormat);
      this->Superclass::SetFormat(buf.str().c_str());
    }
  } else {
    // Recover the current timestep
    double requestedTimeValue =
        inInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());

    // Convert to struct tm
    time_t time = std::chrono::system_clock::to_time_t(
        std::chrono::system_clock::time_point(
            std::chrono::duration_cast<std::chrono::seconds>(
                std::chrono::duration<double>(requestedTimeValue))));
    current_tm = *gmtime(&time);

    // Format and display
    buf << std::put_time(&current_tm, this->TimeFormat);
    this->Superclass::SetFormat(buf.str().c_str());
  }

  // Run RequestData from vtkTimeToTextConvertor
  return this->Superclass::RequestData(request, inputVector, outputVector);
}
