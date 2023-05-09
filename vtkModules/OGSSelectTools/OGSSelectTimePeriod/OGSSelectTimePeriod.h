/*=========================================================================

  Program:   OGSSelectTimePeriod
  Module:    vtkOGSSelectTimePeriod.h

  Copyright (c) 2020 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef OGSSelectTimePeriod_h
#define OGSSelectTimePeriod_h

#include <string>
#include <vector>

#include "OGS/TimeInterval.h"
#include "OGS/TimeList.h"
#include "vtkDataSetAlgorithm.h"
#include "vtkStringArray.h"

class vtkStringArray;

#ifdef PARAVIEW_USE_MPI
class vtkMultiProcessController;
#endif

class VTKCOMMONEXECUTIONMODEL_EXPORT OGSSelectTimePeriod
    : public vtkDataSetAlgorithm {
 public:
  static OGSSelectTimePeriod *New();
  vtkTypeMacro(OGSSelectTimePeriod, vtkDataSetAlgorithm);

  OGSSelectTimePeriod(const OGSSelectTimePeriod &) = delete;
  void operator=(const OGSSelectTimePeriod &) = delete;

  void SetStartTI(const char *tstep);
  void SetEndTI(const char *tstep);

  vtkStringArray *GetTimeValues();

#ifdef PARAVIEW_USE_MPI
  vtkGetObjectMacro(Controller, vtkMultiProcessController);
  virtual void SetController(vtkMultiProcessController *);
#endif

 protected:
  OGSSelectTimePeriod();
  ~OGSSelectTimePeriod() override = default;

  int RequestInformation(vtkInformation *, vtkInformationVector **,
                         vtkInformationVector *) override;
  int RequestData(vtkInformation *, vtkInformationVector **,
                  vtkInformationVector *) override;

#ifdef PARAVIEW_USE_MPI
  vtkMultiProcessController *Controller;
#endif

 private:
  OGS::Time::TimeInterval TI;   // TimeInterval for the generic requestor
  OGS::Time::TimeList TL;       // TimeList containing all the instants

  std::vector<int> instants;    // Instant ID to loop
  std::vector<double> weights;  // Weights for the instants

  bool TL_computed;

  int procId, nProcs, CurrentTimeIndex;
};

#endif
