/*=========================================================================

  Program:   OGSTimeAverage
  Module:    OGSTimeStatsAverage.h

  Copyright (c) 2020 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSTimeAverage_h
#define vtkOGSTimeAverage_h

#include <string>
#include <vector>

#include "OGS/ParaviewImports.h"  // For PARAVIEW_USE_MPI
#include "OGS/TimeInterval.h"
#include "OGS/TimeList.h"
#include "vtkDataSetAlgorithm.h"
#include "vtkStringArray.h"

class vtkStringArray;

#ifdef PARAVIEW_USE_MPI
class vtkMultiProcessController;
#endif

class VTK_EXPORT OGSTimeStatsAverage : public vtkDataSetAlgorithm {
 public:
  static OGSTimeStatsAverage *New();
  vtkTypeMacro(OGSTimeStatsAverage, vtkDataSetAlgorithm);

  // Description:
  // Selection of the algorithm
  vtkGetMacro(use_files, int);
  vtkSetMacro(use_files, int);
  vtkBooleanMacro(use_files, int);

  void SetStartTI(const char *tstep);
  void SetEndTI(const char *tstep);

  vtkStringArray *GetTimeValues();

#ifdef PARAVIEW_USE_MPI
  vtkGetObjectMacro(Controller, vtkMultiProcessController);
  virtual void SetController(vtkMultiProcessController *);
#endif

 protected:
  OGSTimeStatsAverage();
  ~OGSTimeStatsAverage();

  int RequestInformation(vtkInformation *, vtkInformationVector **,
                         vtkInformationVector *) override;
  int RequestUpdateExtent(vtkInformation *, vtkInformationVector **,
                          vtkInformationVector *) override;
  int RequestData(vtkInformation *, vtkInformationVector **,
                  vtkInformationVector *) override;

#ifdef PARAVIEW_USE_MPI
  vtkMultiProcessController *Controller;
#endif

 private:
  OGSTimeStatsAverage(const OGSTimeStatsAverage &) = delete;
  void operator=(const OGSTimeStatsAverage &) = delete;

  int PipelineIterationAlgorithm(vtkInformation *, vtkDataSet *, vtkDataSet *);
  int FileIterationAlgorithm(vtkInformation *, vtkDataSet *, vtkDataSet *);

  OGS::Time::TimeInterval TI;   // TimeInterval for the generic requestor
  OGS::Time::TimeList TL;       // TimeList containing all the instants

  std::vector<int> instants;    // Instant ID to loop
  std::vector<double> weights;  // Weights for the instants

  double sum_weights;
  bool TL_computed, use_files;

  int procId, nProcs, CurrentTimeIndex;
};

#endif
