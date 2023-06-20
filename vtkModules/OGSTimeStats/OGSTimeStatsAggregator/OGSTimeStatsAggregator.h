/*=========================================================================

  Program:   OGSTimeAggregator
  Module:    OGSTimeStatsAggregator.h

  Copyright (c) 2020 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSTimeAggregator_h
#define vtkOGSTimeAggregator_h

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

class VTK_EXPORT OGSTimeStatsAggregator : public vtkDataSetAlgorithm {
 public:
  static OGSTimeStatsAggregator *New();
  vtkTypeMacro(OGSTimeStatsAggregator, vtkDataSetAlgorithm);

  // Description:
  // Selection of the requestor
  vtkGetMacro(ReqType, int);
  vtkSetMacro(ReqType, int);

  // Description:
  // Selection of the weekday for the weekly list
  vtkGetMacro(weekday, int);
  vtkSetMacro(weekday, int);

  // Description:
  // Selection of the algorithm
  vtkGetMacro(use_files, int);
  vtkSetMacro(use_files, int);
  vtkBooleanMacro(use_files, int);

#ifdef PARAVIEW_USE_MPI
  vtkGetObjectMacro(Controller, vtkMultiProcessController);
  virtual void SetController(vtkMultiProcessController *);
#endif

 protected:
  OGSTimeStatsAggregator();
  ~OGSTimeStatsAggregator();

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
  OGSTimeStatsAggregator(const OGSTimeStatsAggregator &) = delete;
  void operator=(const OGSTimeStatsAggregator &) = delete;

  int PipelineIterationAlgorithm(vtkInformation *, vtkDataSet *, vtkDataSet *);
  int FileIterationAlgorithm(vtkInformation *, vtkDataSet *, vtkDataSet *);

  OGS::Time::TimeList TL;       // TimeList containing all the instants
  OGS::Time::REQ_LIST ReqList;  // Requestor list


  std::vector<int> instants;    // Instant ID to loop
  std::vector<double> weights;  // Weights for the instants

  double sum_weights;
  bool TL_computed, use_files;

  int procId, nProcs, weekday, ReqType, CurrentTimeIndex;
};

#endif
