// -*- c++ -*-
/*=========================================================================

  Program:   OGSSpatialStats
  Module:    OGSSpatialStatsFromFile.h

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef OGSSpatialStatsFromFile_h
#define OGSSpatialStatsFromFile_h

#include "OGS/ParaviewImports.h"
#include "vtkDataArraySelection.h"
#include "vtkDataSet.h"
#include "vtkDataSetAlgorithm.h"
#include "vtkRectilinearGridAlgorithm.h"

#ifdef PARAVIEW_USE_MPI
class vtkMultiProcessController;
#endif

class VTKCOMMONEXECUTIONMODEL_EXPORT OGSSpatialStatsFromFile
    : public vtkRectilinearGridAlgorithm {
 public:
  static OGSSpatialStatsFromFile* New();
  vtkTypeMacro(OGSSpatialStatsFromFile, vtkRectilinearGridAlgorithm);

  // Description:
  // Get the folder where the stats are
  vtkSetStringMacro(FolderName);

  // Description:
  // Get the name of the mask fields to operate
  vtkSetStringMacro(bmask_field);
  vtkSetStringMacro(cmask_field);
  vtkSetStringMacro(lmask_field);

  // Description:
  // If true, obtains the statistics per basin and per coast
  vtkGetMacro(per_coast, int);
  vtkSetMacro(per_coast, int);
  vtkBooleanMacro(per_coast, int);

  // Description:
  // The following methods allow selective seleccion of aggregate variables.
  int GetNumberOfStatArrays();
  const char* GetStatArrayName(int index);
  int GetStatArrayIndex(const char* name);
  int GetStatArrayStatus(const char* name);
  void SetStatArrayStatus(const char* name, int status);
  void DisableAllStatArrays();
  void EnableAllStatArrays();

#ifdef PARAVIEW_USE_MPI
  // Description:
  // Set the controller use in compositing (set to
  // the global controller by default)
  // If not using the default, this must be called before any
  // other methods.
  virtual void SetController(vtkMultiProcessController* controller);
#endif

 protected:
  OGSSpatialStatsFromFile();
  ~OGSSpatialStatsFromFile() override;

  int RequestData(vtkInformation*, vtkInformationVector**,
                  vtkInformationVector*) override;

  vtkDataArraySelection* StatDataArraySelection;

#ifdef PARAVIEW_USE_MPI
  vtkMultiProcessController* Controller;
#endif

 private:
  OGSSpatialStatsFromFile(const OGSSpatialStatsFromFile&) = delete;
  void operator=(const OGSSpatialStatsFromFile&) = delete;

  char* FolderName;
  char *bmask_field, *cmask_field, *lmask_field;

  int per_coast, procId, nProcs;
};

#endif
