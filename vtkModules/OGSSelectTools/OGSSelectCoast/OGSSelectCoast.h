// -*- c++ -*-
/*=========================================================================

  Program:   OGSSelectTools
  Module:    OGSSelectCoast.h

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef OGSSelectCoast_h
#define OGSSelectCoast_h

#include "vtkDataSet.h"
#include "vtkThreshold.h"


class vtkDataArraySelection;

#ifdef PARAVIEW_USE_MPI
class vtkMultiProcessController;
#endif

class VTKFILTERSCORE_EXPORT OGSSelectCoast : public vtkThreshold {
 public:
  static OGSSelectCoast* New();
  vtkTypeMacro(OGSSelectCoast, vtkThreshold);

  OGSSelectCoast(const OGSSelectCoast&) = delete;
  void operator=(const OGSSelectCoast&) = delete;

  void PrintSelf(ostream& os, vtkIndent indent) override;

  // Description:
  // Get the name of the mask field to operate
  vtkSetStringMacro(mask_field);

  // Description:
  // The following methods allow selective selection of the basins.
  int GetNumberOfCoastsArrays();
  const char* GetCoastsArrayName(int index);
  int GetCoastsArrayIndex(const char* name);
  int GetCoastsArrayStatus(const char* name);
  void SetCoastsArrayStatus(const char* name, int status);
  void DisableAllCoastsArrays();
  void EnableAllCoastsArrays();

#ifdef PARAVIEW_USE_MPI
  // Description:
  // Set the controller use in compositing (set to
  // the global controller by default)
  // If not using the default, this must be called before any
  // other methods.
  virtual void SetController(vtkMultiProcessController* controller);
#endif

 protected:
  OGSSelectCoast();
  ~OGSSelectCoast() override;

  int RequestData(vtkInformation*, vtkInformationVector**,
                  vtkInformationVector*) override;

  vtkDataArraySelection* CoastsDataArraySelection;

#ifdef PARAVIEW_USE_MPI
  vtkMultiProcessController* Controller;
#endif

 private:
  int procId, nProcs;
  const char* mask_field;
};

#endif
