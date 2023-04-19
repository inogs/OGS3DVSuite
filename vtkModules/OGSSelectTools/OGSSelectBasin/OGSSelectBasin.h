// -*- c++ -*-
/*=========================================================================

  Program:   OGSSelectTools
  Module:    OGSSelectBasin.h

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef OGSSelectBasin_h
#define OGSSelectBasin_h

#include "vtkDataSet.h"
#include "vtkThreshold.h"

class vtkDataArraySelection;

#ifdef PARAVIEW_USE_MPI
class vtkMultiProcessController;
#endif

class VTKFILTERSCORE_EXPORT OGSSelectBasin : public vtkThreshold {
 public:
  static OGSSelectBasin* New();
  vtkTypeMacro(OGSSelectBasin, vtkThreshold);

  OGSSelectBasin(const OGSSelectBasin&) = delete;
  void operator=(const OGSSelectBasin&) = delete;

  void PrintSelf(ostream& os, vtkIndent indent) override;

  // Description:
  // Get the name of the mask field to operate
  vtkSetStringMacro(mask_field);

  // Description:
  // If false, do not include the meshmask. True by default.
  vtkGetMacro(basins_field, int);
  vtkSetMacro(basins_field, int);
  vtkBooleanMacro(basins_field, int);

  // Description:
  // The following methods allow selective selection of the basins.
  int GetNumberOfBasinsArrays();
  const char* GetBasinsArrayName(int index);
  int GetBasinsArrayIndex(const char* name);
  int GetBasinsArrayStatus(const char* name);
  void SetBasinsArrayStatus(const char* name, int status);
  void DisableAllBasinsArrays();
  void EnableAllBasinsArrays();

#ifdef PARAVIEW_USE_MPI
  // Description:
  // Set the controller use in compositing (set to
  // the global controller by default)
  // If not using the default, this must be called before any
  // other methods.
  virtual void SetController(vtkMultiProcessController* controller);
#endif

 protected:
  OGSSelectBasin();
  ~OGSSelectBasin() override;

  int RequestData(vtkInformation*, vtkInformationVector**,
                  vtkInformationVector*) override;

  vtkDataArraySelection* BasinsDataArraySelection;

#ifdef PARAVIEW_USE_MPI
  vtkMultiProcessController* Controller;
#endif

 private:
  int procId, nProcs, basins_field;
  const char* mask_field;
};

#endif
