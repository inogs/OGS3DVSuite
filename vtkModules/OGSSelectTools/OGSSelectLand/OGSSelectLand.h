// -*- c++ -*-
/*=========================================================================

  Program:   OGSSelectTools
  Module:    OGSSelectLand.h

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef OGSSelectLand_h
#define OGSSelectLand_h

#include "vtkDataSet.h"
#include "vtkThreshold.h"

#ifdef PARAVIEW_USE_MPI
class vtkMultiProcessController;
#endif

class VTKFILTERSCORE_EXPORT OGSSelectLand : public vtkThreshold {
 public:
  static OGSSelectLand* New();
  vtkTypeMacro(OGSSelectLand, vtkThreshold);

  OGSSelectLand(const OGSSelectLand&) = delete;
  void operator=(const OGSSelectLand&) = delete;

  void PrintSelf(ostream& os, vtkIndent indent) override;

  // Description:
  // Get the name of the mask field to operate
  vtkSetStringMacro(mask_field);

#ifdef PARAVIEW_USE_MPI
  // Description:
  // Set the controller use in compositing (set to
  // the global controller by default)
  // If not using the default, this must be called before any
  // other methods.
  virtual void SetController(vtkMultiProcessController* controller);
#endif

 protected:
  OGSSelectLand();
  ~OGSSelectLand() override;

  int RequestData(vtkInformation*, vtkInformationVector**,
                  vtkInformationVector*) override;

#ifdef PARAVIEW_USE_MPI
  vtkMultiProcessController* Controller;
#endif

 private:
  int procId, nProcs;
  const char* mask_field;
};

#endif
