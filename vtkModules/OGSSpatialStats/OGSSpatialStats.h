// -*- c++ -*-
/*=========================================================================

  Program:   OGSSpatialStats
  Module:    OGSSpatialStats.h

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef OGSSpatialStats_h
#define OGSSpatialStats_h

#include <vector>

#include "OGS/ParaviewImports.h"
#include "OGS/V3.h"
#include "OGS/field.h"
#include "vtkDataArraySelection.h"
#include "vtkDataSet.h"
#include "vtkDataSetAlgorithm.h"

#ifdef PARAVIEW_USE_MPI
class vtkMultiProcessController;
#endif

class VTKCOMMONEXECUTIONMODEL_EXPORT OGSSpatialStats
    : public vtkDataSetAlgorithm {
 public:
  static OGSSpatialStats* New();
  vtkTypeMacro(OGSSpatialStats, vtkDataSetAlgorithm);

  // Description:
  // Get the epsilon
  vtkSetMacro(epsi, double);
  vtkGetMacro(epsi, double);

  // Description:
  // Flag for changing mesh
  vtkGetMacro(changing_mesh, int);
  vtkSetMacro(changing_mesh, int);
  vtkBooleanMacro(changing_mesh, int);

  // Description:
  // Decide to use the volume for statistics
  vtkGetMacro(useVolume, int);
  vtkSetMacro(useVolume, int);
  vtkBooleanMacro(useVolume, int);

  // Description:
  // Number of user inputed depths
  void SetNumberOfDepthLevels(int n);
  void SetDepthLevels(int i, double value);

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
  OGSSpatialStats();
  ~OGSSpatialStats() override;

  int RequestData(vtkInformation*, vtkInformationVector**,
                  vtkInformationVector*) override;
  int RequestInformation(vtkInformation*, vtkInformationVector**,
                         vtkInformationVector*) override;

  vtkDataArraySelection* StatDataArraySelection;

#ifdef PARAVIEW_USE_MPI
  vtkMultiProcessController* Controller;
#endif

 private:
  OGSSpatialStats(const OGSSpatialStats&) = delete;
  void operator=(const OGSSpatialStats&) = delete;

  // Tolerance for finding the depth levels
  double epsi;

  // Depth levels and number of depth levels
  int ndepths, procId, nProcs;
  std::vector<double> zcoords;

  // Auxiliar variables worth conserving
  // they are read once in the RequestInformation and used in the
  // RequestData
  bool iscelld;                    // Whether we have cell or point data
  OGS::V3::V3v xyz;                // Stores cell/point coordinates
  OGS::field::Field<int> cId2zId;  // Cell to depth level connectivity

  bool isReqInfo, changing_mesh;  // Set true when request information
  int useVolume;  // Use the volume as weight instead of the area
};

#endif
