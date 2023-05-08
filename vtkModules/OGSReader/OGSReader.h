// -*- c++ -*-
/*=========================================================================

  Program:   OGSReader
  Module:    OGSReader.h

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef OGSReader_h
#define OGSReader_h

#include <vtkRectilinearGridAlgorithm.h>
#include <vtkStructuredGridAlgorithm.h>

#include <string>

#include "OGS/OGS3DVSuite.h"
#include "vtkDataSet.h"
#include "vtkStreamingDemandDrivenPipeline.h"

class vtkDataArraySelection;
class vtkStringArray;

//----------------------------------------------------------------------------

typedef vtkStructuredGridAlgorithm AlgorithmType;
typedef vtkStructuredGrid MeshType;
// typedef vtkRectilinearGridAlgorithm AlgorithmType;
// typedef vtkRectilinearGrid MeshType;


// Filter class
class VTKCOMMONEXECUTIONMODEL_EXPORT OGSReader : public AlgorithmType {
 public:
  static OGSReader* New();
  vtkTypeMacro(OGSReader, AlgorithmType)

      OGSReader(const OGSReader&) = delete;

  void operator=(const OGSReader&) = delete;

  // Description:
  // Get the name of the master file to read
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  // Description:
  // If false, do not include the meshmask. True by default.
  vtkGetMacro(RMeshMask, int);
  vtkSetMacro(RMeshMask, int);
  vtkBooleanMacro(RMeshMask, int);

  // Description:
  // Lets the user select a multiplier factor for the depth
  vtkGetMacro(DepthScale, double);
  vtkSetMacro(DepthScale, double);

  // Description:
  // The following methods allow selective reading of the mask variables.
  // By default, ALL variables are read, but this can be modified
  // (e.g. from the ParaView GUI).
  int GetNumberOfMaskArrays();
  const char* GetMaskArrayName(int index);
  int GetMaskArrayIndex(const char* name);
  int GetMaskArrayStatus(const char* name);
  void SetMaskArrayStatus(const char* name, int status);
  void DisableAllMaskArrays();
  void EnableAllMaskArrays();

  // Description:
  // The following methods allow selective reading of the physical variables.
  // By default, ALL variables are read, but this can be modified
  // (e.g. from the ParaView GUI).
  int GetNumberOfAvePhysArrays();
  const char* GetAvePhysArrayName(int index);
  int GetAvePhysArrayIndex(const char* name);
  int GetAvePhysArrayStatus(const char* name);
  void SetAvePhysArrayStatus(const char* name, int status);
  void DisableAllAvePhysArrays();
  void EnableAllAvePhysArrays();

  // Description:
  // The following methods allow selective reading of the biogeochemical
  // variables. By default, ALL variables are read, but this can be modified
  // (e.g. from the ParaView GUI).
  int GetNumberOfAveFreqArrays();
  const char* GetAveFreqArrayName(int index);
  int GetAveFreqArrayIndex(const char* name);
  int GetAveFreqArrayStatus(const char* name);
  void SetAveFreqArrayStatus(const char* name, int status);
  void DisableAllAveFreqArrays();
  void EnableAllAveFreqArrays();

  // Description:
  // The following methods allow selective reading of the forcing variables.
  // By default, ALL variables are read, but this can be modified
  // (e.g. from the ParaView GUI).
  int GetNumberOfForcingArrays();
  const char* GetForcingArrayName(int index);
  int GetForcingArrayIndex(const char* name);
  int GetForcingArrayStatus(const char* name);
  void SetForcingArrayStatus(const char* name, int status);
  void DisableAllForcingArrays();
  void EnableAllForcingArrays();

  // Description:
  // The following methods allow selective reading of the general variables.
  // By default, ALL variables are read, but this can be modified
  // (e.g. from the ParaView GUI).
  int GetNumberOfGeneralArrays();
  const char* GetGeneralArrayName(int index);
  int GetGeneralArrayIndex(const char* name);
  int GetGeneralArrayStatus(const char* name);
  void SetGeneralArrayStatus(const char* name, int status);
  void DisableAllGeneralArrays();
  void EnableAllGeneralArrays();

  vtkStringArray* GetProjections();
  void SetProjection(const char* proj);

  vtkDataArraySelection* MaskDataArraySelection;

#ifdef PARAVIEW_USE_MPI
  // Description:
  // Set the controller use in compositing (set to
  // the global controller by default)
  // If not using the default, this must be called before any
  // other methods.
  virtual void SetController(vtkMultiProcessController* controller);
#endif

 protected:
  OGSReader();
  ~OGSReader() override;

  void initialize();

  int RequestInformation(vtkInformation*, vtkInformationVector**,
                         vtkInformationVector*) override;
  int RequestData(vtkInformation*, vtkInformationVector**,
                  vtkInformationVector*) override;

  // Variables
  char* FileName;

  int RMeshMask;
  double DepthScale;

  vtkDataArraySelection* AvePhysDataArraySelection;
  vtkDataArraySelection* AveFreqDataArraySelection;
  vtkDataArraySelection* ForcingDataArraySelection;
  vtkDataArraySelection* GeneralDataArraySelection;

  vtkStringArray* Projections;

#ifdef PARAVIEW_USE_MPI
  vtkMultiProcessController* Controller;
#endif

 private:
  int abort, procId, nProcs, projId;
  MeshType* Mesh;

  std::string projName = "UNDEFINED";

  std::unique_ptr<OGS::Simulation> ogsdata = nullptr;
};

#endif
