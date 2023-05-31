// -*- c++ -*-
/*=========================================================================

  Program:   OGSDepthProfile
  Module:    OGSDepthProfile.h

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef OGSDepthProfile_h
#define OGSDepthProfile_h

#include "OGS/ParaviewImports.h"
#include "vtkDataSetAlgorithm.h"
#include "vtkDataSetAttributes.h"

class vtkAbstractCellLocator;

#ifdef PARAVIEW_USE_MPI
class vtkMultiProcessController;
#endif

class VTKCOMMONEXECUTIONMODEL_EXPORT OGSDepthProfile
    : public vtkDataSetAlgorithm {
 public:
  static OGSDepthProfile *New();
  vtkTypeMacro(OGSDepthProfile, vtkDataSetAlgorithm);

  OGSDepthProfile(const OGSDepthProfile &) = delete;
  void operator=(const OGSDepthProfile &) = delete;

  // Description:
  // Specify the data set that will be probed at the input points.
  // The Input gives the geometry (the points and cells) for the output,
  // while the Source is probed (interpolated) to generate the scalars,
  // vectors, etc. for the output points based on the point locations.
  void SetSourceData(vtkDataObject *source);
  vtkDataObject *GetSource();
  void SetSourceConnection(vtkAlgorithmOutput *algOutput);

  // Description
  // Set/Get the prototype cell locator to use for probing the source dataset.
  // By default, vtkStaticCellLocator will be used.
  virtual void SetCellLocatorPrototype(vtkAbstractCellLocator *);
  vtkGetObjectMacro(CellLocatorPrototype, vtkAbstractCellLocator);

#ifdef PARAVIEW_USE_MPI
  // Description:
  // Set the controller use in compositing (set to
  // the global controller by default)
  // If not using the default, this must be called before any
  // other methods.
  virtual void SetController(vtkMultiProcessController *controller);
#endif

 protected:
  OGSDepthProfile();
  ~OGSDepthProfile() override;

  int RequestInformation(vtkInformation *, vtkInformationVector **,
                         vtkInformationVector *) override;
  int RequestUpdateExtent(vtkInformation *, vtkInformationVector **,
                          vtkInformationVector *) override;
  int RequestData(vtkInformation *, vtkInformationVector **,
                  vtkInformationVector *) override;

  void Initialize(vtkDataSet *input, vtkDataSet *source, vtkDataSet *output);
  void Interpolate(vtkDataSet *input, vtkDataSet *source, vtkDataSet *output);

#ifdef PARAVIEW_USE_MPI
  vtkMultiProcessController *Controller;
#endif

 private:
  int procId, nProcs;

  vtkAbstractCellLocator *CellLocatorPrototype;

  vtkDataSetAttributes::FieldList *CellList;
  vtkDataSetAttributes::FieldList *PointList;
};

#endif
