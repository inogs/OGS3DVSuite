// -*- c++ -*-
/*=========================================================================

  Program:   OGSSelectTools
  Module:    OGSSelectPolygon.h

  Copyright (c) 2019 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef OGSSelectPolygon_h
#define OGSSelectPolygon_h

#include "vtkDataSet.h"
#include "vtkThreshold.h"

#include <string>
#include "OGS/geometry.h"

#ifdef PARAVIEW_USE_MPI
  class vtkMultiProcessController;
#endif

class VTKFILTERSCORE_EXPORT OGSSelectPolygon : public vtkThreshold {
public:
  static OGSSelectPolygon * New();
  vtkTypeMacro(OGSSelectPolygon, vtkThreshold);

  OGSSelectPolygon(const OGSSelectPolygon &) = delete;
  void operator=(const OGSSelectPolygon &) = delete;

  // Description:
  // Obtain the polygon points from the textbox.
  void GetPolygon(const char *);

  // Description:
  vtkGetMacro(Invert, bool);
  vtkSetMacro(Invert, bool);

protected:
 OGSSelectPolygon();
  ~OGSSelectPolygon() override;

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) override;

private:

  bool Invert;
  int nProcs;
  double dfact;

  std::string projName;
  OGS::Geom::Polygon<double> poly;
};

#endif
