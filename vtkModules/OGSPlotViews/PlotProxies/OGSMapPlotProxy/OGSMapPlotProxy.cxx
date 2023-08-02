/*=========================================================================

  Program:   OGSPlotViews
  Module:    vtkSMOGSMapPlotProxy.cxx

  Copyright (c) 2019 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "OGSMapPlotProxy.h"

#include "OGSMapPlot.h"
#include "vtkClientServerStream.h"
#include "vtkDataArray.h"
#include "vtkMultiProcessController.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkProcessModule.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkSMViewProxyInteractorHelper.h"

vtkStandardNewMacro(OGSMapPlotProxy);

//----------------------------------------------------------------------------
OGSMapPlotProxy::OGSMapPlotProxy() {}

//----------------------------------------------------------------------------
OGSMapPlotProxy::~OGSMapPlotProxy() {}

//----------------------------------------------------------------------------
vtkImageData* OGSMapPlotProxy::CaptureWindowInternal(int magX, int magY) {
  OGSMapPlot* pv = OGSMapPlot::SafeDownCast(this->GetClientSideObject());
  if (pv) pv->SetMagnification(magX, magY);

  vtkImageData* image = this->Superclass::CaptureWindowInternal(magX, magY);
  if (pv) pv->SetMagnification(1, 1);

  return image;
}
