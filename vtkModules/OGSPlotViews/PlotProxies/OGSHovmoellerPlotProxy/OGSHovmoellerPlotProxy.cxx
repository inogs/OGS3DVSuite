/*=========================================================================

  Program:   OGSPlotViews
  Module:    vtkSMOGSHovmoellerPlotProxy.cxx

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "OGSHovmoellerPlotProxy.h"

#include "OGSHovmoellerPlot.h"
#include "vtkClientServerStream.h"
#include "vtkDataArray.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkSMViewProxyInteractorHelper.h"

vtkStandardNewMacro(OGSHovmoellerPlotProxy);

//----------------------------------------------------------------------------
OGSHovmoellerPlotProxy::OGSHovmoellerPlotProxy() {}

//----------------------------------------------------------------------------
OGSHovmoellerPlotProxy::~OGSHovmoellerPlotProxy() {}

//----------------------------------------------------------------------------
vtkImageData* OGSHovmoellerPlotProxy::CaptureWindowInternal(int magX,
                                                            int magY) {
  OGSHovmoellerPlot* pv =
      OGSHovmoellerPlot::SafeDownCast(this->GetClientSideObject());
  if (pv) pv->SetMagnification(magX, magY);

  vtkImageData* image = this->Superclass::CaptureWindowInternal(magX, magY);
  if (pv) pv->SetMagnification(1, 1);

  return image;
}
