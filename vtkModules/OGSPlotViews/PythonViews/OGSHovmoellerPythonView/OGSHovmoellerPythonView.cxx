/*=========================================================================

  Program:   OGSPlotViews
  Module:    OGSHovmoellerPlot.cxx

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "OGSHovmoellerPythonView.h"

#include "../../PlotProxies/OGSHovmoellerPlotProxy/OGSHovmoellerPlotProxy.h"

//-----------------------------------------------------------------------------
OGSHovmoellerPlot::OGSHovmoellerPlot(const QString& group, const QString& name,
                                     vtkSMProxy* renModule, pqServer* server,
                                     QObject* parent)
    : pqPythonView(OGSHovmoellerPlot::OGSHovmoellerPlotType(), group, name,
                   OGSHovmoellerPlotProxy::SafeDownCast(renModule), server,
                   parent)

{}

//-----------------------------------------------------------------------------
OGSHovmoellerPlot::~OGSHovmoellerPlot() {}
