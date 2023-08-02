/*=========================================================================

  Program:   OGSPlotViews
  Module:    OGSSpaghettiPlot.cxx

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "OGSSpaghettiPythonView.h"

#include "OGSSpaghettiPlotProxy.h"

//-----------------------------------------------------------------------------
OGSSpaghettiPlot::OGSSpaghettiPlot(const QString& group, const QString& name,
                                   vtkSMProxy* renModule, pqServer* server,
                                   QObject* parent)
    : pqPythonView(OGSSpaghettiPlot::OGSSpaghettiPlotType(), group, name,
                   OGSSpaghettiPlotProxy::SafeDownCast(renModule), server,
                   parent)

{}

//-----------------------------------------------------------------------------
OGSSpaghettiPlot::~OGSSpaghettiPlot() {}
