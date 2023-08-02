/*=========================================================================

  Program:   OGSPlotViews
  Module:    OGSHovmoellerPlotProxy.h

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkSMOGSHovmoellerPlotProxy_h
#define vtkSMOGSHovmoellerPlotProxy_h

// #include "vtkPVServerManagerRenderingModule.h" //needed for exports

#include "vtkNew.h"  // needed for vtkNew.
#include "vtkSMPythonViewProxy.h"

class vtkImageData;
class vtkSMProxy;

class OGSHovmoellerPlotProxy : public vtkSMPythonViewProxy {
 public:
  static OGSHovmoellerPlotProxy* New();
  vtkTypeMacro(OGSHovmoellerPlotProxy, vtkSMPythonViewProxy);

 protected:
  OGSHovmoellerPlotProxy();
  ~OGSHovmoellerPlotProxy() override;

  /**
   * Subclasses should override this method to do the actual image capture.
   */
  vtkImageData* CaptureWindowInternal(int magX, int magY) override;

 private:
  OGSHovmoellerPlotProxy(const OGSHovmoellerPlotProxy&) = delete;
  void operator=(const OGSHovmoellerPlotProxy&) = delete;
};

#endif
