/*=========================================================================

  Program:   OGSPlotViews
  Module:    OGSMapPlotProxy.h

  Copyright (c) 2019 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkSMOGSMapPlotProxy_h
#define vtkSMOGSMapPlotProxy_h

// #include "vtkPVServerManagerRenderingModule.h" //needed for exports

#include "vtkNew.h"  // needed for vtkNew.
#include "vtkSMPythonViewProxy.h"

class vtkImageData;
class vtkSMProxy;

class OGSMapPlotProxy : public vtkSMPythonViewProxy {
 public:
  static OGSMapPlotProxy* New();
  vtkTypeMacro(OGSMapPlotProxy, vtkSMPythonViewProxy);

 protected:
  OGSMapPlotProxy();
  ~OGSMapPlotProxy() override;

  /**
   * Subclasses should override this method to do the actual image capture.
   */
  vtkImageData* CaptureWindowInternal(int magX, int magY) override;

 private:
  OGSMapPlotProxy(const OGSMapPlotProxy&) = delete;
  void operator=(const OGSMapPlotProxy&) = delete;
};

#endif
