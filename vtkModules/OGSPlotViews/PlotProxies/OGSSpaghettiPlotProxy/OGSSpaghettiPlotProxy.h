/*=========================================================================

  Program:   OGSPlotViews
  Module:    OGSSpaghettiPlotProxy.h

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef OGSSpaghettiPlotProxy_h
#define OGSSpaghettiPlotProxy_h

#include "vtkNew.h"                        // needed for vtkNew.
#include "vtkRemotingViewsPythonModule.h"  // needed for exports
#include "vtkSMPythonViewProxy.h"

class vtkImageData;
class vtkSMProxy;

class VTKREMOTINGVIEWSPYTHON_EXPORT OGSSpaghettiPlotProxy
    : public vtkSMPythonViewProxy {
 public:
  static OGSSpaghettiPlotProxy* New();
  vtkTypeMacro(OGSSpaghettiPlotProxy, vtkSMPythonViewProxy);

 protected:
  OGSSpaghettiPlotProxy();
  ~OGSSpaghettiPlotProxy() override;

  /**
   * Subclasses should override this method to do the actual image capture.
   */
  vtkImageData* CaptureWindowInternal(int magX, int magY) override;

 private:
  OGSSpaghettiPlotProxy(const OGSSpaghettiPlotProxy&) = delete;
  void operator=(const OGSSpaghettiPlotProxy&) = delete;
};

#endif
