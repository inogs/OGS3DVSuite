/*=========================================================================

  Program:   OGSPlotViews
  Module:    vtkOGSHovmoellerPlot.h

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSHovmoellerPlot_h
#define vtkOGSHovmoellerPlot_h

#include <map>

#include "vtkDataSet.h"
#include "vtkPythonView.h"

class OGSHovmoellerPlot : public vtkPythonView {
 public:
  static OGSHovmoellerPlot* New();
  vtkTypeMacro(OGSHovmoellerPlot, vtkPythonView);

  /*
    Get/Set the Python script.
  */
  vtkSetStringMacro(Script);
  vtkGetStringMacro(Script);

  /*
        Set a name-value parameter that will be available to the script
        when it is run
  */
  void SetParameterInternal(const char* name, const char* value);
  void SetParameter(const char* name, const char* value);
  void SetParameter(const char* name, int value);
  void SetParameter(const char* name, double value);
  void SetParameter(const char* name, double value1, double value2);
  void SetParameter(const char* name, double value1, double value2,
                    double value3);

  /*
    Overrides the base class method to request an addition pass that moves data
    from the server to the client.
  */
  void Update() override;

 protected:
  OGSHovmoellerPlot();
  ~OGSHovmoellerPlot() override;

 private:
  OGSHovmoellerPlot(const OGSHovmoellerPlot&) = delete;
  void operator=(const OGSHovmoellerPlot&) = delete;

  char* Script;
  std::map<std::string, std::string> mapParam;
};

#endif