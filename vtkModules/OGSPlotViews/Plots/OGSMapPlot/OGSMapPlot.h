/*=========================================================================

  Program:   OGSPlotViews
  Module:    OGSMapPlot.h

  Copyright (c) 2019 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkOGSMapPlot_h
#define vtkOGSMapPlot_h

#include <map>

#include "vtkDataSet.h"
#include "vtkPythonView.h"

class OGSMapPlot : public vtkPythonView {
 public:
  static OGSMapPlot* New();
  vtkTypeMacro(OGSMapPlot, vtkPythonView);

  /*
    Get/Set the Python script.
  */
  vtkSetStringMacro(Script);
  vtkGetStringMacro(Script);

  /*
    Get/Set the variables.
  */
  vtkSetStringMacro(Variables);
  vtkGetStringMacro(Variables);

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
  OGSMapPlot();
  ~OGSMapPlot() override;

 private:
  OGSMapPlot(const OGSMapPlot&) = delete;
  void operator=(const OGSMapPlot&) = delete;

  char *Script, *Variables;
  std::map<std::string, std::string> mapParam;
};

#endif