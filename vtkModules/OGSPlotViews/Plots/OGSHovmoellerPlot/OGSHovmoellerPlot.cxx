/*=========================================================================

  Program:   OGSPlotViews
  Module:    vtkOGSHovmoellerPlot.cxx

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "OGSHovmoellerPlot.h"

#include <iterator>
#include <sstream>
#include <string>

#include "vtkObjectFactory.h"
#include "vtkPythonRepresentation.h"

vtkStandardNewMacro(OGSHovmoellerPlot);

//----------------------------------------------------------------------------
OGSHovmoellerPlot::OGSHovmoellerPlot() { this->Script = NULL; }

OGSHovmoellerPlot::~OGSHovmoellerPlot() { this->SetScript(NULL); }

//----------------------------------------------------------------------------
void OGSHovmoellerPlot::Update() {
  // Enable all arrays for plotting
  this->Superclass::EnableAllAttributeArrays();
  // First we iterate over the parameters map to create the variable = value
  // pairs
  std::ostringstream buf;
  for (auto it = this->mapParam.begin(); it != mapParam.end(); it++)
    buf << it->first << " = " << it->second << "\n";

  buf << "\n" << this->Script;
  // Update the Script property of the Superclass with the newly generated
  // script
  this->Superclass::SetScript(buf.str().c_str());
  // Finally call the update method of the superclass
  this->Superclass::Update();
}

//----------------------------------------------------------------------------
void OGSHovmoellerPlot::SetParameterInternal(const char* raw_name,
                                             const char* raw_value) {
  const std::string name = raw_name ? raw_name : "";
  const std::string value = raw_value ? raw_value : "";

  if (name.empty()) {
    vtkErrorMacro(<< "cannot set parameter with empty name");
    return;
  }

  // Check if we can insert the parameter in the map
  if (!this->mapParam.insert(std::make_pair(raw_name, raw_value)).second) {
    // We couldn't insert the key, which means it already exists
    // we shall update its value
    mapParam[raw_name] = raw_value;
  }

  this->Modified();
}

void OGSHovmoellerPlot::SetParameter(const char* raw_name, const char* value) {
  std::ostringstream buf;

  buf << "r'" << value << "'";
  this->SetParameterInternal(raw_name, buf.str().c_str());
}

void OGSHovmoellerPlot::SetParameter(const char* raw_name, int value) {
  std::ostringstream buf;
  buf << value;
  this->SetParameterInternal(raw_name, buf.str().c_str());
}

void OGSHovmoellerPlot::SetParameter(const char* raw_name, double value) {
  std::ostringstream buf;
  buf << value;
  this->SetParameterInternal(raw_name, buf.str().c_str());
}

void OGSHovmoellerPlot::SetParameter(const char* raw_name, double value1,
                                     double value2, double value3) {
  std::ostringstream buf;
  buf << "[" << value1 << ", " << value2 << ", " << value3 << "]";
  this->SetParameterInternal(raw_name, buf.str().c_str());
}

void OGSHovmoellerPlot::SetParameter(const char* raw_name, double value1,
                                     double value2) {
  std::ostringstream buf;
  buf << "[" << value1 << ", " << value2 << "]";
  this->SetParameterInternal(raw_name, buf.str().c_str());
}