/*=========================================================================

  Program:   OGSAnnotateDateTime
  Module:    vtkOGSAnnotateDateTime.h

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef OGSAnnotateDateTime_h
#define OGSAnnotateDateTime_h

#include "vtkTimeToTextConvertor.h"


class VTK_EXPORT OGSAnnotateDateTime : public vtkTimeToTextConvertor {
 public:
  static OGSAnnotateDateTime* New();
  vtkTypeMacro(OGSAnnotateDateTime, vtkTimeToTextConvertor);

  OGSAnnotateDateTime(const OGSAnnotateDateTime&) = delete;
  void operator=(const OGSAnnotateDateTime&) = delete;

  vtkSetStringMacro(TimeFormat);
  vtkGetStringMacro(TimeFormat);

  vtkGetMacro(useMetadata, int);
  vtkSetMacro(useMetadata, int);
  vtkBooleanMacro(useMetadata, int);

 protected:
  OGSAnnotateDateTime();
  ~OGSAnnotateDateTime() override;

  int RequestData(vtkInformation*, vtkInformationVector**,
                  vtkInformationVector*) override;

 private:
  char* TimeFormat;
  int useMetadata, procId;
};

#endif
