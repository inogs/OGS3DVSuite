/*=========================================================================

Module:    Utilities

Define macros for all the ParaView filters.

Copyright (c) 2020 Stefano Piani and Arnau Miro, OGS
All rights reserved.

   This software is distributed WITHOUT ANY WARRANTY; without even
   the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
   PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef OGS3DVSUITE_UTILITIES_H
#define OGS3DVSUITE_UTILITIES_H
#include <vtkDataSet.h>
#include <vtkInformation.h>

#include <string>
#include <vector>

#include "TimeList.h"

namespace OGS::utilities {
void strsplit(std::string str, std::string splitBy,
              std::vector<std::string> &tokens);

void BuildTimeList(OGS::Time::TimeList &TL, vtkInformation *Info);
void RecoverMasterFileName(std::string &fname, vtkDataSet *input);

}  // namespace OGS::utilities
#endif
