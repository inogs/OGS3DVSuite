/*=========================================================================

Module:    ParaviewPlaceholders

This module imports all the headers files that we use that are not available
when we rely on a VTK installation that is not part of a ParaView compilation



Copyright (c) 2023 Stefano Piani, OGS
All rights reserved.

   This software is distributed WITHOUT ANY WARRANTY; without even
   the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
   PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef OGS3DVSUITE_PARAVIEWIMPORTS_H
#define OGS3DVSUITE_PARAVIEWIMPORTS_H
#ifdef WITH_PARAVIEW
#include "vtkPVVersion.h"  // For PARAVIEW_USE_MPI
#endif                     // WITH_PARAVIEW
#endif                     // OGS3DVSUITE_PARAVIEWIMPORTS_H
