/*=========================================================================

        Program:   OGSUtils
        Module:    vtkOGSPointSource.cxx

        Copyright (c) 2018 Arnau Miro, OGS
        All rights reserved.

                 This software is distributed WITHOUT ANY WARRANTY; without even
                 the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "OGSPointSource.h"

#include <cfloat>
#include <cmath>

#include "OGS/projection.h"
#include "ProjectionAssociation.h"
#include "vtkCellArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkRandomSequence.h"

vtkStandardNewMacro(OGSPointSource);

//---------------------------------------------------------------------------
// Specify a random sequence, or use the non-threadsafe one in vtkMath by
// default.
vtkCxxSetObjectMacro(OGSPointSource, RandomSequence, vtkRandomSequence);

//----------------------------------------------------------------------------
OGSPointSource::OGSPointSource(vtkIdType numPts) {
  this->NumberOfPoints = (numPts > 0 ? numPts : 1);

  this->Center[0] = 0.0;
  this->Center[1] = 0.0;
  this->Center[2] = 0.0;

  this->Radius = 0.;
  this->DepthScale = 1000.;
  this->Projection = 0;

  this->Distribution = VTK_POINT_UNIFORM;
  this->OutputPointsPrecision = SINGLE_PRECISION;
  this->RandomSequence = nullptr;

  this->SetNumberOfInputPorts(0);
}

//----------------------------------------------------------------------------
OGSPointSource::~OGSPointSource() { this->SetRandomSequence(nullptr); }

//----------------------------------------------------------------------------
int OGSPointSource::RequestData(vtkInformation *vtkNotUsed(request),
                                vtkInformationVector **vtkNotUsed(inputVector),
                                vtkInformationVector *outputVector) {
  // get the info object
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the output
  vtkPolyData *output =
      vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkIdType i;
  double theta, rho, cosphi, sinphi, radius;
  double x[3];
  vtkPoints *newPoints;
  vtkCellArray *newVerts;

  newPoints = vtkPoints::New();

  // Set the desired precision for the points in the output.
  if (this->OutputPointsPrecision == vtkAlgorithm::DOUBLE_PRECISION)
    newPoints->SetDataType(VTK_DOUBLE);
  else
    newPoints->SetDataType(VTK_FLOAT);

  newPoints->Allocate(this->NumberOfPoints);
  newVerts = vtkCellArray::New();
  newVerts->Allocate(newVerts->EstimateSize(1, this->NumberOfPoints));

  newVerts->InsertNextCell(this->NumberOfPoints);

  if (this->Distribution ==
      VTK_POINT_SHELL) {  // only produce points on the surface of the sphere
    for (i = 0; i < this->NumberOfPoints; i++) {
      cosphi = 1 - 2 * this->Random();
      sinphi = sqrt(1 - cosphi * cosphi);
      radius = this->Radius * sinphi;
      theta = 2.0 * vtkMath::Pi() * this->Random();
      x[0] = this->Center[0] + radius * cos(theta);
      x[1] = this->Center[1] + radius * sin(theta);
      x[2] = this->Center[2] + this->Radius * cosphi;
      newVerts->InsertCellPoint(newPoints->InsertNextPoint(x));
    }
  } else {  // uniform distribution throughout the sphere volume
    for (i = 0; i < this->NumberOfPoints; i++) {
      cosphi = 1 - 2 * this->Random();
      sinphi = sqrt(1 - cosphi * cosphi);
      rho = this->Radius * pow(this->Random(), 0.33333333);
      radius = rho * sinphi;
      theta = 2.0 * vtkMath::Pi() * this->Random();
      x[0] = this->Center[0] + radius * cos(theta);
      x[1] = this->Center[1] + radius * sin(theta);
      x[2] = this->Center[2] + rho * cosphi;
      newVerts->InsertCellPoint(newPoints->InsertNextPoint(x));
    }
  }

  // Update ourselves and release memory
  output->SetPoints(newPoints);
  newPoints->Delete();

  output->SetVerts(newVerts);
  newVerts->Delete();

  return 1;
}

//----------------------------------------------------------------------------
double OGSPointSource::Random() {
  if (!this->RandomSequence) return vtkMath::Random();

  this->RandomSequence->Next();
  return this->RandomSequence->GetValue();
}

// ----------------------------------------------------------------------
void OGSPointSource::SetLonLat(double lon, double lat) {
  // Conversion of lon, lat to a projected point
  OGS::PROJ::Projection p("degrees", associate_projection(this->Projection));
  const auto point_coords = p.transform_point(lon, lat);
  // Set center
  this->Center[0] = point_coords[0];
  this->Center[1] = point_coords[1];
  this->Modified();
}

void OGSPointSource::GetLonLat(double &lon, double &lat) {
  // Conversion to lon, lat
  OGS::PROJ::Projection p(associate_projection(this->Projection), "degrees");
  const auto point_coords = p.transform_point(this->Center[0], this->Center[2]);
  lon = point_coords[0];
  lat = point_coords[1];
}

// ----------------------------------------------------------------------
void OGSPointSource::SetDepth(double depth) {
  this->Center[2] = -depth * this->DepthScale;
  this->Modified();
}
void OGSPointSource::GetDepth(double &depth) {
  depth = -this->Center[2] / this->DepthScale;
}