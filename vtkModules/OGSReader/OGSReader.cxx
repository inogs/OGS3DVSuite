/*=========================================================================

  Program:   OGSReader
  Module:    OGSReader.cxx

  Copyright (c) 2018 Arnau Miro, OGS
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "OGSReader.h"

#include <algorithm>
#include <ctime>

#include "vtkCallbackCommand.h"
#include "vtkCellData.h"
#include "vtkCommand.h"
#include "vtkDataArraySelection.h"
#include "vtkFloatArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPoints.h"
#include "vtkRectilinearGrid.h"
#include "vtkStringArray.h"
#include "vtkTypeUInt8Array.h"


vtkStandardNewMacro(OGSReader);

//----------------------------------------------------------------------------
#include "OGS/field.h"
#include "OGS/fieldOperations.h"
#include "OGS/macros.h"
#include "OGS/netcdfio.h"
#include "OGS/utilities.h"
#include "OGS/vtkFields.h"


//----------------------------------------------------------------------------
OGSReader::OGSReader()
    : MaskDataArraySelection(vtkDataArraySelection::New()),
      FileName(nullptr),
      RMeshMask(1),
      DepthScale(1000.),
      AvePhysDataArraySelection(vtkDataArraySelection::New()),
      AveFreqDataArraySelection(vtkDataArraySelection::New()),
      ForcingDataArraySelection(vtkDataArraySelection::New()),
      GeneralDataArraySelection(vtkDataArraySelection::New()),
      abort(0),
      procId(0),
      nProcs(0),
      projId(0),
      Mesh(MeshType::New()),
      Projections(vtkStringArray::New()) {
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
}

//----------------------------------------------------------------------------
OGSReader::~OGSReader() {
  this->FileName = nullptr;

  this->MaskDataArraySelection->Delete();
  this->AvePhysDataArraySelection->Delete();
  this->AveFreqDataArraySelection->Delete();
  this->ForcingDataArraySelection->Delete();
  this->GeneralDataArraySelection->Delete();

  this->Projections->Delete();

  this->Mesh->Delete();
}

void OGSReader::initialize() {
  if (this->FileName == nullptr)
    vtkErrorMacro(
        "No filename has been passed to this reader! Please, set the path of "
        "the main OGS file using the method SetFileName");

  if (this->ogsdata == nullptr) {
    try {
      this->ogsdata = std::make_unique<OGS::Simulation>(this->FileName);
    } catch (...) {
      vtkErrorMacro("Cannot read <" << this->FileName << ">!\nAborting");
      this->abort = 1;
      throw;
    }
    if (this->projName == "UNDEFINED")
      this->projName = this->ogsdata->projection(0);
  }
}

//----------------------------------------------------------------------------
int OGSReader::RequestInformation(
    vtkInformation *vtkNotUsed(request),
    vtkInformationVector **vtkNotUsed(inputVector),
    vtkInformationVector *outputVector) {
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  this->initialize();

  this->UpdateProgress(0.0);

  /* READING THE MESH FILE

      The mesh data (dims, Lon2Meters, Lat2Meters and nav_lev) are inside
      the .ogsmsh file containing the correct resolution of the meshmask for
      the current simulation.

      This is done in the RequestInformation to ensure it is only read once
      and stored in memory.

  */

  // Obtain the projection index
  this->projId =
      static_cast<int>(this->ogsdata->get_mesh_index(this->projName));

  try {
    this->ogsdata->load_mesh(this->projId);
  } catch (...) {
    vtkErrorMacro("Problems reading the mesh!\nAborting.");
    this->abort = 1;
    throw;
  }

  this->UpdateProgress(0.05);

  /* CREATE RECTILINEAR GRID

      The rectilinear grid is created here from the mesh data that has
      been previously read. This is handled in the RequestInformation as it only
     needs to be executed once. Moreover, the rectilinear grid is only created
     when the mesh dimensions are (0,0,0).

  auto horizontal_line =
      this->ogsdata->loaded_mesh()->get_horizontal_grid_line(0, 0);
  auto vertical_line =
      this->ogsdata->loaded_mesh()->get_vertical_grid_line(0, 0);
  auto column = this->ogsdata->loaded_mesh()->get_column_grid_line(0, 0);

  auto lon2meters = &(std::get<0>(*horizontal_line)[0]);
  auto lat2meters = &(std::get<1>(*vertical_line)[0]);
  auto nav_lev = &(std::get<2>(*column)[0]);

  VTK::createRectilinearGrid(this->ogsdata->loaded_mesh()->nlon_vertices(),
                             this->ogsdata->loaded_mesh()->nlat_vertices(),
                             this->ogsdata->loaded_mesh()->nlev_vertices(),
                             lon2meters, lat2meters, nav_lev, this->DepthScale,
                             this->Mesh);

  horizontal_line.reset();
  vertical_line.reset();
  column.reset();
  */

  /* CREATE STRUCTURED GRID */

  unsigned int nx = this->ogsdata->loaded_mesh()->nlon_vertices();
  unsigned int ny = this->ogsdata->loaded_mesh()->nlat_vertices();
  unsigned int nz = this->ogsdata->loaded_mesh()->nlev_vertices();

  // This specifies the valid index inside the grid of the mesh; this
  // information must be returned by the RequestInformation method
  auto whole_extent = new int[6];
  whole_extent[0] = 0;
  whole_extent[1] = static_cast<int>(nx) - 1;
  whole_extent[2] = 0;
  whole_extent[3] = static_cast<int>(ny) - 1;
  whole_extent[4] = 0;
  whole_extent[5] = static_cast<int>(nz) - 1;
  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), whole_extent,
               6);
  delete[] whole_extent;

  vtkNew<vtkPoints> points;
  points->SetNumberOfPoints(nx * ny * nz);

  for (unsigned int k = 0; k < nz; k++) {
    for (unsigned int j = 0; j < ny; j++) {
      auto horizontal_line =
          this->ogsdata->loaded_mesh()->get_horizontal_grid_line(j, k);
      for (unsigned int i = 0; i < nx; i++) {
        double z = std::get<2>(*horizontal_line)[i];
        double y = std::get<1>(*horizontal_line)[i];
        double x = std::get<0>(*horizontal_line)[i];
        points->SetPoint(((k * ny) + j) * nx + i, x, y, -z * this->DepthScale);
      }
    }
  }

  this->Mesh->SetDimensions(static_cast<int>(nx), static_cast<int>(ny),
                            static_cast<int>(nz));
  this->Mesh->SetPoints(points);

  this->UpdateProgress(0.10);

  // Set up the data array selection containing the masks to be loaded.
  this->MaskDataArraySelection->AddArray("Sub-basins");
  this->MaskDataArraySelection->AddArray("Continental shelf");
  this->MaskDataArraySelection->AddArray("Land");

  /* ADD MASK ARRAYS

      Add the masks to the mesh. This process is handled here since it only
     needs to change when the request information executes, and not with the
     time-stepping

      Current masks are:
          > Sub-basins: mask with all the basins of the MED.
          > Continental shelf: area from 0 to 200 m of depth.

  */
  vtkTypeUInt8Array *vtkmask;
  VTKARRAY *vtkarray;

  if (this->GetMaskArrayStatus("Sub-basins")) {
    vtkmask = VTK::createVTKfromField<vtkTypeUInt8Array, uint8_t>(
        "basins mask", this->ogsdata->loaded_mesh()->get_mask(OGS::BASINS));
    this->Mesh->GetCellData()->AddArray(vtkmask);
    vtkmask->Delete();
  } else {
    this->Mesh->GetCellData()->RemoveArray("basins mask");
  }

  // Continental shelf mask ("coast mask")
  if (this->GetMaskArrayStatus("Continental shelf")) {
    vtkmask = VTK::createVTKfromField<vtkTypeUInt8Array, uint8_t>(
        "coast mask", this->ogsdata->loaded_mesh()->get_mask(OGS::COASTS));
    this->Mesh->GetCellData()->AddArray(vtkmask);
    vtkmask->Delete();
  } else {
    this->Mesh->GetCellData()->RemoveArray("coast mask");
  }

  // Continental shelf mask ("coast mask")
  if (this->GetMaskArrayStatus("Land")) {
    vtkmask = VTK::createVTKfromField<vtkTypeUInt8Array, uint8_t>(
        "land mask", this->ogsdata->loaded_mesh()->get_mask(OGS::LAND));
    this->Mesh->GetCellData()->AddArray(vtkmask);
    vtkmask->Delete();
  } else {
    this->Mesh->GetCellData()->RemoveArray("land mask");
  }

  this->UpdateProgress(0.15);

  /* ADD MESHMASK STRETCHING ARRAYS

      Add the stretching arrays e1, e2 and e3 found in the meshmask. These
     arrays are needed to project the velocity from a face centered coordinates
     to cell centered coordinates as well as to compute gradients.

  */

  if (this->RMeshMask) {
    // e1
    vtkarray = VTK::createVTKfromField<VTKARRAY, double>(
        "e1", this->ogsdata->loaded_mesh()->e1());
    this->Mesh->GetCellData()->AddArray(vtkarray);
    vtkarray->Delete();
    // e2
    vtkarray = VTK::createVTKfromField<VTKARRAY, double>(
        "e2", this->ogsdata->loaded_mesh()->e2());
    this->Mesh->GetCellData()->AddArray(vtkarray);
    vtkarray->Delete();
    // e3
    vtkarray = VTK::createVTKfromField<VTKARRAY, double>(
        "e3", this->ogsdata->loaded_mesh()->e3());
    this->Mesh->GetCellData()->AddArray(vtkarray);
    vtkarray->Delete();
  } else {
    this->Mesh->GetCellData()->RemoveArray("e1");
    this->Mesh->GetCellData()->RemoveArray("e2");
    this->Mesh->GetCellData()->RemoveArray("e3");
  }

  this->UpdateProgress(0.20);

  // Detect which physical variables (AVE_PHYS) are present in the master
  // file and list them inside the array selection.
  for (unsigned int ii = 0; ii < this->ogsdata->var_n(OGS::PHYS); ii++)
    this->AvePhysDataArraySelection->AddArray(
        this->ogsdata->var_name(OGS::PHYS, ii).c_str());

  // As we did for the physical variables, now we scan the biogeochemical ones
  for (unsigned int ii = 0; ii < this->ogsdata->var_n(OGS::AVE); ii++)
    this->AveFreqDataArraySelection->AddArray(
        this->ogsdata->var_name(OGS::AVE, ii).c_str());

  // Scan of the forcings variables
  for (unsigned int ii = 0; ii < this->ogsdata->var_n(OGS::FORCINGS); ii++)
    this->ForcingDataArraySelection->AddArray(
        this->ogsdata->var_name(OGS::FORCINGS, ii).c_str());

  // Scan of the generals variables
  for (unsigned int ii = 0; ii < this->ogsdata->var_n(OGS::GENERALS); ii++)
    this->GeneralDataArraySelection->AddArray(
        this->ogsdata->var_name(OGS::GENERALS, ii).c_str());

  /* SET UP THE TIMESTEP

      Time stepping information is contained inside the master file. Here we set
      the time step and the time step range for paraview.

      We store a time_t variable representing the current unix timestamp.

  */
  // Set the time step value
  auto *timeSteps = new double[this->ogsdata->ntsteps()];
  for (int ii = 0; ii < this->ogsdata->ntsteps(); ii++) {
    struct tm tm = {0};
    tm.tm_isdst = 0;
    // Convert to struct tm
    strptime((this->ogsdata->datetime(ii)).c_str(), "%Y%m%d-%H:%M:%S", &tm);
    timeSteps[ii] = difftime(timegm(&tm), 0);
  }
  outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), &timeSteps[0],
               (int)this->ogsdata->ntsteps());

  // Set up the time range
  double timeRange[2] = {timeSteps[0], timeSteps[this->ogsdata->ntsteps() - 1]};
  outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), timeRange, 2);

  delete[] timeSteps;

  this->UpdateProgress(0.25);
  return 1;
}

//----------------------------------------------------------------------------
int OGSReader::RequestData(vtkInformation *vtkNotUsed(request),
                           vtkInformationVector **vtkNotUsed(inputVector),
                           vtkInformationVector *outputVector) {
  if (this->abort) return 0;

  this->initialize();

  // get the info object
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the output
  MeshType *output =
      MeshType::SafeDownCast(outInfo->Get(MeshType::DATA_OBJECT()));

  /* SET THE TIME STEPPING

      The number of time steps are read in the master file. If that number is
      greater than one, then use the ParaView built in functions to move
      through the time steps.
  */
  int ii_tstep = 0;
  if (outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP())) {
    double *timeSteps;
    // Get the requested time step. We only support requests of a single time
    // step in this reader right now
    double requestedTimeValue =
        outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
    output->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(),
                                  requestedTimeValue);
    // Recover the timestep list
    timeSteps = outInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());

    ii_tstep = std::distance(
        timeSteps, std::find(timeSteps, timeSteps + this->ogsdata->ntsteps(),
                             requestedTimeValue));
    if (ii_tstep >= this->ogsdata->ntsteps()) ii_tstep = 0;
  }

  this->UpdateProgress(0.25);

  VTKARRAY *vtkarray;
  OGS::field::Field<FLDARRAY> array, array1;
  int n_vars_loaded = 0;

  /* READING THE PHYSICAL VARIABLES
      Variables inside AVE_PHYS are read here. Velocity is outputed as
      a vector while the others are scalar arrays. User can select which
      variables to load by using a panel.
  */
  for (unsigned int ii = 0; ii < this->ogsdata->var_n(OGS::PHYS); ii++) {
    const std::string var_name = this->ogsdata->var_name(OGS::PHYS, ii);
    const std::string var_vname = this->ogsdata->var_vname(OGS::PHYS, ii);
    const std::string var_path =
        this->ogsdata->var_path(OGS::PHYS, ii, ii_tstep);
    // Test if the variable has been activated
    if (this->GetAvePhysArrayStatus(var_name.c_str())) {
      if (var_name == "Velocity") {
        array.clear();
        array.set_dim(this->ogsdata->loaded_mesh()->ncells(), 3);

        std::vector<std::string> vel_vars;
        OGS::utilities::strsplit(this->ogsdata->var_vname(OGS::PHYS, ii), ",",
                                 vel_vars);

        if (OGS::NetCDF::readNetCDF(var_path.c_str(), vel_vars, array) !=
            NETCDF_OK) {
          vtkErrorMacro("Cannot read NetCDF <Velocity>! Aborting!");
          return 0;
        }

        // We need to project the velocity field from a face centered grid to a
        // cell centered grid
        array1 = OGS::field::UVW2T<FLDARRAY>(
            array, this->ogsdata->loaded_mesh()->e1(),
            this->ogsdata->loaded_mesh()->e2(),
            this->ogsdata->loaded_mesh()->e3(),
            this->ogsdata->loaded_mesh()->nlon_cells(),
            this->ogsdata->loaded_mesh()->nlat_cells(),
            this->ogsdata->loaded_mesh()->nlev_cells());

        vtkarray = VTK::createVTKfromField<VTKARRAY, FLDARRAY>(var_name.c_str(),
                                                               array1);
        this->Mesh->GetCellData()->AddArray(vtkarray);
        vtkarray->Delete();
        array1.clear();

        n_vars_loaded++;
      } else {
        array.clear();
        array.set_dim(this->ogsdata->loaded_mesh()->ncells(), 1);

        if (OGS::NetCDF::readNetCDF(var_path.c_str(), var_vname.c_str(),
                                    array) != NETCDF_OK) {
          vtkErrorMacro("Cannot read NetCDF <" << var_vname << "> inside file: "
                                               << var_path << ">! Aborting!");
          return 0;
        }

        vtkarray = VTK::createVTKfromField<VTKARRAY, FLDARRAY>(
            this->ogsdata->var_name(OGS::PHYS, ii), array);
        this->Mesh->GetCellData()->AddArray(vtkarray);
        vtkarray->Delete();

        n_vars_loaded++;
      }
    } else {
      this->Mesh->GetCellData()->RemoveArray(var_name.c_str());
    }
    this->UpdateProgress(0.25 + ii * (0.175 / this->ogsdata->var_n(OGS::PHYS)));
  }

  /* READING THE BIOGEOCHEMICAL VARIABLES
      Variables inside AVE_FREQ are read here. User can select which
      variables to load by using a panel.
  */
  for (unsigned int ii = 0; ii < this->ogsdata->var_n(OGS::AVE); ii++) {
    const std::string var_name = this->ogsdata->var_name(OGS::AVE, ii);
    const std::string var_vname = this->ogsdata->var_vname(OGS::AVE, ii);

    // Test if the variable has been activated
    if (this->GetAveFreqArrayStatus(var_name.c_str())) {
      array.clear();
      array.set_dim(this->ogsdata->loaded_mesh()->ncells(), 1);

      if (OGS::NetCDF::readNetCDF(
              this->ogsdata->var_path(OGS::AVE, ii, ii_tstep).c_str(),
              var_vname.c_str(), array) != NETCDF_OK) {
        vtkErrorMacro("Cannot read NetCDF <" << var_vname << ">! Aborting!");
        return 0;
      }

      vtkarray = VTK::createVTKfromField<VTKARRAY, FLDARRAY>(var_name, array);
      this->Mesh->GetCellData()->AddArray(vtkarray);
      vtkarray->Delete();

      n_vars_loaded++;
    } else {
      this->Mesh->GetCellData()->RemoveArray(var_name.c_str());
    }
    this->UpdateProgress(0.425 + ii * (0.175 / this->ogsdata->var_n(OGS::AVE)));
  }

  /* READING THE FORCINGS VARIABLES
      Variables inside FORCINCS are read here. User can select which
      variables to load by using a panel.
  */
  for (unsigned int ii = 0; ii < this->ogsdata->var_n(OGS::FORCINGS); ii++) {
    const std::string var_name = this->ogsdata->var_name(OGS::FORCINGS, ii);
    const std::string var_vname = this->ogsdata->var_vname(OGS::FORCINGS, ii);
    // Test if the variable has been activated
    if (this->GetForcingArrayStatus(var_name.c_str())) {
      array.clear();
      array.set_dim(this->ogsdata->loaded_mesh()->ncells(), 1);

      if (OGS::NetCDF::readNetCDF(
              this->ogsdata->var_path(OGS::FORCINGS, ii, ii_tstep).c_str(),
              var_vname.c_str(), array) != NETCDF_OK) {
        vtkErrorMacro("Cannot read NetCDF <" << var_vname << ">! Aborting!");
        return 0;
      }

      vtkarray = VTK::createVTKfromField<VTKARRAY, FLDARRAY>(var_name, array);
      this->Mesh->GetCellData()->AddArray(vtkarray);
      vtkarray->Delete();

      n_vars_loaded++;
    } else {
      this->Mesh->GetCellData()->RemoveArray(var_name.c_str());
    }
    this->UpdateProgress(0.6 +
                         ii * (0.175 / this->ogsdata->var_n(OGS::FORCINGS)));
  }

  /* READING THE GENERAL VARIABLES
      Variables inside GENERAL are read here. User can select which
      variables to load by using a panel.
  */
  for (unsigned int ii = 0; ii < this->ogsdata->var_n(OGS::GENERALS); ii++) {
    const std::string var_name = this->ogsdata->var_name(OGS::GENERALS, ii);
    const std::string var_vname = this->ogsdata->var_vname(OGS::GENERALS, ii);
    // Test if the variable has been activated
    if (this->GetGeneralArrayStatus(var_name.c_str())) {
      array.clear();
      array.set_dim(this->ogsdata->loaded_mesh()->ncells(), 1);

      if (OGS::NetCDF::readNetCDF(
              this->ogsdata->var_path(OGS::GENERALS, ii, ii_tstep).c_str(),
              var_vname.c_str(), array) != NETCDF_OK) {
        vtkErrorMacro("Cannot read NetCDF <" << var_vname << ">! Aborting!");
        return 0;
      }

      vtkarray = VTK::createVTKfromField<VTKARRAY, FLDARRAY>(var_name, array);
      this->Mesh->GetCellData()->AddArray(vtkarray);
      vtkarray->Delete();

      n_vars_loaded++;
    } else {
      this->Mesh->GetCellData()->RemoveArray(var_name.c_str());
    }
    this->UpdateProgress(0.775 +
                         ii * (0.175 / this->ogsdata->var_n(OGS::GENERALS)));
  }

  /* SET THE METADATA ARRAY
      Set the metadata, an array that contains multiple information for
      the performance of the OGS filters. Data:
          0 -> Date
          1 -> Datevec
          2 -> Conversion factors
          3 -> Loaded variables
          4 -> File name
          5 -> Meshmask
          6 -> Meshfile
          7 -> Projection ID
  */
  std::string aux_str;
  vtkStringArray *vtkmetadata = VTK::createVTKstrf("Metadata", 8, NULL);

  // Set the current file date
  vtkmetadata->SetValue(0, this->ogsdata->datetime(ii_tstep));

  // Set the datevec
  aux_str = std::to_string(this->ogsdata->ntsteps()) + std::string(";");
  for (int ii = 0; ii < this->ogsdata->ntsteps(); ii++)
    aux_str += std::string(this->ogsdata->datetime(ii)) + std::string(";");
  vtkmetadata->SetValue(1, aux_str.c_str());

  // Set conversion factors
  aux_str = std::to_string(this->DepthScale);
  vtkmetadata->SetValue(2, aux_str.c_str());

  // Set the number of variables loaded
  aux_str = std::to_string(n_vars_loaded) + std::string(";");
  // AVE_PHYS
  for (int ii = 0; ii < this->ogsdata->var_n(OGS::PHYS); ii++)
    if (this->GetAvePhysArrayStatus(
            this->ogsdata->var_name(OGS::PHYS, ii).c_str()))
      aux_str += std::string(this->ogsdata->var_name(OGS::PHYS, ii)) +
                 std::string(";");
  // AVE_FREQ
  for (int ii = 0; ii < this->ogsdata->var_n(OGS::AVE); ii++)
    if (this->GetAveFreqArrayStatus(
            this->ogsdata->var_name(OGS::AVE, ii).c_str()))
      aux_str +=
          std::string(this->ogsdata->var_name(OGS::AVE, ii)) + std::string(";");
  // FORCINGS
  for (int ii = 0; ii < this->ogsdata->var_n(OGS::FORCINGS); ii++)
    if (this->GetForcingArrayStatus(
            this->ogsdata->var_name(OGS::FORCINGS, ii).c_str()))
      aux_str += std::string(this->ogsdata->var_name(OGS::FORCINGS, ii)) +
                 std::string(";");
  // GENERALS
  for (int ii = 0; ii < this->ogsdata->var_n(OGS::GENERALS); ii++)
    if (this->GetGeneralArrayStatus(
            this->ogsdata->var_name(OGS::GENERALS, ii).c_str()))
      aux_str += std::string(this->ogsdata->var_name(OGS::GENERALS, ii)) +
                 std::string(";");
  vtkmetadata->SetValue(3, aux_str.c_str());

  // Set the file name
  vtkmetadata->SetValue(4, this->FileName);

  // Set the mesh mask
  vtkmetadata->SetValue(5, aux_str.c_str());

  // Set the mesh file
  vtkmetadata->SetValue(6, aux_str.c_str());

  // Set the projection id
  vtkmetadata->SetValue(7, this->projName.c_str());

  // Add array to mesh
  this->Mesh->GetFieldData()->AddArray(vtkmetadata);
  vtkmetadata->Delete();

  /*
      FINAL TOUCHES
  */
  output->ShallowCopy(this->Mesh);
  this->UpdateProgress(1.0);

  // Function exit status successful
  return 1;
}

//----------------------------------------------------------------------------
void OGSReader::DisableAllAvePhysArrays() {
  this->AvePhysDataArraySelection->DisableAllArrays();
}

void OGSReader::EnableAllAvePhysArrays() {
  this->AvePhysDataArraySelection->EnableAllArrays();
}

int OGSReader::GetNumberOfAvePhysArrays() {
  return this->AvePhysDataArraySelection->GetNumberOfArrays();
}

const char *OGSReader::GetAvePhysArrayName(int index) {
  if (index >= (int)this->GetNumberOfAvePhysArrays() || index < 0)
    return NULL;
  else
    return this->AvePhysDataArraySelection->GetArrayName(index);
}

int OGSReader::GetAvePhysArrayIndex(const char *name) {
  return this->AvePhysDataArraySelection->GetArrayIndex(name);
}

int OGSReader::GetAvePhysArrayStatus(const char *name) {
  return this->AvePhysDataArraySelection->ArrayIsEnabled(name);
}

void OGSReader::SetAvePhysArrayStatus(const char *name, int status) {
  if (status)
    this->AvePhysDataArraySelection->EnableArray(name);
  else
    this->AvePhysDataArraySelection->DisableArray(name);

  this->Modified();
}

//----------------------------------------------------------------------------
void OGSReader::DisableAllAveFreqArrays() {
  this->AveFreqDataArraySelection->DisableAllArrays();
}

void OGSReader::EnableAllAveFreqArrays() {
  this->AveFreqDataArraySelection->EnableAllArrays();
}

int OGSReader::GetNumberOfAveFreqArrays() {
  return this->AveFreqDataArraySelection->GetNumberOfArrays();
}

const char *OGSReader::GetAveFreqArrayName(int index) {
  if (index >= (int)this->GetNumberOfAveFreqArrays() || index < 0)
    return nullptr;
  else
    return this->AveFreqDataArraySelection->GetArrayName(index);
}

int OGSReader::GetAveFreqArrayIndex(const char *name) {
  return this->AveFreqDataArraySelection->GetArrayIndex(name);
}

int OGSReader::GetAveFreqArrayStatus(const char *name) {
  return this->AveFreqDataArraySelection->ArrayIsEnabled(name);
}

void OGSReader::SetAveFreqArrayStatus(const char *name, int status) {
  if (status)
    this->AveFreqDataArraySelection->EnableArray(name);
  else
    this->AveFreqDataArraySelection->DisableArray(name);

  this->Modified();
}

//----------------------------------------------------------------------------
void OGSReader::DisableAllMaskArrays() {
  this->MaskDataArraySelection->DisableAllArrays();
}

void OGSReader::EnableAllMaskArrays() {
  this->MaskDataArraySelection->EnableAllArrays();
}

int OGSReader::GetNumberOfMaskArrays() {
  return this->MaskDataArraySelection->GetNumberOfArrays();
}

const char *OGSReader::GetMaskArrayName(int index) {
  if (index >= (int)this->GetNumberOfMaskArrays() || index < 0)
    return NULL;
  else
    return this->MaskDataArraySelection->GetArrayName(index);
}

int OGSReader::GetMaskArrayIndex(const char *name) {
  return this->MaskDataArraySelection->GetArrayIndex(name);
}

int OGSReader::GetMaskArrayStatus(const char *name) {
  return this->MaskDataArraySelection->ArrayIsEnabled(name);
}

void OGSReader::SetMaskArrayStatus(const char *name, int status) {
  if (status)
    this->MaskDataArraySelection->EnableArray(name);
  else
    this->MaskDataArraySelection->DisableArray(name);

  this->Modified();
}

//----------------------------------------------------------------------------
void OGSReader::DisableAllForcingArrays() {
  this->ForcingDataArraySelection->DisableAllArrays();
}

void OGSReader::EnableAllForcingArrays() {
  this->ForcingDataArraySelection->EnableAllArrays();
}

int OGSReader::GetNumberOfForcingArrays() {
  return this->ForcingDataArraySelection->GetNumberOfArrays();
}

const char *OGSReader::GetForcingArrayName(int index) {
  if (index >= (int)this->GetNumberOfForcingArrays() || index < 0)
    return nullptr;
  else
    return this->ForcingDataArraySelection->GetArrayName(index);
}

int OGSReader::GetForcingArrayIndex(const char *name) {
  return this->ForcingDataArraySelection->GetArrayIndex(name);
}

int OGSReader::GetForcingArrayStatus(const char *name) {
  return this->ForcingDataArraySelection->ArrayIsEnabled(name);
}

void OGSReader::SetForcingArrayStatus(const char *name, int status) {
  if (status)
    this->ForcingDataArraySelection->EnableArray(name);
  else
    this->ForcingDataArraySelection->DisableArray(name);

  this->Modified();
}

//----------------------------------------------------------------------------
void OGSReader::DisableAllGeneralArrays() {
  this->GeneralDataArraySelection->DisableAllArrays();
}

void OGSReader::EnableAllGeneralArrays() {
  this->GeneralDataArraySelection->EnableAllArrays();
}

int OGSReader::GetNumberOfGeneralArrays() {
  return this->GeneralDataArraySelection->GetNumberOfArrays();
}

const char *OGSReader::GetGeneralArrayName(int index) {
  if (index >= (int)this->GetNumberOfGeneralArrays() || index < 0)
    return nullptr;
  else
    return this->GeneralDataArraySelection->GetArrayName(index);
}

int OGSReader::GetGeneralArrayIndex(const char *name) {
  return this->GeneralDataArraySelection->GetArrayIndex(name);
}

int OGSReader::GetGeneralArrayStatus(const char *name) {
  return this->GeneralDataArraySelection->ArrayIsEnabled(name);
}

void OGSReader::SetGeneralArrayStatus(const char *name, int status) {
  if (status)
    this->GeneralDataArraySelection->EnableArray(name);
  else
    this->GeneralDataArraySelection->DisableArray(name);

  this->Modified();
}

//----------------------------------------------------------------------------
vtkStringArray *OGSReader::GetProjections() {
  this->Projections->Delete();
  this->Projections = VTK::createVTKstrf(
      "Projections", static_cast<int>(this->ogsdata->n_projections()), nullptr);

  for (unsigned int ii = 0; ii < this->ogsdata->n_projections(); ++ii)
    this->Projections->SetValue(ii, this->ogsdata->projection(ii));

  return this->Projections;
}

void OGSReader::SetProjection(const char *proj) {
  this->projName = std::string(proj);
  this->Modified();
}
