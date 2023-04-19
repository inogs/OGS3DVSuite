/*=========================================================================

  Module:    Simulation Main class

  Library to deal with the reading and writing of Simulation files as well as
  the NETCDF4 interface.

  Copyright (c) 2018 Arnau Miro, Simulation
  All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include <cstdio>
#include <cstring>
#include <filesystem>
#include <netcdf>
#include <sstream>
#include <string>
#include <vector>

#ifdef USE_OMP
#include <omp.h>
#endif

#include "OGS/OGS3DVSuite.h"
#include "OGS/netcdfio.h"

#define CLLIND(ii, jj, kk, nx, ny) \
  ((nx - 1) * (ny - 1) * (kk) + (nx - 1) * (jj) + (ii))
#define PNTIND(ii, jj, kk, nx, ny) ((nx) * (ny) * (kk) + (nx) * (jj) + (ii))

using namespace OGS;

/* PROTOTYPES */

char *reads(char *line, int size, FILE *fin);
char *trim(char *str);

/* CLASS FUNCTIONS */

std::unique_ptr<RectilinearMeshVertexCoords>
RectilinearMeshVertexCoords::create_rectilinear_mesh_vertex_coords(
    const std::string &file_name) {
  // Open file for reading
  FILE *myfile;
  myfile = std::fopen(file_name.c_str(), "rb");
  if (myfile == nullptr)
    throw std::ios_base::failure("Cannot open file " + file_name);

  unsigned int nlon = 0;
  unsigned int nlat = 0;
  unsigned int nlev = 0;

  // Read the dimensions
  if (std::fread(&nlon, sizeof(unsigned int), 1, myfile) != 1)
    throw InvalidRectilinearMeshFile("nlon", file_name);
  if (std::fread(&nlat, sizeof(unsigned int), 1, myfile) != 1)
    throw InvalidRectilinearMeshFile("nlat", file_name);
  if (std::fread(&nlev, sizeof(unsigned int), 1, myfile) != 1)
    throw InvalidRectilinearMeshFile("nlev", file_name);

  // Allocate space (and we put everything to 0 to be on the safe side)
  std::vector<double> lon2m(nlon, 0);
  std::vector<double> lat2m(nlat, 0);
  std::vector<double> nav_lev(nlev, 0);

  // Read the vectors
  if (std::fread(lon2m.data(), DBLSZ, nlon, myfile) != nlon)
    throw InvalidRectilinearMeshFile("longitudes", file_name);

  if (std::fread(lat2m.data(), DBLSZ, nlat, myfile) != nlat)
    throw InvalidRectilinearMeshFile("latitudes", file_name);

  if (std::fread(nav_lev.data(), DBLSZ, nlev, myfile) != nlev)
    throw InvalidRectilinearMeshFile("nav_lev", file_name);

  unsigned int ncells = (nlev - 1) * (nlon - 1) * (nlat - 1);

  // Read the masks
  field::Field<uint8_t> *masks[3];

  for (auto mask_type : AllMaskTypes) {
    int m;
    if (std::fread(&m, INTSZ, 1, myfile) != 1)
      throw InvalidRectilinearMeshFile(
          "size of mask " + masktype2string(mask_type), file_name);
    masks[mask_type] = new field::Field<uint8_t>(ncells, m);
    if (std::fread(masks[mask_type]->data(), UI8SZ, masks[mask_type]->get_sz(),
                   myfile) != masks[mask_type]->get_sz())
      throw InvalidRectilinearMeshFile("mask " + masktype2string(mask_type),
                                       file_name);
  }
  std::fclose(myfile);

  return std::make_unique<RectilinearMeshVertexCoords>(
      std::move(lon2m), std::move(lat2m), std::move(nav_lev),
      std::move(*masks[BASINS]), std::move(*masks[COASTS]),
      std::move(*masks[LAND]));
}



void RectilinearMeshVertexCoords::write(const std::string &file_path) {
  // Open file for writing
  FILE *myfile;
  myfile = std::fopen(file_path.c_str(), "wb");
  if (myfile == nullptr)
    throw std::ios_base::failure("Cannot open file " + file_path);

  // Copy the dimensions in a temporary buffer
  unsigned int dims[3];
  dims[0] = this->nlon();
  dims[1] = this->nlat();
  dims[2] = this->nlev();

  // Write the dimensions
  if (std::fwrite(dims, INTSZ, 3, myfile) != 3)
    throw std::ios_base::failure("Error writing mesh dimensions on file " +
                                 file_path);

  // Write the vectors
  if (std::fwrite(this->_lon2m.data(), DBLSZ, this->nlon(), myfile) !=
      (size_t)this->nlon())
    throw std::ios_base::failure("Error writing longitudes on file " +
                                 file_path);
  if (std::fwrite(this->_lat2m.data(), DBLSZ, this->nlat(), myfile) !=
      (size_t)this->nlat())
    throw std::ios_base::failure("Error writing latitudes on file " +
                                 file_path);
  if (std::fwrite(this->_nav_lev.data(), DBLSZ, this->nlev(), myfile) !=
      (size_t)this->nlev())
    throw std::ios_base::failure("Error writing levels on file " + file_path);

  // Write the masks
  for (auto mask_type : AllMaskTypes) {
    auto &current_mask = this->_masks[mask_type];
    unsigned int m = current_mask.get_m();
    if (std::fwrite(&m, INTSZ, 1, myfile) != 1)
      throw std::ios_base::failure("Error writing mesh size for mesh " +
                                   masktype2string(mask_type) + " on file " +
                                   file_path);
    if (std::fwrite(current_mask.data(), UI8SZ, current_mask.get_sz(),
                    myfile) != current_mask.get_sz())
      throw std::ios_base::failure("Error writing mesh data for mesh " +
                                   masktype2string(mask_type) + " on file " +
                                   file_path);
  }

  std::fclose(myfile);
}



linecoords RectilinearMeshVertexCoords::get_horizontal_grid_line(
    unsigned int line_index, unsigned int column_index) const {
  if (line_index >= this->nlat()) {
    std::ostringstream error_msg;
    error_msg << "Line index outside valid range (from 0 to " << this->nlat()
              << "); received " << line_index;
    throw std::invalid_argument(error_msg.str());
  }
  if (column_index >= this->nlev()) {
    std::ostringstream error_msg;
    error_msg << "Column index outside valid range (from 0 to " << this->nlev()
              << "); received " << column_index;
    throw std::invalid_argument(error_msg.str());
  }

  auto output =
      std::make_unique<std::tuple<std::vector<double>, std::vector<double>,
                                  std::vector<double>>>();
  std::get<0>(*output).resize(this->nlon());
  std::get<1>(*output).resize(this->nlon(), this->_lat2m.at(line_index));
  std::get<2>(*output).resize(this->nlon(), this->_nav_lev.at(column_index));

  auto &lon_array = std::get<0>(*output);
  for (unsigned int i = 0; i < this->nlon(); ++i)
    lon_array[i] = this->_lon2m.at(i);

  return output;
}



linecoords RectilinearMeshVertexCoords::get_vertical_grid_line(
    unsigned int line_index, unsigned int column_index) const {
  if (line_index >= this->nlon()) {
    std::ostringstream error_msg;
    error_msg << "Line index outside valid range (from 0 to " << this->nlon()
              << "); received " << line_index;
    throw std::invalid_argument(error_msg.str());
  }
  if (column_index >= this->nlev()) {
    std::ostringstream error_msg;
    error_msg << "Column index outside valid range (from 0 to " << this->nlev()
              << "); received " << column_index;
    throw std::invalid_argument(error_msg.str());
  }

  auto output =
      std::make_unique<std::tuple<std::vector<double>, std::vector<double>,
                                  std::vector<double>>>();
  std::get<0>(*output).resize(this->nlat(), this->_lon2m.at(line_index));
  std::get<1>(*output).resize(this->nlat());
  std::get<2>(*output).resize(this->nlat(), this->_nav_lev.at(column_index));

  auto &lat_array = std::get<1>(*output);
  for (unsigned int i = 0; i < this->nlat(); ++i)
    lat_array[i] = this->_lat2m.at(i);

  return output;
}



linecoords RectilinearMeshVertexCoords::get_column_grid_line(
    unsigned int vline_index, unsigned int hline_index) const {
  if (vline_index >= this->nlon()) {
    std::ostringstream error_msg;
    error_msg << "Line index outside valid range (from 0 to " << this->nlon()
              << "); received " << vline_index;
    throw std::invalid_argument(error_msg.str());
  }

  if (hline_index >= this->nlat()) {
    std::ostringstream error_msg;
    error_msg << "Line index outside valid range (from 0 to " << this->nlat()
              << "); received " << hline_index;
    throw std::invalid_argument(error_msg.str());
  }

  auto output =
      std::make_unique<std::tuple<std::vector<double>, std::vector<double>,
                                  std::vector<double>>>();
  std::get<0>(*output).resize(this->nlev(), this->_lon2m.at(vline_index));
  std::get<1>(*output).resize(this->nlev(), this->_lat2m.at(hline_index));
  std::get<2>(*output).resize(this->nlev());

  auto &lev_array = std::get<2>(*output);
  for (unsigned int i = 0; i < this->nlev(); ++i)
    lev_array[i] = this->_nav_lev.at(i);

  return output;
}



std::unique_ptr<Struct2DMeshVertexCoords>
Struct2DMeshVertexCoords::create_struct2d_mesh_vertex_coords(
    const std::string &file_name) {
  int fid;
  if (nc_open(file_name.c_str(), NC_NOWRITE, &fid) != NC_NOERR)
    throw std::ios_base::failure("Cannot open file " + file_name);

  int nlon_dim_id = 0, nlat_dim_id = 0, nlev_dim_id = 0;
  size_t nlon, nlat, nlev;

  if (nc_inq_dimid(fid, "nlon_vertices", &nlon_dim_id) != NC_NOERR)
    throw InvalidStruct2DMeshFile(
        "Dimension \"nlon_vertices\" not found in file " + file_name);
  if (nc_inq_dimid(fid, "nlat_vertices", &nlat_dim_id) != NC_NOERR)
    throw InvalidStruct2DMeshFile(
        "Dimension \"nlat_vertices\" not found in file " + file_name);
  if (nc_inq_dimid(fid, "nlev_vertices", &nlev_dim_id) != NC_NOERR)
    throw InvalidStruct2DMeshFile(
        "Dimension \"nlev_vertices\" not found in file " + file_name);

  if (nc_inq_dimlen(fid, nlon_dim_id, &nlon) != NC_NOERR)
    throw InvalidStruct2DMeshFile(
        "Error while reading size of dimension \"nlon_vertices\" in file " +
        file_name);
  if (nc_inq_dimlen(fid, nlat_dim_id, &nlat) != NC_NOERR)
    throw InvalidStruct2DMeshFile(
        "Error while reading size of dimension \"nlat_vertices\" in file " +
        file_name);
  if (nc_inq_dimlen(fid, nlev_dim_id, &nlev) != NC_NOERR)
    throw InvalidStruct2DMeshFile(
        "Error while reading size of dimension \"nlev_vertices\" in file " +
        file_name);

  auto lon2m = field::Field<double>(nlat, nlon);
  auto lat2m = field::Field<double>(nlat, nlon);
  auto nav_lev = std::vector<double>(nlev);

  if (NetCDF::NetCDFGetVar(fid, "lon", nlon * nlat, lon2m.data(), false) !=
      NETCDF_OK)
    throw InvalidStruct2DMeshFile(
        "Error while reading variable \"lon\" in file " + file_name);

  if (NetCDF::NetCDFGetVar(fid, "lat", nlon * nlat, lat2m.data(), false) !=
      NETCDF_OK)
    throw InvalidStruct2DMeshFile(
        "Error while reading variable \"lat\" in file " + file_name);

  if (NetCDF::NetCDFGetVar(fid, "nav_lev", nlev, &nav_lev[0], false) !=
      NETCDF_OK)
    throw InvalidStruct2DMeshFile(
        "Error while reading variable \"nav_lev\" in file " + file_name);

  size_t nlon_cells, nlat_cells, nlev_cells;

  if (nc_inq_dimid(fid, "nlon_cells", &nlon_dim_id) != NC_NOERR)
    throw InvalidStruct2DMeshFile(
        "Dimension \"nlon_cells\" not found in file " + file_name);
  if (nc_inq_dimid(fid, "nlat_cells", &nlat_dim_id) != NC_NOERR)
    throw InvalidStruct2DMeshFile(
        "Dimension \"nlat_cells\" not found in file " + file_name);
  if (nc_inq_dimid(fid, "nlev_cells", &nlev_dim_id) != NC_NOERR)
    throw InvalidStruct2DMeshFile(
        "Dimension \"nlev_cells\" not found in file " + file_name);

  if (nc_inq_dimlen(fid, nlon_dim_id, &nlon_cells) != NC_NOERR)
    throw InvalidStruct2DMeshFile(
        "Error while reading size of dimension \"nlon_cells\" in file " +
        file_name);
  if (nc_inq_dimlen(fid, nlat_dim_id, &nlat_cells) != NC_NOERR)
    throw InvalidStruct2DMeshFile(
        "Error while reading size of dimension \"nlat_cells\" in file " +
        file_name);
  if (nc_inq_dimlen(fid, nlev_dim_id, &nlev_cells) != NC_NOERR)
    throw InvalidStruct2DMeshFile(
        "Error while reading size of dimension \"nlev_cells\" in file " +
        file_name);

  if (nlon - nlon_cells != 1)
    throw InvalidStruct2DMeshFile(
        "nlon_vertices and nlon_cells should differ by one in a valid Struct2D "
        "mesh file; error in file: " +
        file_name);
  if (nlat - nlat_cells != 1)
    throw InvalidStruct2DMeshFile(
        "nlat_vertices and nlat_cells should differ by one in a valid Struct2D "
        "mesh file; error in file: " +
        file_name);

  unsigned int ncells = nlon_cells * nlat_cells * nlev_cells;

  field::Field<uint8_t> masks[3];
  masks[BASINS] = field::Field<uint8_t>(ncells, 16);
  masks[COASTS] = field::Field<uint8_t>(ncells, 1);
  masks[LAND] = field::Field<uint8_t>(ncells, 1);

  if (NetCDF::NetCDFGetVar(fid, "coasts_mask", ncells, masks[COASTS].data(),
                           false) != NETCDF_OK)
    throw InvalidStruct2DMeshFile(
        "Error while reading variable \"coasts_mask\" in file " + file_name);

  if (NetCDF::NetCDFGetVar(fid, "land_mask", ncells, masks[LAND].data(),
                           false) != NETCDF_OK)
    throw InvalidStruct2DMeshFile(
        "Error while reading variable \"land_mask\" in file " + file_name);

  auto *temp_array = new uint16_t[ncells];
  if (NetCDF::NetCDFGetVar(fid, "basins_mask", ncells, temp_array, false) !=
      NETCDF_OK)
    throw InvalidStruct2DMeshFile(
        "Error while reading variable \"coasts_mask\" in file " + file_name);

  auto &basins_mask = masks[BASINS];
  basins_mask.set_val(static_cast<uint8_t>(0));
  for (unsigned int i = 0; i < ncells; ++i) {
    auto j = temp_array[i];
    unsigned int k = 0;
    while (j > 0) {
      if (j % 2 == 1) basins_mask[i][k] = 1;
      j = j >> 1;
      k += 1;
    }
  }
  delete[] temp_array;

  nc_close(fid);
  return std::make_unique<Struct2DMeshVertexCoords>(
      std::move(lon2m), std::move(lat2m), std::move(nav_lev),
      std::move(masks[BASINS]), std::move(masks[COASTS]),
      std::move(masks[LAND]));
}



linecoords Struct2DMeshVertexCoords::get_horizontal_grid_line(
    unsigned int line_index, unsigned int column_index) const {
  if (line_index >= this->nlat()) {
    std::ostringstream error_msg;
    error_msg << "Line index outside valid range (from 0 to " << this->nlat()
              << "); received " << line_index;
    throw std::invalid_argument(error_msg.str());
  }
  if (column_index >= this->nlev()) {
    std::ostringstream error_msg;
    error_msg << "Column index outside valid range (from 0 to " << this->nlev()
              << "); received " << column_index;
    throw std::invalid_argument(error_msg.str());
  }

  auto output =
      std::make_unique<std::tuple<std::vector<double>, std::vector<double>,
                                  std::vector<double>>>();
  std::get<0>(*output).resize(this->nlon());
  std::get<1>(*output).resize(this->nlon());
  std::get<2>(*output).resize(this->nlon(), this->_nav_lev.at(column_index));

  const auto lon_data = this->_lon2m[line_index];
  const auto lat_data = this->_lat2m[line_index];

  auto &lon_array = std::get<0>(*output);
  auto &lat_array = std::get<1>(*output);
  for (unsigned int i = 0; i < this->nlon(); ++i) {
    lon_array[i] = lon_data[i];
    lat_array[i] = lat_data[i];
  }

  return output;
}



linecoords Struct2DMeshVertexCoords::get_vertical_grid_line(
    unsigned int line_index, unsigned int column_index) const {
  if (line_index >= this->nlon()) {
    std::ostringstream error_msg;
    error_msg << "Line index outside valid range (from 0 to " << this->nlon()
              << "); received " << line_index;
    throw std::invalid_argument(error_msg.str());
  }
  if (column_index >= this->nlev()) {
    std::ostringstream error_msg;
    error_msg << "Column index outside valid range (from 0 to " << this->nlev()
              << "); received " << column_index;
    throw std::invalid_argument(error_msg.str());
  }

  auto output =
      std::make_unique<std::tuple<std::vector<double>, std::vector<double>,
                                  std::vector<double>>>();
  std::get<0>(*output).resize(this->nlat());
  std::get<1>(*output).resize(this->nlat());
  std::get<2>(*output).resize(this->nlat(), this->_nav_lev.at(column_index));

  auto &lon_array = std::get<0>(*output);
  auto &lat_array = std::get<1>(*output);
  for (unsigned int i = 0; i < this->nlat(); ++i) {
    lon_array[i] = this->_lon2m[i][line_index];
    lat_array[i] = this->_lat2m[i][line_index];
  }
  return output;
}



linecoords Struct2DMeshVertexCoords::get_column_grid_line(
    unsigned int vline_index, unsigned int hline_index) const {
  if (vline_index >= this->nlon()) {
    std::ostringstream error_msg;
    error_msg << "Line index outside valid range (from 0 to " << this->nlon()
              << "); received " << vline_index;
    throw std::invalid_argument(error_msg.str());
  }

  if (hline_index >= this->nlat()) {
    std::ostringstream error_msg;
    error_msg << "Line index outside valid range (from 0 to " << this->nlat()
              << "); received " << hline_index;
    throw std::invalid_argument(error_msg.str());
  }

  auto output =
      std::make_unique<std::tuple<std::vector<double>, std::vector<double>,
                                  std::vector<double>>>();
  std::get<0>(*output).resize(this->nlev(),
                              this->_lon2m[hline_index][vline_index]);
  std::get<1>(*output).resize(this->nlev(),
                              this->_lat2m[hline_index][vline_index]);
  std::get<2>(*output).resize(this->nlev());

  auto &lev_array = std::get<2>(*output);
  for (unsigned int i = 0; i < this->nlev(); ++i)
    lev_array[i] = this->_nav_lev.at(i);

  return output;
}



GeoCoords::GeoCoords(const std::string &file_name) {
  // Open the netCDF file that describe the meshmask
  int fid;
  if (nc_open(file_name.c_str(), NC_NOWRITE, &fid) != NC_NOERR)
    throw std::ios_base::failure("Cannot open file " + file_name);

  int nlon_dim_id = 0, nlat_dim_id = 0, nlev_dim_id = 0;
  size_t nlon, nlat, nlev;

  if (nc_inq_dimid(fid, "x", &nlon_dim_id) != NC_NOERR)
    throw InvalidMeshmaskFile("Dimension \"x\" not found in file " + file_name);
  if (nc_inq_dimid(fid, "y", &nlat_dim_id) != NC_NOERR)
    throw InvalidMeshmaskFile("Dimension \"y\" not found in file " + file_name);
  if (nc_inq_dimid(fid, "z", &nlev_dim_id) != NC_NOERR)
    throw InvalidMeshmaskFile("Dimension \"z\" not found in file " + file_name);

  if (nc_inq_dimlen(fid, nlon_dim_id, &nlon) != NC_NOERR)
    throw InvalidMeshmaskFile(
        "Error while reading size of dimension \"x\" in file " + file_name);
  if (nc_inq_dimlen(fid, nlat_dim_id, &nlat) != NC_NOERR)
    throw InvalidMeshmaskFile(
        "Error while reading size of dimension \"y\" in file " + file_name);
  if (nc_inq_dimlen(fid, nlev_dim_id, &nlev) != NC_NOERR)
    throw InvalidMeshmaskFile(
        "Error while reading size of dimension \"z\" in file " + file_name);

  // Read the fields as double arrays
  unsigned int ncells2D = nlon * nlat;
  unsigned int ncells = ncells2D * nlev;

  std::vector<double> e1t(ncells2D), e1u(ncells2D), e1v(ncells2D),
      e1f(ncells2D);

  if (NetCDF::NetCDFGetVar(fid, "e1t", ncells2D, &(e1t[0]), true) != NETCDF_OK)
    throw InvalidMeshmaskFile("Error while reading variable \"e1t\" in file " +
                              file_name);
  if (NetCDF::NetCDFGetVar(fid, "e1u", ncells2D, &(e1u[0]), true) != NETCDF_OK)
    throw InvalidMeshmaskFile("Error while reading variable \"e1u\" in file " +
                              file_name);
  if (NetCDF::NetCDFGetVar(fid, "e1v", ncells2D, &(e1v[0]), true) != NETCDF_OK)
    throw InvalidMeshmaskFile("Error while reading variable \"e1v\" in file " +
                              file_name);
  if (NetCDF::NetCDFGetVar(fid, "e1f", ncells2D, &(e1f[0]), true) != NETCDF_OK)
    throw InvalidMeshmaskFile("Error while reading variable \"e1f\" in file " +
                              file_name);

  std::vector<double> e2t(ncells2D), e2u(ncells2D), e2v(ncells2D),
      e2f(ncells2D);
  if (NetCDF::NetCDFGetVar(fid, "e2t", ncells2D, &(e2t[0]), true) != NETCDF_OK)
    throw InvalidMeshmaskFile("Error while reading variable \"e2t\" in file " +
                              file_name);
  if (NetCDF::NetCDFGetVar(fid, "e2u", ncells2D, &(e2u[0]), true) != NETCDF_OK)
    throw InvalidMeshmaskFile("Error while reading variable \"e2u\" in file " +
                              file_name);
  if (NetCDF::NetCDFGetVar(fid, "e2v", ncells2D, &(e2v[0]), true) != NETCDF_OK)
    throw InvalidMeshmaskFile("Error while reading variable \"e2v\" in file " +
                              file_name);
  if (NetCDF::NetCDFGetVar(fid, "e2f", ncells2D, &(e2f[0]), true) != NETCDF_OK)
    throw InvalidMeshmaskFile("Error while reading variable \"e2f\" in file " +
                              file_name);

  std::vector<double> e3t(ncells), e3u(ncells), e3v(ncells), e3w(ncells);
  if (NetCDF::NetCDFGetVar(fid, "e3t_0", ncells, &(e3t[0]), true) != NETCDF_OK)
    throw InvalidMeshmaskFile(
        "Error while reading variable \"e3t_0\" in file " + file_name);
  if (NetCDF::NetCDFGetVar(fid, "e3u_0", ncells, &(e3u[0]), true) != NETCDF_OK)
    throw InvalidMeshmaskFile(
        "Error while reading variable \"e3u_0\" in file " + file_name);
  if (NetCDF::NetCDFGetVar(fid, "e3v_0", ncells, &(e3v[0]), true) != NETCDF_OK)
    throw InvalidMeshmaskFile(
        "Error while reading variable \"e3v_0\" in file " + file_name);
  if (NetCDF::NetCDFGetVar(fid, "e3w_0", ncells, &(e3w[0]), true) != NETCDF_OK)
    throw InvalidMeshmaskFile(
        "Error while reading variable \"e3w_0\" in file " + file_name);

  // Allocate the fields
  this->_e1.set_dim(ncells, 4);
  this->_e2.set_dim(ncells, 4);
  this->_e3.set_dim(ncells, 4);

  // Fill the fields
  for (unsigned int kk = 0; kk < nlev; ++kk) {
    for (unsigned int jj = 0; jj < nlat; ++jj) {
      for (unsigned int ii = 0; ii < nlon; ++ii) {
        // e1
        const auto i3d = CLLIND(ii, jj, kk, nlon, nlat);
        const auto i2d = CLLIND(ii, jj, 0, nlon, nlat);
        this->_e1[i3d][0] = e1t[i2d];
        this->_e1[i3d][1] = e1u[i2d];
        this->_e1[i3d][2] = e1v[i2d];
        this->_e1[i3d][3] = e1f[i2d];
        // e2
        this->_e2[i3d][0] = e2t[i2d];
        this->_e2[i3d][1] = e2u[i2d];
        this->_e2[i3d][2] = e2v[i2d];
        this->_e2[i3d][3] = e2f[i2d];
      }
    }
  }
  for (unsigned int i = 0; i < ncells; ++i) {
    this->_e3[i][0] = e3t[i];
    this->_e3[i][1] = e3u[i];
    this->_e3[i][2] = e3v[i];
    this->_e3[i][3] = e3w[i];
  }

  this->_nlon = nlon;
  this->_nlat = nlat;
  this->_nlev = nlev;
}



int Simulation::readMainFile() {
  char line[BUFFSZ];
  char *linep;

  int sec = 0, n_var_r = -1, n_time_r = -1;
  int n_meshes = 0;

  bool bworkdir = false, bmesh = false, btime = false;
  bool bavephys = false, bavefreq = false, bforcings = false, bgenerals = false;
  bool bmeshtype = false;

  // Get the directory of the main file (that we will use later when we create
  // the meshes)
  std::filesystem::path main_file(this->_ogsfile);
  std::filesystem::path main_file_dir =
      std::filesystem::absolute(main_file).parent_path();

  // Open file for reading
  FILE *myfile;
  myfile = std::fopen(this->_ogsfile.c_str(), "r");
  if (myfile == nullptr)
    throw std::ios_base::failure("Cannot open file " + this->_ogsfile);

  // Read file line by line
  while (reads(line, sizeof(line), myfile)) {
    /* IGNORE COMMENTS */
    if (line[0] == '#') continue;

    /* DEFINITION OF SECTIONS */
    if (std::string(line) == "WRKDIR") {
      sec++;
      continue;
    }  // Section containing the working directory
    if (std::string(line) == "MESH") {
      sec++;
      continue;
    }  // Section containing the mesh information
    if (std::string(line) == "AVE_PHYS") {
      sec++;
      continue;
    }  // Section containing the physical variables
    if (std::string(line) == "AVE_FREQ") {
      sec++;
      continue;
    }  // Section containing the biogeochemical variables
    if (std::string(line) == "FORCINGS") {
      sec++;
      continue;
    }  // Section containing the forcings variables
    if (std::string(line) == "GENERALS") {
      sec++;
      continue;
    }  // Section containing general variables (e.g., downloaded from web)
    if (std::string(line) == "TIME") {
      sec++;
      continue;
    }  // Section containing the time-stepping information
    if (std::string(line) == "MESHTYPE") {
      sec++;
      continue;
    }  // Section containing if the mesh is rectilinear or structured

    /* WORKDIR SECTION */
    if (sec == 1 && !bworkdir) {
      auto wrkdir = std::string(trim(line));
      while (wrkdir.find("${THIS_FILE_DIR}") != std::string::npos) {
        const unsigned int str_position = wrkdir.find("${THIS_FILE_DIR}");
        wrkdir = wrkdir.replace(str_position, 16, main_file_dir);
      }
      this->_wrkdir = wrkdir;
      bworkdir = true;
    }

    /* MESH SECTION */
    if (sec == 2 && !bmesh) {
      if (n_var_r >= 0) {
        // If we have read already the total number of meshes, all the other
        // lines are ignored until the next section
        if (n_var_r > n_meshes - 1) {
          bmesh = true;
          n_var_r = -1;
          continue;
        }
        // Read projection name
        std::string mesh_name, mesh_meshfile, mesh_meshmask;
        strtok(line, ":");
        mesh_name = trim(line);

        // Read mesh file
        linep = strtok(NULL, ":");
        mesh_meshfile = trim(linep);

        // Read meshmask file
        linep = strtok(NULL, " ");
        mesh_meshmask = trim(linep);

        // Here we glue together the path of the workdir with the name of the
        // files for obtaining the total path of the files
        std::filesystem::path workdir(this->_wrkdir);
        std::filesystem::path meshfile_filename(mesh_meshfile);
        std::filesystem::path meshmask_filename(mesh_meshmask);

        std::filesystem::path meshfile_filepath = workdir / meshfile_filename;
        std::filesystem::path meshmask_filepath = workdir / meshmask_filename;

        // Add object to array
        this->meshes.emplace_back(mesh_name, meshfile_filepath,
                                  meshmask_filepath);
        n_var_r++;
      } else {
        n_meshes = std::stoi(line);
        this->meshes.clear();
        n_var_r++;
      }
    }

    /* AVE PHYS SECTION */
    if (sec == 3 && !bavephys) {
      if (n_var_r >= 0) {
        if (n_var_r > this->_vars[0].get_nvars() - 1) {
          bavephys = true;
          n_var_r = -1;
          continue;
        }
        // Split the string, read variable name
        strtok(line, ":");
        this->_vars[PHYS].set_name(n_var_r, trim(line));
        // Split the string, read netcfd name
        linep = strtok(NULL, ":");
        this->_vars[PHYS].set_vname(n_var_r, trim(linep));
        // Split the string, read variable path
        linep = strtok(NULL, " ");
        this->_vars[PHYS].set_path(n_var_r, trim(linep));
        n_var_r++;
      } else {
        this->_vars[PHYS].set_nvars(std::stoi(line));
        n_var_r++;
      }
    }

    /* AVE FREQ SECTION */
    if (sec == 4 && !bavefreq) {
      if (n_var_r >= 0) {
        if (n_var_r > this->_vars[1].get_nvars() - 1) {
          bavefreq = true;
          n_var_r = -1;
          continue;
        }
        // Split the string, read variable name
        strtok(line, ":");
        this->_vars[AVE].set_name(n_var_r, trim(line));
        // Split the string, read netcfd name
        linep = strtok(NULL, ":");
        this->_vars[AVE].set_vname(n_var_r, trim(linep));
        // Split the string, read variable path
        linep = strtok(NULL, " ");
        this->_vars[AVE].set_path(n_var_r, trim(linep));
        n_var_r++;
      } else {
        this->_vars[AVE].set_nvars(std::stoi(line));
        n_var_r++;
      }
    }

    /* FORCINGS SECTION */
    if (sec == 5 && !bforcings) {
      if (n_var_r >= 0) {
        if (n_var_r > this->_vars[2].get_nvars() - 1) {
          bforcings = true;
          n_var_r = -1;
          continue;
        }
        // Split the string, read variable name
        strtok(line, ":");
        this->_vars[FORCINGS].set_name(n_var_r, trim(line));
        // Split the string, read netcfd name
        linep = strtok(NULL, ":");
        this->_vars[FORCINGS].set_vname(n_var_r, trim(linep));
        // Split the string, read variable path
        linep = strtok(NULL, " ");
        this->_vars[FORCINGS].set_path(n_var_r, trim(linep));
        n_var_r++;
      } else {
        this->_vars[FORCINGS].set_nvars(std::stoi(line));
        n_var_r++;
      }
    }

    /* GENERALS SECTION */
    if (sec == 6 && !bgenerals) {
      if (n_var_r >= 0) {
        if (n_var_r > this->_vars[3].get_nvars() - 1) {
          bgenerals = true;
          n_var_r = -1;
          continue;
        }
        // Split the string, read variable name
        strtok(line, ":");
        this->_vars[GENERALS].set_name(n_var_r, trim(line));
        // Split the string, read netcfd name
        linep = strtok(NULL, ":");
        this->_vars[GENERALS].set_vname(n_var_r, trim(linep));
        // Split the string, read variable path
        linep = strtok(NULL, " ");
        this->_vars[GENERALS].set_path(n_var_r, trim(linep));
        n_var_r++;
      } else {
        this->_vars[GENERALS].set_nvars(std::stoi(line));
        n_var_r++;
      }
    }

    /* TIMESTEP SECTION */
    if (sec == 7 && !btime) {
      if (n_time_r >= 0) {
        if (n_time_r > this->_ntsteps - 1) {
          btime = true;
          continue;
        }
        // Read datetime
        this->_datetime[n_time_r] = std::string(trim(line));
        n_time_r++;
      } else {
        // Read number of timesteps
        this->_ntsteps = std::stoi(line);
        // Allocate
        this->_datetime.resize(this->_ntsteps);
        n_time_r++;
      }
    }

    /* MESHTYPE SECTION */
    if (sec == 8 && !bmeshtype) {
      auto me = std::string(trim(line));
      if (me == "<RECTILINEAR>")
        this->mesh_type = RECTILINEAR;
      else if (me == "<STRUCT2D>")
        this->mesh_type = STRUCT2D;
      else
        throw InvalidMainFile("Invalid mesh type:" + me);
      bmeshtype = true;
    }
  }

  std::fclose(myfile);
  return 1;
}



void Simulation::print() const {
  std::cout << "Simulation CLASS" << std::endl;
  std::cout << "---------" << std::endl;
  std::cout << "Main file:          " << this->_ogsfile << std::endl;
  std::cout << "Working directory:  " << this->_wrkdir << std::endl
            << std::endl;

  std::cout << "Mesh information" << std::endl;
  std::cout << "* Number of projections: " << this->n_meshes() << std::endl;
  std::cout << std::endl;
  for (Mesh mesh : this->meshes) {
    std::cout << "    Projection:     " << mesh.get_name() << std::endl;
    std::cout << "    Mesh file:      " << mesh.get_meshfile() << std::endl;
    std::cout << "    Meshmask file:  " << mesh.get_meshmask() << std::endl
              << std::endl;
  }

  std::cout << "Variables information" << std::endl;
  for (const auto var_type : AllVarTypes) {
    std::cout << "    Variables " << vartype2string(var_type) << ":    "
              << this->var_n(var_type) << std::endl;
    for (int jj = 0; jj < this->var_n(var_type); jj++)
      std::cout << "      " << jj
                << " Name:       " << this->var_name(var_type, jj) << " ("
                << this->var_vname(var_type, jj) << ")" << std::endl;
  }
}



unsigned int Simulation::get_mesh_index(const std::string &proj_name) const {
  for (unsigned int ii = 0; ii < this->n_meshes(); ++ii)
    if (this->meshes.at(ii).get_name() == proj_name) return ii;
  throw std::invalid_argument("Projection " + proj_name + " not found!");
}



void Simulation::load_mesh(unsigned int mesh_index) {
  // If this mesh is already loaded, there is nothing to do, and we can rely on
  // the already loaded mesh
  if (this->mesh_loaded == (int)mesh_index) return;

  // Ensure that the mesh_index is valid
  if (mesh_index >= this->n_meshes()) {
    std::ostringstream error_msg;
    error_msg << "Can not load mesh number " << mesh_index
              << "; this simulation only contains " << this->n_meshes()
              << " meshes (counting from 0)!" << std::endl;
    throw std::invalid_argument(error_msg.str());
  }

  const auto &mesh_pointer = this->meshes.at(mesh_index);

  std::unique_ptr<MeshVertexCoords> vertex_coords;
  switch (this->mesh_type) {
    case (RECTILINEAR):
      vertex_coords =
          RectilinearMeshVertexCoords::create_rectilinear_mesh_vertex_coords(
              mesh_pointer.get_meshfile());
      break;
    case (STRUCT2D):
      vertex_coords =
          Struct2DMeshVertexCoords::create_struct2d_mesh_vertex_coords(
              mesh_pointer.get_meshfile());
      break;
    default:
      throw std::range_error("Unexpected kind of mesh found");
  }

  auto geo_coords = std::make_unique<GeoCoords>(mesh_pointer.get_meshmask());

  this->mesh_data = std::make_shared<MeshData>(std::move(vertex_coords),
                                               std::move(geo_coords));
  this->mesh_loaded = (int)mesh_index;
}



std::string Simulation::var_path(const VAR_TYPE vrt, const unsigned int jj,
                                 const unsigned int t) const {
  // First of all, we need to replace the "*" in the path with the
  // right datetime
  const auto datetime_str = this->datetime(t);
  std::string file_namemask = this->_vars[vrt].get_path(jj);
  const unsigned int star_pos = file_namemask.find('*');
  std::filesystem::path filename =
      file_namemask.replace(star_pos, 1, datetime_str);

  // Now we can add the previous part of the path to the filename
  std::filesystem::path workdir = this->_wrkdir;
  return workdir / filename;
}



/* INTERFACE C FUNCTIONS FOR PYTHON LIBRARY */

extern "C" {
int OGSWriteMesh(const char *fname, const int nlon, const int nlat,
                 const int nlev, double *lon2m, double *lat2m, double *nav_lev,
                 uint8_t *bmask, uint8_t *cmask, uint8_t *lmask) {
  // Unfortunately, here we are forced to copy the data, to convert C pointers
  // to C++ standard arrays
  auto *lon2m_array = new std::vector<double>(lon2m, lon2m + nlon);
  auto *lat2m_array = new std::vector<double>(lat2m, lat2m + nlat);
  auto *nav_lev_array = new std::vector<double>(nav_lev, nav_lev + nlev);

  unsigned int ncells = (nlon - 1) * (nlat - 1) * (nlev - 1);
  auto *basins_mask = new field::Field<uint8_t>(ncells, 16, bmask);
  auto *coasts_mask = new field::Field<uint8_t>(ncells, 1, cmask);
  auto *land_mask = new field::Field<uint8_t>(ncells, 1, lmask);

  // Create a new instance of the class
  RectilinearMeshVertexCoords mesh(
      std::move(*lon2m_array), std::move(*lat2m_array),
      std::move(*nav_lev_array), std::move(*basins_mask),
      std::move(*coasts_mask), std::move(*land_mask));

  mesh.write(fname);

  return 0;
}
}



/* AUXILIARY FUNCTIONS */

char *reads(char *line, int size, FILE *fin) {
  char *l = std::fgets(line, size, fin);
  if (l != NULL)
    for (int ii = 0; line[ii] != '\0'; ii++)
      if (line[ii] == 10) {
        line[ii] = '\0';
        break;
      }
  return l;
}



char *trim(char *str) {
  char *end;
  while (std::isspace(*str)) str++;  // Trim leading space
  if (*str == 0) return str;         // All spaces?
  // Trim trailing space
  end = str + std::strlen(str) - 1;
  while (end > str && std::isspace(*end)) end--;
  // Write new null terminator
  *(end + 1) = 0;
  return str;
}
