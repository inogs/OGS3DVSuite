/*=========================================================================

  Module:    Simulation Main class

  Library to deal with the reading and writing of simulation files as well as
  the NETCDF4 interface.

  Copyright (c) 2018 Arnau Miro, Simulation
  All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef OGS_H
#define OGS_H

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "field.h"

/* SIZES FOR MALLOCS */

#define INTSZ sizeof(int)
#define UI8SZ sizeof(uint8_t)
#define DBLSZ sizeof(double)

/* BUFFER SIZE FOR FILE READING */

#define VARSZ 80
#define BUFFSZ 256
#define LBUFFSZ 512

namespace OGS {

// The variables of the model are divided in four categories, that are
// represented by the following enum
enum VAR_TYPE { PHYS = 0, AVE = 1, FORCINGS = 2, GENERALS = 3 };

// This allows to loop over all the var types
const VAR_TYPE AllVarTypes[4] = {PHYS, AVE, FORCINGS, GENERALS};

// Convenient function for translating a var type to a string
inline std::string vartype2string(VAR_TYPE var_type) {
  switch (var_type) {
    case PHYS:
      return "PHYS";
    case AVE:
      return "AVE";
    case FORCINGS:
      return "FORCINGS";
    case GENERALS:
      return "GENERALS";
    default:
      return "UNKNOWN VAR TYPE";
  }
}

// These values, instead, identify the kind of mask that a mesh stores
enum MASK_TYPE { BASINS = 0, COASTS = 1, LAND = 2 };

// Convenient function for translating a mask to a string
inline std::string masktype2string(MASK_TYPE mask_type) {
  switch (mask_type) {
    case BASINS:
      return "basins";
    case COASTS:
      return "coasts";
    case LAND:
      return "land";
    default:
      return "UNKNOWN MASK";
  }
}

// This allows to loop over all the mask types
const MASK_TYPE AllMaskTypes[3] = {BASINS, COASTS, LAND};

// This last enum is for the different kinds of meshes
enum MESH_TYPE { RECTILINEAR = 0, STRUCT2D = 1 };

inline std::string meshtype2string(MESH_TYPE mesh_type) {
  switch (mesh_type) {
    case RECTILINEAR:
      return "RECTILINEAR";
    case STRUCT2D:
      return "STRUCT2D";
    default:
      return "UNKNOWN MESH TYPE";
  }
}

// This is the error that we will raise when we found invalid values inside the
// OGS file, i.e., the ASCII file that describes the structure of a Simulation
class InvalidMainFile : std::ios_base::failure {
 public:
  explicit InvalidMainFile(const std::string &failure_) : failure(failure_) {}
};

// This is the error that we will raise when we can not read a RectilinearMesh
// file (which is a custom binary file)
class InvalidRectilinearMeshFile : std::ios_base::failure {
 public:
  InvalidRectilinearMeshFile(std::string field_, std::string file_path_)
      : failure("Error reading field <" + field_ + "> of file " + file_path_),
        field(std::move(field_)),
        file_path(std::move(file_path_)) {}

 private:
  const std::string field, file_path;
};

// This is the error that we will raise when we can not read a Struct2DMesh
// netCDF4 file
class InvalidStruct2DMeshFile : std::ios_base::failure {
 public:
  explicit InvalidStruct2DMeshFile(const std::string &failure_)
      : failure(failure_) {}
};

// This is the error that we will raise when we can not read a meshmask netCDF4
// file (i.e., when some expected fields are missing)
class InvalidMeshmaskFile : std::ios_base::failure {
 public:
  explicit InvalidMeshmaskFile(const std::string &failure_)
      : failure(failure_) {}
};

/* Simulation VARIABLE

        Stores an array with the name of the variable, its name inside
        the NETCDF and the relative path to the netcdf file.

*/
class OGS_VAR {
 public:
  OGS_VAR() = default;
  ~OGS_VAR() = default;

  [[nodiscard]] inline int get_nvars() const { return this->_n; }
  inline void set_nvars(int n) {
    this->_n = n;
    this->allocate();
  }

  inline void set_name(unsigned int i, char *str) {
    this->_name[i] = std::string(str);
  }
  inline void set_vname(unsigned int i, char *str) {
    this->_vname[i] = std::string(str);
  }
  inline void set_path(unsigned int i, char *str) {
    this->_path[i] = std::string(str);
  }

  [[nodiscard]] inline const char *get_name(unsigned int i) const {
    return this->_name.at(i).c_str();
  }
  [[nodiscard]] inline const char *get_vname(unsigned int i) const {
    return this->_vname.at(i).c_str();
  }
  [[nodiscard]] inline const char *get_path(unsigned int i) const {
    return this->_path.at(i).c_str();
  }

  int find_name(const char *str) const {
    for (unsigned int ii = 0; ii < this->_name.size(); ++ii)
      if (this->_name[ii] == std::string(str)) return (int)ii;
    return -1;
  }

 private:
  int _n = 0;
  std::vector<std::string> _name, _vname, _path;

  inline void allocate() {
    this->_name.resize(this->_n);
    this->_vname.resize(this->_n);
    this->_path.resize(this->_n);
  }
};

typedef std::unique_ptr<
    std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>>
    linecoords;

/* RectilinearMeshVertexCoords

      A RectilinearMeshVertexCoords contains all the information about the
   geometry of a mesh whose domain is a rectangle made of rectangular cells for
   each level. On a fixed level, the mesh geometry is obtained by subdividing
   the domain with horizontal and vertical lines (the position of the vertical
      lines is saved inside the lon2m array, the position of the horizontal ones
      is in the lat2m array)

 */
class MeshVertexCoords {
 public:
  MeshVertexCoords(const field::Field<uint8_t> &basins_mask,
                   const field::Field<uint8_t> &coasts_mask,
                   const field::Field<uint8_t> &land_mask) {
    this->_masks[BASINS] = basins_mask;
    this->_masks[COASTS] = coasts_mask;
    this->_masks[LAND] = land_mask;
  }
  MeshVertexCoords(field::Field<uint8_t> &&basins_mask,
                   field::Field<uint8_t> &&coasts_mask,
                   field::Field<uint8_t> &&land_mask) {
    this->_masks[BASINS] = std::move(basins_mask);
    this->_masks[COASTS] = std::move(coasts_mask);
    this->_masks[LAND] = std::move(land_mask);
  }

  // Ensure this class is never copied
  MeshVertexCoords(const MeshVertexCoords &other) = delete;
  MeshVertexCoords &operator=(const MeshVertexCoords &) = delete;

  // But we want it to be moved
  MeshVertexCoords(MeshVertexCoords &&) = default;

  virtual ~MeshVertexCoords() = default;

  [[nodiscard]] virtual unsigned int nlon() const = 0;
  [[nodiscard]] virtual unsigned int nlat() const = 0;
  [[nodiscard]] virtual unsigned int nlev() const = 0;

  [[nodiscard]] virtual unsigned int ncells() const {
    return (this->nlon() - 1) * (this->nlat() - 1) * (this->nlev() - 1);
  }

  [[nodiscard]] virtual unsigned int nvertices() const {
    return this->nlon() * this->nlat() * this->nlev();
  }

  [[nodiscard]] inline const field::Field<uint8_t> &get_mask(
      const MASK_TYPE mask_type) const {
    return this->_masks[mask_type];
  }

  [[nodiscard]] virtual linecoords get_horizontal_grid_line(
      unsigned int line_index, unsigned int column_index) const = 0;
  [[nodiscard]] virtual linecoords get_vertical_grid_line(
      unsigned int line_index, unsigned int column_index) const = 0;
  [[nodiscard]] virtual linecoords get_column_grid_line(
      unsigned int vline_index, unsigned int hline_index) const = 0;

 protected:
  field::Field<uint8_t> _masks[3];
};

class RectilinearMeshVertexCoords : public MeshVertexCoords {
 public:
  RectilinearMeshVertexCoords(const std::vector<double> &lon2m,
                              const std::vector<double> &lat2m,
                              const std::vector<double> &nav_lev,
                              const field::Field<uint8_t> &basins_mask,
                              const field::Field<uint8_t> &coasts_mask,
                              const field::Field<uint8_t> &land_mask)
      : MeshVertexCoords(basins_mask, coasts_mask, land_mask),
        _lon2m(lon2m),
        _lat2m(lat2m),
        _nav_lev(nav_lev) {}

  RectilinearMeshVertexCoords(std::vector<double> &&lon2m,
                              std::vector<double> &&lat2m,
                              std::vector<double> &&nav_lev,
                              const field::Field<uint8_t> &basins_mask,
                              const field::Field<uint8_t> &coasts_mask,
                              const field::Field<uint8_t> &land_mask)
      : MeshVertexCoords(basins_mask, coasts_mask, land_mask),
        _lon2m(std::move(lon2m)),
        _lat2m(std::move(lat2m)),
        _nav_lev(std::move(nav_lev)) {}

  RectilinearMeshVertexCoords(const std::vector<double> &lon2m,
                              const std::vector<double> &lat2m,
                              const std::vector<double> &nav_lev,
                              field::Field<uint8_t> &&basins_mask,
                              field::Field<uint8_t> &&coasts_mask,
                              field::Field<uint8_t> &&land_mask)
      : MeshVertexCoords(std::move(basins_mask), std::move(coasts_mask),
                         std::move(land_mask)),
        _lon2m(lon2m),
        _lat2m(lat2m),
        _nav_lev(nav_lev) {}

  RectilinearMeshVertexCoords(std::vector<double> &&lon2m,
                              std::vector<double> &&lat2m,
                              std::vector<double> &&nav_lev,
                              field::Field<uint8_t> &&basins_mask,
                              field::Field<uint8_t> &&coasts_mask,
                              field::Field<uint8_t> &&land_mask)
      : MeshVertexCoords(std::move(basins_mask), std::move(coasts_mask),
                         std::move(land_mask)),
        _lon2m(std::move(lon2m)),
        _lat2m(std::move(lat2m)),
        _nav_lev(std::move(nav_lev)) {}

  ~RectilinearMeshVertexCoords() override = default;

  static std::unique_ptr<RectilinearMeshVertexCoords>
  create_rectilinear_mesh_vertex_coords(const std::string &file_name);

  [[nodiscard]] inline unsigned int nlon() const override {
    return this->_lon2m.size();
  }
  [[nodiscard]] inline unsigned int nlat() const override {
    return this->_lat2m.size();
  }
  [[nodiscard]] inline unsigned int nlev() const override {
    return this->_nav_lev.size();
  }

  [[nodiscard]] linecoords get_horizontal_grid_line(
      unsigned int line_index, unsigned int column_index) const override;
  [[nodiscard]] linecoords get_vertical_grid_line(
      unsigned int line_index, unsigned int column_index) const override;
  [[nodiscard]] linecoords get_column_grid_line(
      unsigned int vline_index, unsigned int hline_index) const override;

  void write(const std::string &file_path);

 private:
  std::vector<double> _lon2m, _lat2m, _nav_lev;
};

class Struct2DMeshVertexCoords : public MeshVertexCoords {
 public:
  Struct2DMeshVertexCoords(const field::Field<double> &lon2m,
                           const field::Field<double> &lat2m,
                           const std::vector<double> &nav_lev,
                           const field::Field<uint8_t> &basins_mask,
                           const field::Field<uint8_t> &coasts_mask,
                           const field::Field<uint8_t> &land_mask)
      : MeshVertexCoords(basins_mask, coasts_mask, land_mask),
        _lon2m(lon2m),
        _lat2m(lat2m),
        _nav_lev(nav_lev) {}

  Struct2DMeshVertexCoords(field::Field<double> &&lon2m,
                           field::Field<double> &&lat2m,
                           std::vector<double> &&nav_lev,
                           const field::Field<uint8_t> &basins_mask,
                           const field::Field<uint8_t> &coasts_mask,
                           const field::Field<uint8_t> &land_mask)
      : MeshVertexCoords(basins_mask, coasts_mask, land_mask),
        _lon2m(std::move(lon2m)),
        _lat2m(std::move(lat2m)),
        _nav_lev(std::move(nav_lev)) {}

  Struct2DMeshVertexCoords(const field::Field<double> &lon2m,
                           const field::Field<double> &lat2m,
                           const std::vector<double> &nav_lev,
                           field::Field<uint8_t> &&basins_mask,
                           field::Field<uint8_t> &&coasts_mask,
                           field::Field<uint8_t> &&land_mask)
      : MeshVertexCoords(std::move(basins_mask), std::move(coasts_mask),
                         std::move(land_mask)),
        _lon2m(lon2m),
        _lat2m(lat2m),
        _nav_lev(nav_lev) {}

  Struct2DMeshVertexCoords(field::Field<double> &&lon2m,
                           field::Field<double> &&lat2m,
                           std::vector<double> &&nav_lev,
                           field::Field<uint8_t> &&basins_mask,
                           field::Field<uint8_t> &&coasts_mask,
                           field::Field<uint8_t> &&land_mask)
      : MeshVertexCoords(std::move(basins_mask), std::move(coasts_mask),
                         std::move(land_mask)),
        _lon2m(std::move(lon2m)),
        _lat2m(std::move(lat2m)),
        _nav_lev(std::move(nav_lev)) {}

  ~Struct2DMeshVertexCoords() override = default;

  static std::unique_ptr<Struct2DMeshVertexCoords>
  create_struct2d_mesh_vertex_coords(const std::string &file_name);

  [[nodiscard]] inline unsigned int nlon() const override {
    return this->_lon2m.get_m();
  }
  [[nodiscard]] inline unsigned int nlat() const override {
    return this->_lon2m.get_n();
  }
  [[nodiscard]] inline unsigned int nlev() const override {
    return this->_nav_lev.size();
  }

  [[nodiscard]] linecoords get_horizontal_grid_line(
      unsigned int line_index, unsigned int column_index) const override;
  [[nodiscard]] linecoords get_vertical_grid_line(
      unsigned int line_index, unsigned int column_index) const override;
  [[nodiscard]] linecoords get_column_grid_line(
      unsigned int vline_index, unsigned int hline_index) const override;

 private:
  field::Field<double> _lon2m, _lat2m;
  std::vector<double> _nav_lev;
};

class GeoCoords {
 public:
  explicit GeoCoords(const std::string &file_name);

  [[nodiscard]] inline unsigned int nlon_cells() const { return this->_nlon; }
  [[nodiscard]] inline unsigned int nlat_cells() const { return this->_nlat; }
  [[nodiscard]] inline unsigned int nlev_cells() const { return this->_nlev; }
  [[nodiscard]] inline unsigned int ncells() const {
    return this->_nlon * this->_nlat * this->_nlev;
  }

  [[nodiscard]] const field::Field<double> &e1() const { return this->_e1; }
  [[nodiscard]] const field::Field<double> &e2() const { return this->_e2; }
  [[nodiscard]] const field::Field<double> &e3() const { return this->_e3; }

 private:
  unsigned int _nlon, _nlat, _nlev;

  field::Field<double> _e1, _e2, _e3;
};

class MeshData {
 public:
  MeshData(std::unique_ptr<MeshVertexCoords> vcoords,
           std::unique_ptr<GeoCoords> gcoords)
      : vertex_coords(std::move(vcoords)), geo_coords(std::move(gcoords)) {
    assert(vertex_coords->nlon() == geo_coords->nlon_cells() + 1);
    assert(vertex_coords->nlat() == geo_coords->nlat_cells() + 1);
    assert(vertex_coords->nlev() == geo_coords->nlev_cells() + 1);
  }

  [[nodiscard]] inline const field::Field<uint8_t> &get_mask(
      const MASK_TYPE mask_type) const {
    return this->vertex_coords->get_mask(mask_type);
  }

  [[nodiscard]] unsigned int nlon_vertices() const {
    return this->vertex_coords->nlon();
  }
  [[nodiscard]] unsigned int nlat_vertices() const {
    return this->vertex_coords->nlat();
  }
  [[nodiscard]] unsigned int nlev_vertices() const {
    return this->vertex_coords->nlev();
  }
  [[nodiscard]] unsigned int nvertices() const {
    return this->vertex_coords->nvertices();
  }

  [[nodiscard]] inline linecoords get_horizontal_grid_line(
      unsigned int line_index, unsigned int column_index) const {
    return this->vertex_coords->get_horizontal_grid_line(line_index,
                                                         column_index);
  }
  [[nodiscard]] inline linecoords get_vertical_grid_line(
      unsigned int line_index, unsigned int column_index) const {
    return this->vertex_coords->get_vertical_grid_line(line_index,
                                                       column_index);
  }
  [[nodiscard]] linecoords get_column_grid_line(
      unsigned int vline_index, unsigned int hline_index) const {
    return this->vertex_coords->get_column_grid_line(vline_index, hline_index);
  }

  [[nodiscard]] unsigned int nlon_cells() const {
    return this->geo_coords->nlon_cells();
  }
  [[nodiscard]] unsigned int nlat_cells() const {
    return this->geo_coords->nlat_cells();
  }
  [[nodiscard]] unsigned int nlev_cells() const {
    return this->geo_coords->nlev_cells();
  }
  [[nodiscard]] unsigned int ncells() const {
    return this->geo_coords->ncells();
  }

  [[nodiscard]] inline const field::Field<double> &e1() const {
    return this->geo_coords->e1();
  }
  [[nodiscard]] inline const field::Field<double> &e2() const {
    return this->geo_coords->e2();
  }
  [[nodiscard]] inline const field::Field<double> &e3() const {
    return this->geo_coords->e3();
  }

 private:
  std::unique_ptr<MeshVertexCoords> vertex_coords;
  std::unique_ptr<GeoCoords> geo_coords;
};

/* Simulation MESH

        Stores information on the mesh and the projection used.
        This class does not hold the data, but it is only a pointer to the
        files on the disk that actually describe the mesh. If you want to access
        to the real data, use the "get_mesh_data" method for obtaining a
        MeshData object

*/
class Mesh {
 public:
  Mesh(std::string name, std::string meshfile, std::string meshmask)
      : _name(std::move(name)),
        _meshfile(std::move(meshfile)),
        _meshmask(std::move(meshmask)) {}
  ~Mesh() = default;

  [[nodiscard]] std::string get_name() const { return this->_name; }
  [[nodiscard]] std::string get_meshfile() const { return this->_meshfile; }
  [[nodiscard]] std::string get_meshmask() const { return this->_meshmask; }

 private:
  std::string _name, _meshfile, _meshmask;
};

/* Simulation CLASS

        Stores the necessary information from the Simulation master file and
   helps in reading/writing the variables.

*/
class Simulation {
 public:
  // Constructors
  explicit Simulation(const char *fname) : _ogsfile(fname) { readMainFile(); }
  explicit Simulation(std::string fname) : _ogsfile(std::move(fname)) {
    readMainFile();
  }

  // Destructor
  ~Simulation() = default;

  [[nodiscard]] inline int var_n(VAR_TYPE vrt) const {
    return this->_vars[vrt].get_nvars();
  };
  [[nodiscard]] inline std::string var_name(VAR_TYPE vrt,
                                            unsigned int jj) const {
    return this->_vars[vrt].get_name(jj);
  };
  [[nodiscard]] inline std::string var_vname(VAR_TYPE vrt,
                                             unsigned int jj) const {
    return this->_vars[vrt].get_vname(jj);
  };

  [[nodiscard]] std::optional<std::pair<VAR_TYPE, unsigned int>>
  get_var_indices(const std::string &var_name) const;

  // Return the path of the file that contains the data for the t timestep of
  // the jj variable of type vrt
  [[nodiscard]] std::string var_path(VAR_TYPE vrt, unsigned int jj,
                                     unsigned int t) const;

  [[nodiscard]] inline unsigned int n_meshes() const {
    return this->meshes.size();
  }

  [[nodiscard]] unsigned int get_mesh_index(const std::string &proj_name) const;

  [[nodiscard]] inline unsigned int n_projections() const {
    return this->n_meshes();
  }

  std::string projection(unsigned int i) { return meshes.at(i).get_name(); };

  [[nodiscard]] inline unsigned int ntsteps() const { return this->_ntsteps; }
  [[nodiscard]] inline std::string datetime(unsigned int t) const {
    return this->_datetime.at(t);
  }

  void load_mesh(unsigned int mesh_index);

  [[nodiscard]] std::shared_ptr<const MeshData> loaded_mesh() const {
    if (this->mesh_loaded == -1)
      throw std::logic_error(
          "No mesh has been loaded inside this object at the moment");
    return this->mesh_data;
  }

  void print() const;

 private:
  int readMainFile();

  std::string _ogsfile, _wrkdir;
  OGS_VAR _vars[4];  // Variables information
  int _ntsteps = -1;
  std::vector<std::string> _datetime;
  std::vector<Mesh> meshes;

  MESH_TYPE mesh_type = RECTILINEAR;

  int mesh_loaded = -1;
  std::shared_ptr<MeshData> mesh_data = nullptr;
};
}  // namespace OGS
#endif
