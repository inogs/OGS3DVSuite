/*=========================================================================

  Module:    Projection

  Map projection conversions using the proj API.

  Copyright (c) 2019 Arnau Miro, OGS
  All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef PROJECTIONS_H
#define PROJECTIONS_H

#include <proj.h>
#include "OGS/geometry.h"

#include <memory>
#include <map>
#include <string>

#define PROJ_ERR 0
#define PROJ_OK 1

namespace OGS::PROJ {
inline constexpr double rad2deg(const double a) {
  return (180.0 / 3.14159265359) * a;
}
inline constexpr double deg2rad(const double a) {
  return (3.14159265359 / 180.0) * a;
}


// A unique pointer that points to a projection of PROJ
typedef std::unique_ptr<PJ, decltype(&proj_destroy)> PJpointer;


// When an invalid projection is created, we raise the following error
class InvalidProjection : std::invalid_argument {
 public:
  InvalidProjection(std::string from_, std::string to_)
      : from(std::move(from_)),
        to(std::move(to_)),
        std::invalid_argument(create_error_message(from_, to_)) {}

  const std::string from, to;

 private:
  static std::string create_error_message(const std::string &from_,
                                          const std::string &to_) {
    std::ostringstream error_msg;
    error_msg << "Invalid projection specified!" << std::endl
              << "FROM:  " << from_ << std::endl
              << "TO:  " << to_ << std::endl;
    return error_msg.str();
  }
};

// map containing all the default projections from the Suite
const std::map<const std::string, const std::string> PROJS = {
    {"degrees",
     "+ellps=WGS84 +a=57.29577951308232 +proj=eqc +lon_0=0.0 "
     "+no_defs"},  // Plate Carree
    {"mercator",
     "+ellps=WGS84 +proj=merc +lon_0=0.0 +x_0=989634.3811336625 "
     "+y_0=-3512473.95569 +units=m +no_defs"},  // Mercator centered
                                                // on MED
    {"cylindrical",
     "+datum=WGS84 +ellps=WGS84 +proj=eqc +lat_ts=0 +lat_0=0 "
     "+lon_0=0 +x_0=0 +y_0=0 +units=m +no_defs +no_defs"},  // Cylindrical
                                                            // using
                                                            // EPSG:4087
    {"google",
     "+a=6378137.0 +b=6378137.0 +nadgrids=@null +proj=merc "
     "+lon_0=0.0 +x_0=0.0 +y_0=0.0 +units=m +no_defs"},  // Google
                                                         // Mercator
    {"mollweide", "+a=6378137.0 +proj=moll +lon_0=0 +no_defs"},
    {"orthographic", "+ellps=WGS84 +proj=ortho +lon_0=0.0 +lat_0=0.0 +no_defs"},
    {"robinson", "+a=6378137.0 +proj=robin +lon_0=0 +no_defs"},
    {"satellite",
     "+a=6378137.0 +proj=nsper +lon_0=17.5 +lat_0=36.4 +h=6779000 "
     "+x_0=0 +y_0=0 +units=m +no_defs"},  // Nearside Perspective
                                          // centered on MED at ISS
                                          // altitude
    {"eckert iv", "+a=6378137.0 +proj=eck4 +lon_0=0 +no_defs"},
    {"equal earth", "+ellps=WGS84 +proj=eqearth +lon_0=0 +no_defs"},
    {"epsg 3857",
     "+ellps=WGS84 +a=6378137 +b=6378137 +nadgrids=@null +proj=merc "
     "+lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +wktext "
     "+ +no_defs +no_defs"}};

static std::map<const std::string, const std::string> get_projection_map() {
  return PROJS;
}


class Projection {
  /* PROJECTION

          Class interface with the proj library
  */
 public:
  Projection(
      const std::string &from, const std::string &to,
      const std::map<const std::string, const std::string> &projection_map)
      : projection(proj_create_crs_to_crs(
                       PJ_DEFAULT_CTX, projection_map.at(from).c_str(),
                       projection_map.at(to).c_str(), nullptr),
                   &proj_destroy) {
    if (projection == nullptr) throw InvalidProjection(from, to);
  }

  Projection(const std::string &from, const std::string &to)
      : Projection(from, to, PROJS) {}

  ~Projection() = default;

  OGS::Geom::Point<double> transform_point(double x, double y, double z) {
    PJ_COORD c, c_out;
    c.uvwt.u = x;
    c.uvwt.v = y;
    c.uvwt.w = z;
    c.uvwt.t = HUGE_VAL;
    c_out = proj_trans(&(*(this->projection)), PJ_FWD, c);

    return {c_out.uvwt.u, c_out.uvwt.v, c_out.uvwt.w};
  }

  Geom::Point<double> transform_point(double x, double y) {
    return this->transform_point(x, y, 0.);
  }

  /*inline int transform_point(const char *from, const char *to,
                             const double xy[], double out[]) {
    double x = xy[0], y = xy[1];
    int errcode = transform_point(from, to, x, y);
    out[0] = x;
    out[1] = y;
    return errcode;
  }

  inline int transform_point(const char *from, const std::string &to,
                             double xy[], double out[]) {
    double x = xy[0], y = xy[1];
    int errcode = transform_point(from, to, x, y);
    out[0] = x;
    out[1] = y;
    return errcode;
  }

  inline int transform_point(const std::string &from, const char *to,
                             double xy[], double out[]) {
    double x = xy[0], y = xy[1];
    int errcode = transform_point(from, to, x, y);
    out[0] = x;
    out[1] = y;
    return errcode;
  }

  inline int transform_point(const std::string &from, const std::string &to,
                             double xy[], double out[]) {
    double x = xy[0], y = xy[1];
    int errcode = transform_point(from, to, x, y);
    out[0] = x;
    out[1] = y;
    return errcode;
  }

  inline int proj_transform_point(const char *from, const char *to, double *x,
                                  double *y) {
    //if (from == NULL || to == NULL) return PROJ_ERR;

    source = pj_init_plus(from);
    target = pj_init_plus(to);
    alloc = true;
    if (source == NULL || target == NULL) return PROJ_ERR;

    int retval =
        pj_transform(source, target, 1, 1, x, y, NULL);  // 0 on success
    free();
    return (retval == 0) ? PROJ_OK : PROJ_ERR;
  }*/

 private:
  PJpointer projection;
};
}  // namespace OGS::PROJ

#endif
