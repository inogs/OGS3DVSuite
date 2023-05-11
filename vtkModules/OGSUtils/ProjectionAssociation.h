//
// Created by Stefano Piani on 5/10/23.
//

#ifndef OGS3DVSUITE_PROJECTIONASSOCIATION_H
#define OGS3DVSUITE_PROJECTIONASSOCIATION_H

#include <stdexcept>
#include <string>

class InvalidProjectionIndex : std::out_of_range {
 public:
  explicit InvalidProjectionIndex(const int i)
      : std::out_of_range("No projection associated to index " +
                          std::to_string(i)) {}
};


static std::string associate_projection(int i) {
  switch (i) {
    case 0:
      return "mercator";
    case 1:
      return "cylindrical";
    case 2:
      return "google";
    case 3:
      return "mollweide";
    case 4:
      return "orthographic";
    case 5:
      return "robinson";
    case 6:
      return "satellite";
    case 7:
      return "eckert iv";
    case 8:
      return "equal earth";
    case 9:
      return "epsg 3857";
    default:
      throw InvalidProjectionIndex(i);
  }
}

#endif  // OGS3DVSUITE_PROJECTIONASSOCIATION_H
