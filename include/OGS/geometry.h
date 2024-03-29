/*
    GEOMETRY

        Definition of a generic geometric tools such as polygons.

        Arnau Miro (OGS) (c) 2019
*/

#ifndef Geometry_h
#define Geometry_h

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace OGS::Geom {

// Preamble definition of all classes
template <class T>
class Vector;  // Class describing a vector
template <class T>
class Ball;  // Ball bounding box
template <class T>
class Polygon;

// Classes definition

template <class T>
class Point {
 public:
  // Constructors and destructors
  inline Point() { set(0, 0, 0); }
  inline Point(const T x, const T y, const T z) { set(x, y, z); }
  inline explicit Point(const T *v) { set(v[0], v[1], v[2]); }
  inline Point(const Point<T> &pp) { set(pp[0], pp[1], pp[2]); }

  template <class U>
  inline explicit Point(const Point<U> &pp)
      : p({static_cast<T>(pp.p[0]), static_cast<T>(pp.p[1])},
          static_cast<T>(pp.p[2])) {}

  inline ~Point() = default;

  // Functions
  inline void set(const T x, const T y, const T z) {
    this->p[0] = x;
    this->p[1] = y;
    this->p[2] = z;
  }
  inline void set(const T *v) { set(v[0], v[1], v[2]); }
  inline void set(const Point<T> &pp) { set(pp[0], pp[1], pp[2]); }
  inline T *data() { return p; }
  inline T x() { return p[0]; }
  inline T y() { return p[1]; }
  inline T z() { return p[2]; }
  inline void print() const { std::printf("[ %f %f %f ]", p[0], p[1], p[2]); }

  inline T dist(const Point<T> &pp) const {
    Vector<T> d = (*this) - pp;
    return d.norm();
  }
  inline T dist2(const Point<T> &pp) const {
    Vector<T> d = (*this) - pp;
    return d.norm2();
  }

  /* ISLEFT

          Tests if a point is Left|On|Right of an infinite line.

          Input:  two points P0, P1; defining a line
          Return: >0 for P2 left of the line through P0 and P1
                  =0 for P2  on the line
                  <0 for P2  right of the line

          from: http://geomalgorithms.com/a03-_inclusion.html

          Copyright 2001 softSurfer, 2012 Dan Sunday
          This code may be freely used and modified for any purpose
          providing that this copyright notice is included with it.
          SoftSurfer makes no warranty for this code, and cannot be held
          liable for any real or imagined damage resulting from its use.
          Users of this code must verify correctness for their application.
  */
  inline T isLeft(const Point<T> P0, const Point<T> P1) const {
    return ((P1[0] - P0[0]) * (p[1] - P0[1]) -
            (p[0] - P0[0]) * (P1[1] - P0[1]));
  }

  // Operators
  inline T operator[](int i) const { return p[i]; }
  inline T &operator[](int i) { return p[i]; }
  inline Point<T> &operator=(const Point<T> &pp) {
    const double v0 = pp[0];
    const double v1 = pp[1];
    const double v2 = pp[2];
    set(v0, v1, v2);
    return (*this);
  }

  template <class U>
  inline Point<T> &operator=(const Point<U> &pp) {
    this->p[0] = static_cast<T>(pp.p[0]);
    this->p[1] = static_cast<T>(pp.p[1]);
    this->p[2] = static_cast<T>(pp.p[2]);
    return (*this);
  }

  inline Point<T> operator+(const Vector<T> &vv) const {
    return Point<T>(p[0] + vv[0], p[1] + vv[1], p[2] + vv[2]);
  }
  inline Point<T> operator-(const Vector<T> &vv) const {
    return Point<T>(p[0] - vv[0], p[1] - vv[1], p[2] - vv[2]);
  }
  inline Vector<T> operator-(const Point<T> &pp) const {
    return Vector<T>(p[0] - pp[0], p[1] - pp[1], p[2] - pp[2]);
  }
  inline bool operator==(const Point<T> &v) const {
    return ((p[0] == v[0]) && (p[1] == v[1]) && (p[2] == v[2]));
  }
  inline bool operator!=(const Point<T> &v) const {
    return ((!(*this)) == (v));
  }

 private:
  T p[3];
};



template <class T>
class Vector {
 public:
  // Constructors and destructors
  inline Vector() { set(0, 0, 0); }
  inline Vector(const T x, const T y, const T z) { set(x, y, z); }
  inline explicit Vector(const T *p) { set(p[0], p[1], p[2]); }
  inline Vector(const Vector<T> &vv) { set(vv[0], vv[1], vv[2]); }
  inline ~Vector() = default;

  // Functions
  inline void set(const T x, const T y, const T z) {
    v[0] = x;
    v[1] = y;
    v[2] = z;
  }
  inline void set(const T *p) { set(p[0], p[1], p[2]); }
  inline void set(const Vector<T> &vv) { set(vv[0], vv[1], vv[2]); }
  inline T *data() { return v; }
  inline T x() { return v[0]; }
  inline T y() { return v[1]; }
  inline T z() { return v[2]; }

  inline T dot(const Vector<T> &vv) const {
    return ((T)(v[0] * vv[0] + v[1] * vv[1] + v[2] * vv[2]));
  }
  inline Vector<T> cross(const Vector<T> &vv) const {
    return (Vector<T>(v[1] * vv[2] - v[2] * vv[1], -v[0] * vv[2] + v[2] * vv[0],
                      v[0] * vv[1] - v[1] * vv[0]));
  }
  inline T norm() const { return ((T)(std::sqrt(norm2()))); }
  inline T norm2() const { return (dot((*this))); }
  inline void print() const { std::printf("( %f %f %f )", v[0], v[1], v[2]); }

  // Operators
  inline T operator[](int i) const { return v[i]; }
  inline T &operator[](int i) { return v[i]; }
  inline Vector<T> &operator=(const Vector<T> &vv) {
    const double p0 = vv[0];
    const double p1 = vv[1];
    const double p2 = vv[2];

    set(p0, p1, p2);
    return (*this);
  }
  inline Vector<T> operator+(const Vector<T> &vv) const {
    return Vector<T>(v[0] + vv[0], v[1] + vv[1], v[2] + vv[2]);
  }
  inline Vector<T> operator-(const Vector<T> &vv) const {
    return Vector<T>(v[0] - vv[0], v[1] - vv[1], v[2] - vv[2]);
  }
  inline Vector<T> operator*(const T a) const {
    return Vector<T>(v[0] * a, v[1] * a, v[2] * a);
  }
  inline Vector<T> operator/(const T a) const {
    return Vector<T>(v[0] / a, v[1] / a, v[2] / a);
  }
  inline T operator*(const Vector<T> &vv) const { return dot(vv); }
  inline Vector<T> operator^(const Vector<T> &vv) const { return cross(vv); }

 private:
  T v[3];
};

template <class T>
Vector<T> operator*(const T a, const Vector<T> &v) {
  return Vector<T>(v[0] * a, v[1] * a, v[2] * a);
}



template <class T>
class Ball {
 public:
  // Constructor and destructors
  inline Ball() { set(Point<T>(0, 0, 0), 0); }
  inline Ball(const T r, const Point<T> &p) { set(r, p); }
  inline Ball(const Point<T> &p, const T r) { set(r, p); }
  inline explicit Ball(const Polygon<T> &p) { fastBall(p); }
  inline Ball(const Ball<T> &b) { set(b.get_center(), b.get_radius()); }
  inline ~Ball() = default;

  // Functions
  inline void set(const T r, const Point<T> &p) {
    set_center(p);
    set_radius(r);
  }
  inline void set(const Point<T> &p, const T r) {
    set_center(p);
    set_radius(r);
  }
  inline void set_center(const Point<T> &p) { center = p; }
  inline void set_radius(const T r) { radius = r; }
  inline Point<T> get_center() const { return center; }
  inline T get_radius() const { return radius; }

  inline void print() const {
    std::printf("center = ");
    center.print();
    printf(" radius = %f\n", radius);
  }

  [[nodiscard]] inline bool isempty() const { return (radius == 0); }
  inline bool contains(const Point<T> &p) const {
    return (!isempty() && p.dist(center) < static_cast<double>(radius));
  }
  inline bool isdisjoint(const Ball<T> &b) const {
    return b.get_center().dist(this->center) <
           static_cast<double>(radius + b.get_radius());
  }

  inline void fastBall(const Polygon<T> &p);

  // Operators
  inline Ball<T> &operator=(const Ball<T> &b) {
    set(b.get_center(), b.get_radius());
    return (*this);
  }
  inline bool operator==(const Ball<T> &b) const {
    return (center == b.get_center() && radius == b.get_radius());
  }
  inline bool operator>(const Point<T> &p) const {
    return (contains(p));
  }  // Return true if Point inside Ball
  inline bool operator<(const Point<T> &p) const {
    return (!contains(p));
  }  // Return true if Point outside Ball

 private:
  Point<T> center;
  T radius;
};

/* FASTBALL

        Get a fast approximation for the 2D bounding ball
        (based on the algorithm given by [Jack Ritter, 1990]).

        Input:  A polygon
        Output: Nothing, sets the ball class

        from: http://geomalgorithms.com/a08-_containers.html

        Copyright 2001 softSurfer, 2012 Dan Sunday
        This code may be freely used and modified for any purpose
        providing that this copyright notice is included with it.
        SoftSurfer makes no warranty for this code, and cannot be held
        liable for any real or imagined damage resulting from its use.
        Users of this code must verify correctness for their application.
*/
template <class T>
inline void Ball<T>::fastBall(const Polygon<T> &p) {
  T xmin, xmax, ymin, ymax;        // bounding box extremes
  int Pxmin, Pxmax, Pymin, Pymax;  // index of  P[] at box extreme

  // Find a large diameter to start with
  // first get the bounding box and P[] extreme points for it
  xmin = xmax = p[0][0];
  ymin = ymax = p[0][1];
  Pxmin = Pxmax = Pymin = Pymax = 0;

  for (int ii = 1; ii < p.get_npoints(); ++ii) {
    if (p[ii][0] < xmin) {
      xmin = p[ii][0];
      Pxmin = ii;
    } else if (p[ii][0] > xmax) {
      xmax = p[ii][0];
      Pxmax = ii;
    }
    if (p[ii][1] < ymin) {
      ymin = p[ii][1];
      Pymin = ii;
    } else if (p[ii][1] > ymax) {
      ymax = p[ii][1];
      Pymax = ii;
    }
  }

  // Select the largest extent as an initial diameter for the  ball
  Point<T> C;
  Vector<T> dPx = p[Pxmax] - p[Pxmin], dPy = p[Pymax] - p[Pymin];
  T rad2, dx2 = dPx.norm2(), dy2 = dPy.norm2();

  if (dx2 >= dy2) {  // x direction is largest extent
    C = p[Pxmin] + (dPx / 2.);
    rad2 = p[Pxmax].dist2(C);
  } else {
    C = p[Pymin] + (dPy / 2.);
    rad2 = p[Pymax].dist2(C);
  }
  T rad = std::sqrt(rad2);

  // Now check that all points points[i] are in the ball
  // and if not, expand the ball just enough to include them
  Vector<T> dP;
  T dist, dist2;

  for (int ii = 0; ii < p.get_npoints(); ++ii) {
    dP = p[ii] - C;
    dist2 = dP.norm2();
    if (dist2 <= rad2) continue;  // points[i] is inside the ball already
    // points[i] not in ball, so expand ball  to include it
    dist = sqrt(dist2);
    rad = (rad + dist) / 2.;
    rad2 = rad * rad;                    // enlarge radius just enough
    C = C + ((dist - rad) / dist) * dP;  // shift Center toward points[i]
  }

  // Set ball parameters
  this->set_center(C);
  this->set_radius(rad);
}

template <class T>
class Polygon {
 public:
  // Constructors and destructors
  inline Polygon() : points(0), bbox() {}
  inline Polygon(const int nn, const Point<T> *v) : points(nn), bbox() {
    for (unsigned int i = 0; i < nn; ++i) this->points[i] = v[i];
  }
  inline explicit Polygon(const std::vector<Point<T>> &v) : points(v), bbox() {}
  inline explicit Polygon(const std::vector<Point<T>> &&v)
      : points(std::move(v)), bbox() {}

  inline ~Polygon() = default;

  void clear() { this->points.resize(0); }

  void set_points(const unsigned int nn, const Point<T> *v) {
    this->points.resize(nn);
    for (unsigned int i = 0; i < nn; ++i) this->points[i] = v[i];
    this->bbox_set = false;
  }

  void set_points(const std::vector<Point<T>> &v) {
    this->points = v;
    this->bbox_set = false;
  }

  void set_points(const std::vector<Point<T>> &&v) {
    this->points = v;
    this->bbox_set = false;
  }

  inline const std::vector<Point<T>> &get_points() const { return points; }

  inline Point<T> get_point(const unsigned int i) const { return points.at(i); }

  [[nodiscard]] inline size_t get_npoints() const { return points.size(); }

  inline Ball<T> get_bbox() {
    if (!this->bbox_set) {
      this->bbox.fastBall(*this);
      this->bbox_set = true;
    }
    return bbox;
  }

  [[nodiscard]] inline bool is_empty() const { return points.empty(); }

  inline bool is_inside(const Point<T> &v) {
    if (!this->bbox_set) {
      this->bbox.fastBall(*this);
      this->bbox_set = true;
    }
    if (this->bbox.contains(v))
      return (wn_PinPoly(v) != 0);
    else
      return false;
  }
  inline bool is_inside_cn(const Point<T> &v) const {
    return cn_PinPoly(v) == 1;
  }
  inline bool is_inside_wn(const Point<T> &v) const {
    return wn_PinPoly(v) != 0;
  }

  template <class U>
  inline Polygon<U> convert() {
    std::vector<Point<U>> poly_data(this->points);
    Polygon<U> pp(std::move(poly_data));
    return pp;
  }

  inline void print() const {
    unsigned int i = 0;
    for (auto &p : this->points) {
      std::cout << "Point " << i++ << ":   " << p.print() << std::endl;
    }
  }

  inline int cn_PinPoly(
      const Point<T> &P) const;  // Return:  0 = outside, 1 = inside
  inline int wn_PinPoly(
      const Point<T> &P) const;  // Return:  =0 only when P is outside

  // Operators
  inline Point<T> operator[](int i) const {
    return (i >= 0) ? points.at(i) : points.at(this->points.size() + i);
  }
  inline Point<T> &operator[](int i) {
    return (i >= 0) ? points.at(i) : points.at(this->points.size() + i);
  }
  inline Point<T> operator[](unsigned int i) const { return points.at(i); }
  inline Point<T> &operator[](unsigned int i) { return points.at(i); }
  inline Point<T> operator[](size_t i) const { return points.at(i); }
  inline Point<T> &operator[](size_t i) { return points.at(i); }

  inline bool operator==(const Polygon<T> &pp) const;
  inline bool operator!=(const Polygon<T> &pp) const {
    return (!(*this) == pp);
  }
  inline bool operator>(const Point<T> &v) {
    return (is_inside(v));
  }  // Return true if Point inside Polygon
  inline bool operator<(const Point<T> &v) {
    return (!is_inside(v));
  }  // Return true if Point outside Polygon

 private:
  std::vector<Point<T>> points;
  Ball<T> bbox;
  bool bbox_set = false;
};

template <class T>
inline bool Polygon<T>::operator==(const Polygon<T> &pp) const {
  // Check if polygons have the same number of points
  if (this->n != pp.get_npoints()) return false;

  // Loop points and check whether they are equal
  for (int ii = 0; ii < this->n; ++ii) {
    if (this->points[ii] != pp[ii]) return false;
  }

  return true;
}

/* CN_PINPOLY

        Crossing number test for a point in a polygon.

        Input:   P = a point,
        Return:  0 = outside, 1 = inside

        This code is patterned after [Franklin, 2000]
        from: http://geomalgorithms.com/a03-_inclusion.html

        Copyright 2001 softSurfer, 2012 Dan Sunday
        This code may be freely used and modified for any purpose
        providing that this copyright notice is included with it.
        SoftSurfer makes no warranty for this code, and cannot be held
        liable for any real or imagined damage resulting from its use.
        Users of this code must verify correctness for their application.
*/
template <class T>
inline int Polygon<T>::cn_PinPoly(const Point<T> &P) const {
  int cn = 0;  // The crossing number counter

  // Loop through all edges of the Polygon
  for (int ii = 0; ii < this->points.size(); ++ii) {
    if (((this->points[ii][1] <= P[1]) &&
         (this->points[ii + 1][1] > P[1]))  // an upward crossing
        || ((this->points[ii][1] > P[1]) &&
            (this->points[ii + 1][1] <= P[1]))) {  // a downward crossing

      // Compute  the actual edge-ray intersect x-coordinate
      auto vt =
          static_cast<double>((P[1] - this->points[ii][1]) /
                              (this->points[ii + 1][1] - this->points[ii][1]));
      if (P[0] <
          this->points[ii][0] + vt * (this->points[ii + 1][0] -
                                      this->points[ii][0]))  // P.x < intersect
        ++cn;  // A valid crossing of y=P.y right of P.x
    }
  }
  return (cn & 1);  // 0 if even (out), and 1 if  odd (in)
}

/* WN_PINPOLY

        Winding number test for a point in a polygon.

        Input:   P = a point,
        Return:  wn = the winding number (=0 only when P is outside)

        from: http://geomalgorithms.com/a03-_inclusion.html

        Copyright 2001 softSurfer, 2012 Dan Sunday
        This code may be freely used and modified for any purpose
        providing that this copyright notice is included with it.
        SoftSurfer makes no warranty for this code, and cannot be held
        liable for any real or imagined damage resulting from its use.
        Users of this code must verify correctness for their application.
*/
template <class T>
inline int Polygon<T>::wn_PinPoly(const Point<T> &P) const {
  int wn = 0;  // The  winding number counter

  // Loop through all edges of the polygon
  for (int ii = 0; ii < this->points.size();
       ++ii) {                             // edge from V[i] to  V[i+1]
    if (this->points[ii][1] <= P[1]) {     // start y <= P.y
      if (this->points[ii + 1][1] > P[1])  // an upward crossing
        if (P.isLeft(this->points[ii], this->points[ii + 1]) >
            0)                              // P left of  edge
          ++wn;                             // have  a valid up intersect
    } else {                                // start y > P.y (no test needed)
      if (this->points[ii + 1][1] <= P[1])  // a downward crossing
        if (P.isLeft(this->points[ii], this->points[ii + 1]) <
            0)   // P left of  edge
          --wn;  // have  a valid down intersect
    }
  }
  return wn;
}

}  // namespace OGS::Geom

#endif
