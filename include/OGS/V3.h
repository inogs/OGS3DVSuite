/*
    VECTOR 3D

        Defines a 3 dimensional vector composed of |x,y,z|.
        Allows operations between vectors.

        Let v1 and v2 be vectors of |x,y,z| and k a random value...

        Construction:
                > V3()      empty construction (sets to zero)
                > V3(x,y,z) build from a series of x,y,z points
                > V3( V3 )  build from an existing V3 vector

        Usage:
                > v1 = v2   sets the contents of v2 to v1
                > k  + v1   adds a V3 vector by k in all positions
                > v1 + k    adds a V3 vector by k in all positions
                > k  - v1   substracts a V3 vector by k in all positions
                > v1 - k    substracts a V3 vector by k in all positions
                > k  * v1   multiplies a V3 vector by k in all positions
                > v1 * k    multiplies a V3 vector by k in all positions
                > v1 / k    divides a V3 vector by k in all positions
                > v1 + v2   adds two V3 vectors
                > v1 - v2   subtracts two V3 vectors
                > v1 * v2   scalar product between two V3 vectors
                > v1 ^ v2   cross product between two V3 vectors
                > v1 += k   adds a V3 vector by k in all positions
                > v1 -= k   substracts a V3 vector by k in all positions
                > v1 *= k   multiplies a V3 vector by k in all positions
                > v1 /= k   divides a V3 vector by k in all positions
                > v1 += v2  adds two V3 vectors
                > v1 -= v2  subtracts two V3 vectors
                > v1 *= v2  point to point product between two V3 vectors
                > v1 == v2  checks equalty using the euclidean norm

        Functions:
                > norm2:    performs the square of the vector norm
                > print:    prints the vector in a nice format

        Arnau Miro (UPC-ESEIAAT) (c) 2018
*/

#ifndef V3_h
#define V3_h

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

namespace OGS::V3 {
class V3 {
  friend inline V3 operator+(double val, const V3 &v);
  friend inline V3 operator-(double val, const V3 &v);
  friend inline V3 operator*(double val, const V3 &v);

 public:
  // Constructor
  V3(double x, double y, double z) : val{x, y, z} {}
  V3() : V3(0., 0., 0.) {}
  V3(const V3 &v) = default;

  // Destructor
  ~V3() = default;

  // Functions
  [[nodiscard]] double norm2() const {
    return val[0] * val[0] + val[1] * val[1] + val[2] * val[2];
  }
  void print() const {
    std::cout << "|" << val[0] << ", " << val[1] << ", " << val[2] << "|"
              << std::endl;
  }

  // Operators
  V3 &operator=(const V3 &v) = default;

  // Operations with double
  inline V3 operator+(double v) const {
    return {val[0] + v, val[1] + v, val[2] + v};
  }
  inline V3 operator-(double v) const {
    return {val[0] - v, val[1] - v, val[2] - v};
  }
  inline V3 operator*(double v) const {
    return {val[0] * v, val[1] * v, val[2] * v};
  }
  inline V3 operator/(double v) const {
    return {val[0] / v, val[1] / v, val[2] / v};
  }

  // Inplace operations with double
  inline void operator+=(double v) {
    val[0] += v;
    val[1] += v;
    val[2] += v;
  }
  inline void operator-=(double v) {
    val[0] -= v;
    val[1] -= v;
    val[2] -= v;
  }
  inline void operator*=(double v) {
    val[0] *= v;
    val[1] *= v;
    val[2] *= v;
  }
  inline void operator/=(double v) {
    val[0] /= v;
    val[1] /= v;
    val[2] /= v;
  }

  // Operators with other vectors
  inline V3 operator+(const V3 &v) const {
    return {val[0] + v[0], val[1] + v[1], val[2] + v[2]};
  }
  inline V3 operator-(const V3 &v) const {
    return {val[0] - v[0], val[1] - v[1], val[2] - v[2]};
  }
  inline double operator*(const V3 &v) const {
    return val[0] * v[0] + val[1] * v[1] + val[2] * v[2];
  }
  inline V3 operator^(const V3 &v) const {
    return {val[1] * v[2] - val[2] * v[1], -val[0] * v[2] + val[2] * v[0],
            val[0] * v[1] - val[1] * v[0]};
  }

  // Inplace operators with other vectors
  inline void operator+=(const V3 &v) {
    this->val[0] += v[0];
    this->val[1] += v[1];
    this->val[2] += v[2];
  }
  inline void operator-=(const V3 &v) {
    this->val[0] -= v[0];
    this->val[1] -= v[1];
    this->val[2] -= v[2];
  }
  inline void operator*=(const V3 &v) {
    this->val[0] *= v[0];
    this->val[1] *= v[1];
    this->val[2] *= v[2];
  }

  // Tolerance of the order of 10^-10
  inline bool operator==(const V3 &v) const {
    return ((*this - v).norm2() < 1e-10);
  }

  // We define that one vector is smaller than the other if the l2 norm of the
  // first vector is smaller than the other
  inline bool operator<(const V3 &v) const {
    return (this->norm2() < v.norm2());
  }

  inline double operator[](int i) const { return val[i]; }
  inline double &operator[](int i) { return val[i]; }

 private:
  double val[3];
};


// Operators
inline V3 operator+(double val, const V3 &v) {
  return {val + v[0], val + v[1], val + v[2]};
}
inline V3 operator-(double val, const V3 &v) {
  return {val - v[0], val - v[1], val - v[2]};
}
inline V3 operator*(double val, const V3 &v) {
  return {val * v[0], val * v[1], val * v[2]};
}


// V3v:  Vector of V3
class V3v {
 public:
  // Constructors and destructors
  V3v() : n(0), v(nullptr), alloc(false) {}
  explicit V3v(unsigned int nn);
  V3v(unsigned int nn, const float *val);
  V3v(unsigned int nn, const double *val);
  V3v(unsigned int nn, const V3 *val);
  V3v(const V3v &v);
  ~V3v();

  // Functions
  inline void size(unsigned int nn);
  inline void size(unsigned int nn, const V3 *val);
  inline void size(const V3v &v);
  inline void size(unsigned int nn, const float *val);
  inline void size(unsigned int nn, const double *val);

  [[nodiscard]] unsigned int len() const { return this->n; }
  inline void clear();
  inline V3 *data();

  template <typename T=double>
  T *toarray() const;

  inline bool isempty();

  // Operators
  inline V3v &operator=(const V3v &v);

  inline V3 operator[](int i) const { return (i >= 0) ? v[i] : v[n + i]; }
  inline V3 &operator[](int i) { return (i >= 0) ? v[i] : v[n + i]; }

  inline V3 operator[](unsigned int i) const { return v[i]; }
  inline V3 &operator[](unsigned int i) { return v[i]; }

  inline V3 operator[](size_t i) const { return v[i]; }
  inline V3 &operator[](size_t i) { return v[i]; }

  // Iterator
  class iterator {
   public:
    typedef V3 value_type;
    typedef std::ptrdiff_t difference_type;
    typedef V3 *pointer;
    typedef V3 &reference;
    typedef std::input_iterator_tag iterator_category;

    inline iterator() : v(nullptr), i(0) {}
    inline iterator(V3v *vv, unsigned int ii) : v(vv), i(ii) {}

    inline V3 &operator*() { return (*v)[i]; }
    inline const V3 &operator*() const { return (*v)[i]; }
    inline double &operator[](int j) { return (*v)[i][j]; }
    inline const double &operator[](int j) const { return (*v)[i][j]; }

    inline iterator &operator++() {
      ++i;
      return *this;
    }
    inline iterator &operator--() {
      --i;
      return *this;
    }
    inline iterator operator++(int) {
      iterator r(*this);
      ++i;
      return r;
    }
    inline iterator operator--(int) {
      iterator r(*this);
      --i;
      return r;
    }

    inline iterator &operator+=(int nn) {
      i += nn;
      return *this;
    }
    inline iterator &operator-=(int nn) {
      i -= nn;
      return *this;
    }

    inline iterator operator+(int nn) const {
      iterator r(*this);
      return r += nn;
    }
    inline iterator operator-(int nn) const {
      iterator r(*this);
      return r -= nn;
    }

    inline difference_type operator-(iterator const &r) const {
      return i - r.i;
    }

    inline bool operator<(iterator const &r) const { return i < r.i; }
    inline bool operator<=(iterator const &r) const { return i <= r.i; }
    inline bool operator>(iterator const &r) const { return i > r.i; }
    inline bool operator>=(iterator const &r) const { return i >= r.i; }
    inline bool operator!=(const iterator &r) const { return i != r.i; }
    inline bool operator==(const iterator &r) const { return i == r.i; }

    inline unsigned int ind() const { return i; }

   private:
    V3v *v;
    unsigned int i;
  };

  inline iterator begin() { return V3v::iterator{this, 0}; }
  inline iterator end() { return V3v::iterator{this, n}; }

 private:
  unsigned int n;
  V3 *v;
  bool alloc;
};

// Constructor and destructor definitions
inline V3v::V3v(const unsigned int nn) {
  alloc = false;
  n = 0;
  size(nn);
}
inline V3v::V3v(const unsigned int nn, const float *val) {
  alloc = false;
  n = 0;
  size(nn, val);
}
inline V3v::V3v(const unsigned int nn, const double *val) {
  alloc = false;
  n = 0;
  size(nn, val);
}
inline V3v::V3v(const unsigned int nn, const V3 *val) {
  alloc = false;
  n = 0;
  size(nn, val);
}
inline V3v::V3v(const V3v &val) {
  alloc = false;
  n = 0;
  size(val);
}
inline V3v::~V3v() { clear(); }

// Operators
inline V3v &V3v::operator=(const V3v &other) {
  size(other.n, other.v);
  return (*this);
}

// Functions
inline void V3v::size(const unsigned int nn) {
  n = nn;
  v = new V3[n];
  alloc = true;
}
inline void V3v::size(const unsigned int nn, const V3 *val) {
  size(nn);
  std::memcpy(v, val, n * sizeof(V3));
}
inline void V3v::size(const V3v &val) { size(val.n, val.v); }
inline void V3v::size(const unsigned int nn, const float *val) {
  size(nn);
#pragma omp parallel for
  for (unsigned int i = 0; i < n; i++)
    v[i] = V3((double)(val[3 * i + 0]), (double)(val[3 * i + 1]),
              (double)(val[3 * i + 2]));
}
inline void V3v::size(const unsigned int nn, const double *val) {
  size(nn);
#pragma omp parallel for
  for (unsigned int i = 0; i < n; i++)
    v[i] = V3(val[3 * i + 0], val[3 * i + 1], val[3 * i + 2]);
}

inline void V3v::clear() {
  n = 0;
  if (alloc) {
    delete[] v;
  }
  alloc = false;
}
inline V3 *V3v::data() { return v; }
inline bool V3v::isempty() { return n == 0; }

template <typename T>
T *V3v::toarray() const {
  auto *out = new T[3 * n];
  for (unsigned int i = 0; i < n; i++) {
    out[3 * i + 0] = static_cast<T>(v[i][0]);
    out[3 * i + 1] = static_cast<T>(v[i][1]);
    out[3 * i + 2] = static_cast<T>(v[i][2]);
  }
  return out;
}

}  // namespace OGS::V3

#endif
