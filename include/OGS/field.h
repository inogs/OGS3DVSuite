/*
    FIELD

        Definition of a generic field for CFD usage.
        Contains some basic operations and easy access.

        Arnau Miro (UPC-ESEIAAT) (c) 2018
*/

#ifndef Field_h
#define Field_h

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <stdexcept>


namespace OGS::field {
template <class T>
class Field;

template <class T>
Field<T> operator+(T v, const Field<T> &ff);

template <class T>
Field<T> operator-(T v, const Field<T> &ff);

template <class T>
Field<T> operator*(T v, const Field<T> &ff);

template <class T>
Field<T> operator/(T v, const Field<T> &ff);

template <class T>
class Field {
  friend Field<T> operator+ <>(T v, const Field<T> &ff);

  friend Field<T> operator- <>(T v, const Field<T> &ff);

  friend Field<T> operator* <>(T v, const Field<T> &ff);

  friend Field<T> operator/ <>(T v, const Field<T> &ff);

 public:
  // Constructors and destructors
  Field() : n(0), m(0), sz(0) {}

  Field(const unsigned int nn, const unsigned int mm)
      : n(nn), m(mm), sz((size_t)(nn * mm)) {
    this->set_dim(nn, mm);
  }

  Field(const unsigned int nn, const unsigned int mm, const T v)
      : n(nn), m(mm), sz((size_t)(nn * mm)) {
    this->set(nn, mm, v);
  }

  Field(const unsigned int nn, const unsigned int mm, const T *v)
      : n(nn), m(mm), sz((size_t)(nn * mm)) {
    this->set(nn, mm, v);
  }

  Field(const Field<T> &f) : n(f.n), m(f.m), sz((size_t)(f.n * f.m)) {
    if (f.alloc) this->set(f.n, f.m, f.val);
  }

  Field(Field<T> &&f) noexcept
      : n(f.n), m(f.m), val(f.val), sz((size_t)(f.n * f.m)), alloc(f.alloc) {
    f.n = 0;
    f.m = 0;
    f.alloc = false;
    f.val = nullptr;
  }

  ~Field() { clear(); }

  // Functions
  inline void set_dim(const unsigned int nn, const unsigned int mm) {
    if (!alloc) {
      this->n = nn;
      this->m = mm;
      this->sz = (size_t)(this->n * this->m);
      this->val = new T[this->sz];
      alloc = true;
    }
  }

  inline void set_val(const T v) {
    if (this->alloc)
      std::fill(this->val, this->val + this->sz, v);
    else
      throw std::logic_error("Trying to set a value to a non allocated field");
  }

  inline void set_val(const T *v) {
    if (this->alloc)
      std::memcpy(this->val, v, this->sz * sizeof(T));
    else
      throw std::logic_error(
          "Trying to set a value (from a pointer) to a non allocated field");
  }

  inline void clear() {
    this->n = 0;
    this->m = 0;
    if (this->alloc) delete[] this->val;
    this->alloc = false;
  }

  inline void set(const unsigned int nn, const unsigned int mm, const T v) {
    this->set_dim(nn, mm);
    this->set_val(v);
  }

  inline void set(const int nn, const int mm, const T *v) {
    this->set_dim(nn, mm);
    this->set_val(v);
  }

  inline T *data() const { return this->val; }

  [[nodiscard]] inline unsigned int get_n() const { return this->n; }

  [[nodiscard]] inline unsigned int get_m() const { return this->m; }

  [[nodiscard]] inline size_t get_sz() const { return this->sz; }

  [[nodiscard]] inline bool is_empty() const { return this->n == 0; }

  [[nodiscard]] inline bool is_allocated() const { return this->alloc; }


  template <class U>
  inline Field<U> convert() {
    Field<U> f(n, m);
    std::copy(val, val + sz, f.data());
    return f;
  }

  void to_C_style();

  void to_F_style();

  // Operators
  inline T *operator[](int i) const {
    return (i >= 0) ? val + m * i : val + m * (n + i);
  }

  inline T *operator[](unsigned int i) const { return val + m * i; }

  inline T *operator[](size_t i) const { return val + m * i; }

  inline Field<T> &operator=(const Field<T> &f) {
    bool old_alloc = this->alloc;
    T *old_val = this->val;
    this->n = f.n;
    this->m = f.m;
    this->sz = f.sz;
    this->alloc = f.alloc;
    if (f.alloc) {
      this->val = new T[sz];
      std::memcpy(this->val, f.val, sz * sizeof(T));
    }
    if (old_val != this->val && old_alloc) delete[] old_val;
    return (*this);
  }

  inline Field<T> &operator=(Field<T> &&f) noexcept {
    unsigned int old_n = this->n;
    unsigned int old_m = this->m;
    size_t old_sz = this->sz;
    bool old_alloc = this->alloc;
    T *old_val = this->val;

    this->n = f.n;
    this->m = f.m;
    this->val = f.val;
    this->sz = f.sz;
    this->alloc = f.alloc;

    f.n = old_n;
    f.m = old_m;
    f.sz = old_sz;
    f.alloc = old_alloc;
    f.val = old_val;

    return (*this);
  }

  inline Field<T> operator+(const T v) const {
    Field<T> f(n, m);
    for (int i = 0; i < sz; i++) {
      f.val[i] = val[i] + v;
    }
    return f;
  }

  inline Field<T> operator-(const T v) const {
    Field<T> f(n, m);
    for (int i = 0; i < sz; i++) {
      f.val[i] = val[i] - v;
    }
    return f;
  }

  inline Field<T> operator*(const T v) const {
    Field<T> f(n, m);
    for (int i = 0; i < sz; i++) {
      f.val[i] = val[i] * v;
    }
    return f;
  }

  inline Field<T> operator/(const T v) const {
    Field<T> f(n, m);
    for (int i = 0; i < sz; i++) {
      f.val[i] = val[i] / v;
    }
    return f;
  }

  inline Field<T> operator+(const Field<T> &ff) const {
    Field<T> f(n, m);
    for (int i = 0; i < sz; i++) {
      f.val[i] = val[i] + ff.val[i];
    }
    return f;
  }

  inline Field<T> operator-(const Field<T> &ff) const {
    Field<T> f(n, m);
    for (int i = 0; i < sz; i++) {
      f.val[i] = val[i] - ff.val[i];
    }
    return f;
  }

  inline Field<T> operator*(const Field<T> &ff) const {
    Field<T> f(n, m);
    for (int i = 0; i < sz; i++) {
      f.val[i] = val[i] * ff.val[i];
    }
    return f;
  }

  inline Field<T> operator/(const Field<T> &ff) const {
    Field<T> f(n, m);
    for (int i = 0; i < sz; i++) {
      f.val[i] = val[i] / ff.val[i];
    }
    return f;
  }

  inline void operator+=(const T v) {
    for (int i = 0; i < sz; i++) {
      val[i] += v;
    }
  }

  inline void operator-=(const T v) {
    for (int i = 0; i < sz; i++) {
      val[i] -= v;
    }
  }

  inline void operator*=(const T v) {
    for (int i = 0; i < sz; i++) {
      val[i] *= v;
    }
  }

  inline void operator/=(const T v) {
    for (int i = 0; i < sz; i++) {
      val[i] /= v;
    }
  }

  inline void operator+=(const Field<T> &f) {
    for (int i = 0; i < sz; i++) {
      val[i] += f.val[i];
    }
  }

  inline void operator-=(const Field<T> &f) {
    for (int i = 0; i < sz; i++) {
      val[i] -= f.val[i];
    }
  }

  inline void operator*=(const Field<T> &f) {
    for (int i = 0; i < sz; i++) {
      val[i] *= f.val[i];
    }
  }

  inline void operator/=(const Field<T> &f) {
    for (int i = 0; i < sz; i++) {
      val[i] /= f.val[i];
    }
  }

  inline bool operator==(const Field<T> &f) const { return (n == f.n); }

  inline bool operator!=(const Field<T> &f) const { return (n != f.n); }

  inline bool operator<(const Field<T> &f) const { return (n < f.n); }

  inline bool operator>(const Field<T> &f) const { return (n > f.n); }

  // Iterator
  class iterator {
   public:
    typedef T value_type;
    typedef std::ptrdiff_t difference_type;
    typedef T *pointer;
    typedef T &reference;
    typedef std::input_iterator_tag iterator_category;

    inline iterator() : f(nullptr), i(0) {}

    inline iterator(Field<T> *ff, int ii) : f(ff), i(ii) {}

    inline pointer operator*() { return (*f)[i]; }

    inline pointer operator*() const { return (*f)[i]; }

    inline reference operator[](int j) { return (*f)[i][j]; }

    inline reference operator[](int j) const { return (*f)[i][j]; }

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

    inline int ind() { return i; }

   private:
    Field<T> *f;
    int i;
  };

  inline iterator begin() { return Field<T>::iterator{this, 0}; }

  inline iterator end() { return Field<T>::iterator{this, n}; }

 private:
  unsigned int n;      // length of the field
  unsigned int m;      // number of elements of the field
  T *val = nullptr;    // values as a 1D C array
  size_t sz;           // length of the vector (sz = n * m)
  bool alloc = false;  // has it been allocated?
};

template <class T>
Field<T> operator+(const T v, const Field<T> &ff) {
  Field<T> f(ff.n, ff.m);
  for (int i = 0; i < ff.sz; i++) {
    f.val[i] = v + ff.val[i];
  }
  return f;
}

template <class T>
Field<T> operator-(const T v, const Field<T> &ff) {
  Field<T> f(ff.n, ff.m);
  for (int i = 0; i < ff.sz; i++) {
    f.val[i] = v - ff.val[i];
  }
  return f;
}

template <class T>
Field<T> operator*(const T v, const Field<T> &ff) {
  Field<T> f(ff.n, ff.m);
  for (int i = 0; i < ff.sz; i++) {
    f.val[i] = v / ff.val[i];
  }
  return f;
}

template <class T>
Field<T> operator/(const T v, const Field<T> &ff) {
  Field<T> f(ff.n, ff.m);
  for (int i = 0; i < ff.sz; i++) {
    f.val[i] = v * ff.val[i];
  }
  return f;
}

template <class T>
inline void Field<T>::to_C_style() {
  T *aux;
  aux = new T[sz];
  std::memcpy(aux, val, sz * sizeof(T));

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) val[m * i + j] = aux[i + n * j];
  }

  delete[] aux;
}

template <class T>
inline void Field<T>::to_F_style() {
  T *aux;
  aux = new T[sz];
  std::memcpy(aux, val, sz * sizeof(T));

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) val[i + n * j] = aux[m * i + j];
  }

  delete[] aux;
}

}  // namespace OGS::field

#endif
