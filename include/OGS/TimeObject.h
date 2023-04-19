/*=========================================================================

  Module:    Time Object

  Class to manage struct tm.

  Copyright (c) 2019 Arnau Miro, OGS
  All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef TimeObject_h
#define TimeObject_h

#include <algorithm>
#include <cstring>
#include <ctime>
#include <stdexcept>
#include <string>
#include <vector>

namespace OGS::Time {

class TimeObject {
 public:
  inline TimeObject() : t(time(nullptr)), ltm(*localtime(&t)) {
    ltm.tm_isdst = 0;
  }
  inline explicit TimeObject(const tm &in) : t(time(nullptr)), ltm(in) {
    update();
  }
  inline TimeObject(const TimeObject &TO) : TimeObject(TO.ltm) {}
  inline TimeObject(const char *timestr, const char *fmt)
      : t(time(nullptr)), ltm(*localtime(&t)) {
    from_string(timestr, fmt);
  }
  inline TimeObject(const std::string &timestr, const char *fmt)
      : TimeObject(timestr.c_str(), fmt) {}
  inline TimeObject(const std::string &timestr, const std::string &fmt)
      : TimeObject(timestr, fmt.c_str()) {}

  [[nodiscard]] inline tm get_tm() const { return ltm; }

  inline void set_tm(const tm &in) {
    ltm = in;
    update();
  }

  [[nodiscard]] inline time_t epoch() const { return t; }

  [[nodiscard]] inline int isoweekday() const {
    return std::stoi(as_string("%u"));
  }

  inline void from_string(const char *str, const char *fmt) {
    strptime(str, fmt, &ltm);
    update();
  }

  inline void from_string(const std::string &str, const std::string &fmt) {
    from_string(str.c_str(), fmt.c_str());
  }

  inline std::string as_string(const char *fmt) const {
    char buff[256];
    std::strftime(buff, 256, fmt, &ltm);
    return {buff};
  }

  inline std::string as_string(std::string &fmt) const {
    return as_string(fmt.c_str());
  }

  inline void increment_sec(const int i) {
    ltm.tm_sec += i;
    update();
  }
  inline void increment_min(const int i) {
    ltm.tm_min += i;
    update();
  }
  inline void increment_hour(const int i) {
    ltm.tm_hour += i;
    update();
  }
  inline void increment_day(const int i) {
    ltm.tm_mday += i;
    update();
  }
  inline void increment_month(const int i) {
    ltm.tm_mon += i;
    update();
  }
  inline void increment_year(const int i) {
    ltm.tm_year += i;
    update();
  }
  inline void decrement_sec(const int i) {
    ltm.tm_sec -= i;
    update();
  }
  inline void decrement_min(const int i) {
    ltm.tm_min -= i;
    update();
  }
  inline void decrement_hour(const int i) {
    ltm.tm_hour -= i;
    update();
  }
  inline void decrement_day(const int i) {
    ltm.tm_mday -= i;
    update();
  }
  inline void decrement_month(const int i) {
    ltm.tm_mon -= i;
    update();
  }
  inline void decrement_year(const int i) {
    ltm.tm_year -= i;
    update();
  }

  inline TimeObject &operator=(const tm &in) {
    set_tm(in);
    return (*this);
  }
  inline TimeObject &operator=(const TimeObject &TO) {
    set_tm(TO.ltm);
    return (*this);
  }
  inline time_t operator+(const TimeObject &TO) const { return t + TO.t; }
  inline time_t operator-(const TimeObject &TO) const { return t - TO.t; }

  inline bool operator==(const TimeObject &TO) const { return t == TO.t; }
  inline bool operator<(const TimeObject &TO) const { return t < TO.t; }
  inline bool operator<=(const TimeObject &TO) const { return t <= TO.t; }
  inline bool operator>(const TimeObject &TO) const { return t > TO.t; }
  inline bool operator>=(const TimeObject &TO) const { return t >= TO.t; }

 private:
  time_t t;
  std::tm ltm;

  inline void update() {
    t = timegm(&ltm);
    ltm = *gmtime(&t);
    ltm.tm_isdst = 0;
  }
};

class TimeObjectList {
 public:
  TimeObjectList() : n(0), alloc(false), TOList(nullptr) {}
  explicit TimeObjectList(const unsigned int nn)
      : n(nn), alloc(true), TOList(new TimeObject[n]) {}
  TimeObjectList(const unsigned int nn, TimeObject *TO) : TimeObjectList(nn) {
    fill(TO);
  }
  TimeObjectList(const TimeObjectList &TOLi)
      : TimeObjectList(TOLi.n, TOLi.TOList) {}
  ~TimeObjectList() { clear(); }

  TimeObjectList(const TimeObject &start, const TimeObject &end,
                 const char *delta);

  [[nodiscard]] inline unsigned int len() const { return n; }

  inline void clear() {
    if (this->alloc) {
      delete[] this->TOList;
      this->alloc = false;
      this->n = 0;
    }
  }

  inline TimeObject *data() { return TOList; }

  [[nodiscard]] inline bool isempty() const { return n == 0; }

  std::string as_string(const char *fmt) {
    std::string str = std::string("[") + TOList[0].as_string(fmt);
    for (unsigned int ii = 1; ii < n; ++ii)
      str += std::string(", ") + TOList[ii].as_string(fmt);
    str += std::string("]");
    return str;
  }

  void sort() {
    if (this->alloc) std::sort(begin(), end());
  }

  int find(const TimeObject &TO);

  TimeObjectList merge(TimeObjectList &TOLi, bool exact);

  inline TimeObject operator[](int i) const {
    return (i >= 0) ? TOList[i] : TOList[n + i];
  }
  inline TimeObject &operator[](int i) {
    return (i >= 0) ? TOList[i] : TOList[n + i];
  }

  inline TimeObject operator[](unsigned int i) const { return TOList[i]; }
  inline TimeObject &operator[](unsigned int i) { return TOList[i]; }

  inline TimeObject operator[](size_t i) const { return TOList[i]; }
  inline TimeObject &operator[](size_t i) { return TOList[i]; }

  TimeObjectList &operator=(const TimeObjectList &TOLi) {
    if (!TOLi.alloc) {
      if (this->alloc) this->clear();
      return (*this);
    }

    if (this->n != TOLi.n) {
      if (this->alloc) clear();
      this->n = TOLi.n;
      this->allocate();
    }

    if (this->TOList != TOLi.TOList) this->fill(TOLi.TOList);

    return (*this);
  }

  class iterator {
   public:
    typedef TimeObject value_type;
    typedef std::ptrdiff_t difference_type;
    typedef TimeObject *pointer;
    typedef TimeObject &reference;
    typedef std::input_iterator_tag iterator_category;

    iterator(TimeObjectList *TOLi, unsigned int ii) : TOL(TOLi), i(ii) {}
    iterator(const iterator &other) = default;
    iterator &operator=(const iterator &TOLi) = default;

    inline TimeObject &operator*() { return (*TOL)[i]; }
    inline const TimeObject &operator*() const { return (*TOL)[i]; }

    inline iterator &operator++() {
      ++i;
      return *this;
    }
    inline iterator &operator--() {
      --i;
      return *this;
    }

    inline iterator operator++(int) & {
      iterator r(*this);
      ++i;
      return r;
    }
    inline iterator operator--(int) & {
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

    [[nodiscard]] inline unsigned int ind() const { return i; }

   private:
    TimeObjectList *TOL;
    unsigned int i;
  };

  inline iterator begin() { return iterator{this, 0}; }
  inline iterator end() { return iterator{this, n}; }

 private:
  unsigned int n;
  bool alloc;
  TimeObject *TOList = nullptr;

  inline void allocate() {
    if (!this->alloc) {
      TOList = new TimeObject[n];
      alloc = true;
    } else
      throw std::runtime_error(
          "Trying to allocate an already allocated TimeObjectList");
  }

  inline void fill(TimeObject *TO) {
    if (alloc)
      std::copy(TO, TO + this->n, this->TOList);
    else
      throw std::runtime_error("Filling a non initialized TimeObjectList");
  }
};


TimeObjectList::TimeObjectList(const TimeObject &start, const TimeObject &end,
                               const char *delta) {
  /*
          Given a start and an end point, construct a time object list spaced
     every delta. Options for delta are: > "secs   = x" : every x seconds; >
     "mins   = x" : every x minutes; > "hours  = x" : every x hours; > "days   =
     x" : every x days; > "months = x" : every x months; > "years  = x" : every
     x years;
  */
  // First parse the string
  std::string str(delta);
  size_t pos = str.find('=');
  std::string delta_type = str.substr(0, pos);
  int d = std::stoi(str.substr(pos + 1, str.size()));

  TimeObject TO_next(start), TO_end(end);
  std::vector<TimeObject> TO_list;
  this->n = 0;
  this->alloc = false;

  do {
    TO_list.emplace_back(TO_next);
    this->n++;
    if (delta_type == std::string("secs"))
      TO_next.increment_sec(d);
    else if (delta_type == std::string("mins"))
      TO_next.increment_min(d);
    else if (delta_type == std::string("hours"))
      TO_next.increment_hour(d);
    else if (delta_type == std::string("days"))
      TO_next.increment_day(d);
    else if (delta_type == std::string("months"))
      TO_next.increment_month(d);
    else if (delta_type == std::string("years"))
      TO_next.increment_year(d);
  } while (TO_next <= TO_end);

  allocate();
  fill(&TO_list[0]);
}


int TimeObjectList::find(const TimeObject &TO) {
  /*
          Finds the nearest
          Argument:
                  a datetime object
          Returns the index of the nearest element
  */
  int posmin = -1;
  time_t diffmin = 9999999;
  for (int ii = 0; ii < n; ++ii) {
    time_t diff = TOList[ii] - TO;
    if (diff < diffmin) {
      diffmin = diff;
      posmin = ii;
    }
  }
  return posmin;
}


TimeObjectList TimeObjectList::merge(TimeObjectList &TOLi, bool exact = false) {
  /*
          Merge with a datetimelist
          Returns a datetimelist with all the dates
          from both lists without repetitions
  */
  std::vector<TimeObject> list1(this->n + TOLi.len());
  std::vector<TimeObject> list2(this->n + TOLi.len());

  // from 0 to n copy this->TOList
  std::copy(this->TOList, this->TOList + n, list1.data());
  // from n to TOLi.len() copy TOLi
  std::copy(TOLi.data(), TOLi.data() + TOLi.len(), list1.data() + n);

  std::sort(list1.begin(), list1.end());  // sort the list

  // Now remove duplicates from list1 to list2
  std::vector<TimeObject>::iterator iter;
  if (exact)
    iter = std::unique_copy(list1.begin(), list1.end(), list2.begin());
  else
    iter = std::unique_copy(list1.begin(), list1.end(), list2.begin(),
                            [](const TimeObject &T1, const TimeObject &T2) {
                              return labs(T1 - T2) < 24 * 3600;
                            });

  // We remove the last elements of list2, i.e., the ones that where allocated
  // but we've never used (because the unique_copy stopped before reaching them)
  list2.resize(std::distance(list2.begin(), iter));

  // List2 contains the unique values of the two lists
  return {static_cast<unsigned int>(list2.size()), list2.data()};
}
}  // namespace OGS::Time

#endif
