/*=========================================================================

  Module:    Time Interval

  Implementation of bit.sea TimeInterval in C++.

  Copyright (c) 2019 Arnau Miro, OGS
  All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef TimeInterval_h
#define TimeInterval_h

#include <cassert>
#include <string>

#include "TimeObject.h"

namespace OGS::Time {

class TimeInterval {
 public:
  TimeInterval() = default;
  TimeInterval(const TimeInterval &TI) = default;
  TimeInterval(const TimeObject &start_time_, const TimeObject &end_time_)
      : start_time(start_time_), end_time(end_time_) {
    assert(this->start_time <= this->end_time);
  }
  TimeInterval(const char *starttime, const char *endtime,
               const char *dateformat);
  TimeInterval(std::string &starttime, std::string &endtime,
               std::string &dateformat);

  [[nodiscard]] inline TimeObject get_start_time() const {
    return this->start_time;
  }
  [[nodiscard]] inline TimeObject get_end_time() const {
    return this->end_time;
  }
  inline void set_start_time(const TimeObject &start_time_) {
    this->start_time = start_time_;
  }
  inline void set_end_time(const TimeObject &end_time_) {
    this->end_time = end_time_;
  }

  std::string as_string(const char *fmt) const;

  inline bool contains(const TimeObject &specific_time);
  inline time_t overlapTime(const TimeInterval &T2);
  inline bool isInWindow(const TimeInterval &T2);
  inline bool isInside(const TimeInterval &T2);
  [[nodiscard]] inline time_t length() const;

  inline TimeInterval &operator=(const TimeInterval &TI) = default;

 private:
  TimeObject start_time, end_time;
};

inline TimeInterval::TimeInterval(const char *starttime, const char *endtime,
                                  const char *dateformat) {
  this->start_time = TimeObject(starttime, dateformat);
  this->end_time = TimeObject(endtime, dateformat);
  assert(this->start_time <= this->end_time);
}
inline TimeInterval::TimeInterval(std::string &starttime, std::string &endtime,
                                  std::string &dateformat) {
  this->start_time = TimeObject(starttime, dateformat);
  this->end_time = TimeObject(endtime, dateformat);
  assert(this->start_time <= this->end_time);
}

inline std::string TimeInterval::as_string(const char *fmt) const {
  return this->start_time.as_string(fmt) + std::string(" , ") +
         this->end_time.as_string(fmt);
}

inline bool TimeInterval::contains(const TimeObject &specific_time) {
  /*
  Argument:
      * specific_time * : a datetime object
  Returns:
      True if specific_time is inside the time interval
  */
  return (specific_time >= this->start_time && specific_time <= this->end_time);
}



inline time_t TimeInterval::overlapTime(const TimeInterval &T2) {
  /*
  Argument:
          T2 : a TimeInterval object
  Returns:
          the number of seconds of the overlapping interval
*---s1++++++++++e1------------ * Window of temporal average of model output file
*               |
*-----------s2++++++++++e2---- * Window of aggregation defined in PLOT_LIST
*            |  |
*-----------<++++>------------ * overlapping time window = theWindow
  */
  time_t theWindow = std::min(this->end_time.epoch(), T2.end_time.epoch()) -
                     std::max(this->start_time.epoch(), T2.start_time.epoch());
  return std::max((time_t)(0), theWindow);
}



inline bool TimeInterval::isInWindow(const TimeInterval &T2) {
  /*
  T2 : a TimeInterval object
  Returns: True if there is an overlapping interval, False elsewhere.
  */
  return this->overlapTime(T2) > 0;
}



inline bool TimeInterval::isInside(const TimeInterval &T2) {
  (void)T2;
  fprintf(stderr, "Not Implemented!\n");
  exit(-1);
}



inline time_t TimeInterval::length() const {
  /*
  Returns lengths of time interval in seconds
  */
  return this->end_time - this->start_time;
}
}  // namespace OGS::Time

#endif
