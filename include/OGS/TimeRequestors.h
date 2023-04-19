/*=========================================================================

  Module:    Time Requestors

  Implementation of bit.sea TimeRequestors in C++.

  Copyright (c) 2019 Arnau Miro, OGS
  All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef TimeRequestors_h
#define TimeRequestors_h

#include <cassert>
#include <string>

#include "Season.h"
#include "TimeInterval.h"
#include "TimeObject.h"

namespace OGS::Time {
class Requestor {
  /*
  Generic requestor (base class)
  Based on timeinterval object, for non standard times.
  For example, an year defined between
  */
 public:
  inline Requestor() = default;
  inline explicit Requestor(const TimeInterval &TI_,
                            const char *fmt = "%Y%m%d-%H:%M:%S")
      : TI(TI_), request_str(TI_.as_string(fmt)) {}
  [[nodiscard]] virtual inline std::string as_string() const {
    return (std::string("Req: ") + this->request_str);
  }
  [[nodiscard]] virtual inline std::string type() const { return "generic"; }
  virtual inline bool contains(const TimeObject &TO) {
    return (this->TI.contains(TO));
  }
  [[nodiscard]] virtual inline TimeInterval interval() const {
    return (this->TI);
  }

 protected:
  TimeInterval TI;
  std::string request_str;
};

class Decadal_req : public Requestor {
  /*
  Decadal requestor
  Decadal_req(2010) requests years 2010,2011,..., 2019
  Decadal_req(2011) requests years 2011,2012,..., 2020
  */
 public:
  inline explicit Decadal_req(const int year) : Requestor() {
    // Decades are defined to start at 0 up to 9 (e.g.,2010-2019)
    // or start at 1 up to 10 (e.g., 2011-2020)
    // anything else is considered a generic interval.
    assert(year % 10 == 0 || year % 10 == 1);

    char buff[256];
    std::sprintf(buff, "%d0101-00:00:00", year);

    TimeObject start_time(buff, "%Y%m%d-%H:%M:%S"),
        end_time(buff, "%Y%m%d-%H:%M:%S");
    end_time.increment_year(10);  // maybe 9 here...

    this->TI = TimeInterval(start_time, end_time);
    this->request_str = start_time.as_string("%Y") + std::string(" ... ") +
                        end_time.as_string("%Y");
  }

  [[nodiscard]] inline std::string as_string() const override {
    return (std::string("Decadal_req: ") + this->request_str);
  }
  [[nodiscard]] inline std::string type() const override { return "decadal"; }
};

class Yearly_req : public Requestor {
  /*
  Requestor object - for specific year - used in Timelist.select() method.
  */
 public:
  inline explicit Yearly_req(const int year) : Requestor() {
    char buff[256];
    std::sprintf(buff, "%d0101-00:00:00", year);

    TimeObject start_time(buff, "%Y%m%d-%H:%M:%S"),
        end_time(buff, "%Y%m%d-%H:%M:%S");
    end_time.increment_year(1);

    this->TI = TimeInterval(start_time, end_time);
    this->request_str = start_time.as_string("%Y");
  }

  [[nodiscard]] inline std::string as_string() const override {
    return (std::string("Yearly_req: ") + this->request_str);
  }
  [[nodiscard]] inline std::string type() const override { return "yearly"; }
};

class Season_req : public Requestor {
  /*
  Requestor object - for a specific Season - used in Timelist.select() method.
  */
 public:
  inline Season_req(const int year, const int iseason, Season &s_obj)
      : Requestor() {
    /*
    Arguments:
            * year      * : integer
            * Season    * : integer
            * seasonobj * : an instance of class Season, where features of
    Season are defined

    Then seasonObj defines by default seasons like this:
            Winter = 0
            Spring = 1
            Summer = 2
            Fall   = 3
    */
    std::string longname;
    TimeInterval TI_ref = s_obj.get_season_dates(iseason, longname);

    int d_years = year - s_obj.get_reference_year();
    TimeObject start_time = TI_ref.get_start_time();
    start_time.increment_year(d_years);
    TimeObject end_time = TI_ref.get_end_time();
    end_time.increment_year(d_years);

    this->TI = TimeInterval(start_time, end_time);
    this->request_str =
        start_time.as_string("%Y") + std::string(" ") + longname.substr(0, 3);
  }

  [[nodiscard]] inline std::string as_string() const override {
    return (std::string("Season_req: ") + this->request_str);
  }
  [[nodiscard]] inline std::string type() const override { return "Season"; }
};

class Monthly_req : public Requestor {
 public:
  inline Monthly_req(const int year, const int month) : Requestor() {
    char buff[256];

    sprintf(buff, "%d%02d01-00:00:00", year, month);

    TimeObject start_time(buff, "%Y%m%d-%H:%M:%S"),
        end_time(buff, "%Y%m%d-%H:%M:%S");
    end_time.increment_month(1);

    this->TI = TimeInterval(start_time, end_time);
    this->request_str = start_time.as_string("%Y%m");
  }

  [[nodiscard]] inline std::string as_string() const override {
    return (std::string("Monthly_req: ") + this->request_str);
  }
  [[nodiscard]] inline std::string type() const override { return "monthly"; }
};

class Weekly_req : public Requestor {
  /*
  Requestor object - for a specific week - used in Timelist.select() method.
  The week is indicated by its central day.

  Example:
          req=Weekly_req(2015,3,5)
          print r.isoweekday
  */
 public:
  inline Weekly_req(const int year, const int month, const int day)
      : Requestor() {
    char buff[256];

    sprintf(buff, "%d%02d%02d-12:00:00", year, month, day);

    TimeObject centertime(buff, "%Y%m%d-%H:%M:%S");

    int delta = (int)(3.5 * 24. * 3600.);  // half a week, 3.5 days to seconds
    auto start_time = TimeObject(centertime);
    start_time.increment_sec(-delta);
    auto end_time = TimeObject(centertime);
    end_time.increment_sec(delta - 1);

    this->TI = TimeInterval(start_time, end_time);
    this->request_str = centertime.as_string("%Y%m%d");
  }

  [[nodiscard]] inline std::string as_string() const override {
    return (std::string("Weekly_req: ") + this->request_str);
  }
  [[nodiscard]] inline std::string type() const override { return "weekly"; }
};

class Interval_req : public Requestor {
  /*
  Useful for 10 days averages
  Interval_req(2005,1,25,days=10)
  */
 public:
  inline Interval_req(const int year, const int month, const int day,
                      const int hour, const int days)
      : Requestor() {
    char buff[256];

    sprintf(buff, "%d%02d%02d-%02d:00:00", year, month, day, hour);

    TimeObject centertime(buff, "%Y%m%d-%H:%M:%S");

    int delta = static_cast<int>(days * 24. * 3600. / 2.);
    auto start_time = TimeObject(centertime);
    start_time.increment_sec(-delta);
    auto end_time = TimeObject(centertime);
    end_time.increment_sec(delta);

    this->TI = TimeInterval(start_time, end_time);
    this->request_str = centertime.as_string("%Y%m%d") +
                        std::string(" Delta: ") + std::to_string(days) +
                        std::string(" days");
  }

  [[nodiscard]] inline std::string as_string() const override {
    return (std::string("Interval_req: ") + this->request_str);
  }
  [[nodiscard]] inline std::string type() const override { return "interval"; }
};

class Daily_req : public Requestor {
  /*
  Requestor object - for a specific day - used in Timelist.select() method.

  Example:
          req=Daily_req(2015,3,5)
  */
 public:
  inline Daily_req(const int year, const int month, const int day)
      : Requestor() {
    char buff[256];

    sprintf(buff, "%d%02d%02d-00:00:00", year, month, day);

    TimeObject start_time(buff, "%Y%m%d-%H:%M:%S"),
        end_time(buff, "%Y%m%d-%H:%M:%S");
    end_time.increment_day(1);

    this->TI = TimeInterval(start_time, end_time);
    this->request_str = start_time.as_string("%Y%m%d");
  }

  [[nodiscard]] inline std::string as_string() const override {
    return (std::string("Daily_req: ") + this->request_str);
  }
  [[nodiscard]] inline std::string type() const override { return "daily"; }
};

class Hourly_req : public Requestor {
  /*
  Requestor object - for a specific hour interval
  */
 public:
  inline Hourly_req(const int year, const int month, const int day,
                    const int hour, const int hours)
      : Requestor() {
    char buff[256];

    sprintf(buff, "%d%02d%02d-%02d:00:00", year, month, day, hour);

    TimeObject centertime(buff, "%Y%m%d-%H:%M:%S");

    int delta = (int)(hours * 3600. / 2.);
    auto start_time = TimeObject(centertime);
    start_time.increment_sec(-delta);
    auto end_time = TimeObject(centertime);
    end_time.increment_sec(delta);

    this->TI = TimeInterval(start_time, end_time);
    this->request_str = centertime.as_string("%Y%m%d-%H:%M:%S") +
                        std::string(" Delta: ") + std::to_string(hours) +
                        std::string(" hours");
  }

  [[nodiscard]] inline std::string as_string() const override {
    return (std::string("Hourly_req: ") + this->request_str);
  }
  [[nodiscard]] inline std::string type() const override { return "hourly"; }
};

class Seconds_req : public Requestor {
  /*
  Requestor object - for a specific seconds interval
  */
 public:
  inline Seconds_req(const int year, const int month, const int day,
                     const int hour, const int minute, const int delta_seconds)
      : Requestor() {
    char buff[256];

    sprintf(buff, "%d%02d%02d-%02d:%02d:00", year, month, day, hour, minute);

    TimeObject centertime(buff, "%Y%m%d-%H:%M:%S");

    int delta = delta_seconds / 2;
    auto start_time = TimeObject(centertime);
    start_time.increment_sec(-delta);
    auto end_time = TimeObject(centertime);
    end_time.increment_sec(delta);

    this->TI = TimeInterval(start_time, end_time);
    this->request_str = centertime.as_string("%Y%m%d-%H:%M:%S") +
                        std::string(" Delta: ") +
                        std::to_string(delta_seconds) + std::string(" seconds");
  }

  [[nodiscard]] inline std::string as_string() const override {
    return (std::string("Seconds_req: ") + this->request_str);
  }
  [[nodiscard]] inline std::string type() const override { return "seconds"; }
};

class Clim_season : public Requestor {
  /*
  Requestor object  - for a climatologit Season - used in Timelist.select()
  method.
  */
 public:
  inline Clim_season(const int iseason, Season &s_obj) : Requestor() {
    /*
    Arguments:
            * season    * : integer
            * seasonobj * : an instance of class Season, where features of
    Season are defined

    Then seasonObj defines by default seasons like this:
            Winter = 0
            Spring = 1
            Summer = 2
            Fall   = 3
    */
    assert(iseason < s_obj.get_seasons_number());

    this->ref_year = s_obj.get_reference_year();
    this->TI = s_obj.get_season_dates(iseason, this->request_str);
  }

  [[nodiscard]] inline std::string as_string() const override {
    return (std::string("Clim_season: ") + this->request_str);
  }
  [[nodiscard]] inline std::string type() const override {
    return "clim_season";
  }

  inline bool contains(const TimeObject &TO) override {
    /*
    Argument:
            * time * : a datetime object
    Returns:
            True if time is inside every time interval of the Season
    */
    auto aux = TimeObject(TO);
    int d_years = this->ref_year - std::stoi(aux.as_string("%Y"));
    aux.increment_year(d_years);
    return this->TI.contains(aux);
  }

 private:
  int ref_year;
};

class Clim_month : public Requestor {
  /*
  Requestor object - for a climatologic month - used in Timelist.select()
  method.
  */
 public:
  inline explicit Clim_month(const int month_) : Requestor(), month(month_) {
    assert(month_ < 13 && month_ > 0);

    char buff[256];

    sprintf(buff, "2000%02d01-00:00:00", month_);

    TimeObject start_time(buff, "%Y%m%d-%H:%M:%S"),
        end_time(buff, "%Y%m%d-%H:%M:%S");
    end_time.increment_month(1);

    this->TI = TimeInterval(start_time, end_time);
    this->request_str = start_time.as_string("%b");
  }

  [[nodiscard]] inline std::string as_string() const override {
    return (std::string("Clim_month: ") + this->request_str);
  }
  [[nodiscard]] inline std::string type() const override {
    return "clim_month";
  }

  inline bool contains(const TimeObject &TO) override {
    /*
    Argument:
            * time * : a datetime object
    Returns:
            True if time is inside every time interval of the month
    */
    auto aux = TimeObject(TO);
    int d_years = 2000 - std::stoi(aux.as_string("%Y"));
    aux.increment_year(d_years);
    return this->TI.contains(aux);
  }

 private:
  const int month;
};

class Clim_day : public Requestor {
  /*
  Requestor object for - climatologic day -used in Timelist.select() method.
  */
 public:
  inline Clim_day(int month_, int day_) : Requestor(), day(day_) {
    assert(month_ < 13 && month_ > 0);
    assert(day_ < 31 && day_ > 0);

    char buff[256];
    sprintf(buff, "2001%02d%02d-00:00:00", month_, day_);

    TimeObject start_time(buff, "%Y%m%d-%H:%M%S"),
        end_time(buff, "%Y%m%d-%H:%M%S");
    end_time.increment_day(1);

    this->TI = TimeInterval(start_time, end_time);
    this->request_str = std::to_string(day_);
  }

  [[nodiscard]] inline std::string as_string() const override {
    return (std::string("Clim_day: ") + this->request_str);
  }
  [[nodiscard]] inline std::string type() const override { return "clim_day"; }

  inline bool contains(const TimeObject &TO) override {
    /*
    Argument:
            * time * : a datetime object
    Returns:
            True if time is inside every time interval of the day
    */
    return Daily_req(std::stoi(TO.as_string("%Y")),
                     std::stoi(TO.as_string("%m")), this->day)
        .contains(TO);
  }

 private:
  const int day;
};

class Clim_hour : public Requestor {
  /*
   */
 public:
  inline Clim_hour(const int hour, const int hours) : Requestor() {
    this->hour = hour;
    this->hours = hours;

    this->request_str =
        std::to_string(hour) + std::string(" Delta: ") + std::to_string(hours);
  }

  [[nodiscard]] inline std::string as_string() const override {
    return (std::string("Clim_hour: ") + this->request_str);
  }
  [[nodiscard]] inline std::string type() const override { return "clim_hour"; }

  inline bool contains(const TimeObject &TO) override {
    /*
    Argument:
            * time * : a datetime object
    Returns:
            True if time is inside every time interval of the month
    */
    return Hourly_req(std::stoi(TO.as_string("%Y")),
                      std::stoi(TO.as_string("%m")),
                      std::stoi(TO.as_string("%d")), this->hour, this->hours)
        .contains(TO);
  }

 private:
  int hour, hours;
};
}  // namespace OGS::Time

#endif
