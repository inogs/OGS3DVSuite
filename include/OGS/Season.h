/*=========================================================================

  Module:    Season

  Implementation of bit.sea Season in C++.

  Copyright (c) 2019 Arnau Miro, OGS
  All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef season_h
#define season_h

#include <string>
#include <vector>

#include "TimeInterval.h"
#include "TimeObject.h"

namespace OGS::Time {

static const std::vector<std::string> DEF_START_SEASON{
    std::string("0101"), std::string("0401"), std::string("0701"),
    std::string("1001")};

static const std::vector<std::string> DEF_NAME_SEASON{
    std::string("winter"), std::string("spring"), std::string("summer"),
    std::string("fall")};


class Season {
 public:
  explicit Season(const int reference_year)
      : numbers_season(0), reference_year(reference_year), alloc(false) {}

  Season() : Season(2000) {}

  Season(int reference_year, const std::vector<std::string> &start_season,
         const std::vector<std::string> &name_season);

  Season(const std::vector<std::string> &start_season,
         const std::vector<std::string> &name_season)
      : Season(2000, start_season, name_season) {}


  ~Season() { this->clear(); }

  [[nodiscard]] inline int get_reference_year() const {
    return this->reference_year;
  };

  inline void set_reference_year(int reference_year_) {
    this->reference_year = reference_year_;
  };

  inline void set_seasons(const std::vector<std::string> &start_season,
                          const std::vector<std::string> &name_season);

  [[nodiscard]] inline unsigned int get_seasons_number() const {
    return this->numbers_season;
  }

  inline TimeInterval get_season_dates(int season_num,
                                       std::string &season_name);

  inline int find_season(const TimeObject &TO);

 private:
  unsigned int numbers_season;
  int reference_year;
  bool alloc;
  std::vector<TimeObject> SEASON_LIST;
  std::vector<std::string> SEASON_LIST_NAME;

  inline void allocate();
  inline void clear();
};



inline void Season::allocate() {
  if (this->alloc)
    throw std::runtime_error("Allocating an already allocated Season object");

  this->SEASON_LIST.clear();
  this->SEASON_LIST_NAME.clear();

  this->SEASON_LIST.resize(this->numbers_season);
  this->SEASON_LIST_NAME.resize(this->numbers_season);
  this->alloc = true;
}



inline void Season::clear() {
  this->SEASON_LIST.clear();
  this->SEASON_LIST_NAME.clear();
  this->alloc = false;
}



inline Season::Season(const int reference_year,
                      const std::vector<std::string> &start_season,
                      const std::vector<std::string> &name_season) {
  this->numbers_season = 0;
  this->reference_year = reference_year;
  this->alloc = false;
  this->set_seasons(start_season, name_season);
}


inline void Season::set_seasons(const std::vector<std::string> &start_season,
                                const std::vector<std::string> &name_season) {
  /*
Given two arrays where is defined the date when Season starts and its name
the subroutine generate the season list. In input take two arrays of string:

- Example of start_season : ["1221","0321","0622","0921"]
- Example of name_season  : ["winter","spring","summer","fall"]

In the previous example we shown that winter start to 21 december,
spring to 21 march, summer to 22 june and fall on 21 september.
  */
  this->numbers_season = start_season.size();
  this->allocate();

  if (start_season.size() != name_season.size()) {
    fprintf(stderr, "ERROR : arrays definitions mismatch!\n");
    exit(-1);
  }

  int ref_year = (start_season[0] == std::string("0101"))
                     ? this->reference_year
                     : this->reference_year - 1;
  std::string time_str =
      std::to_string(ref_year) + start_season[0] + std::string("-00:00:00");
  this->SEASON_LIST.emplace(this->SEASON_LIST.begin(), time_str.c_str(),
                            "%Y%m%d-%H:%M:%S");
  this->SEASON_LIST_NAME.insert(this->SEASON_LIST_NAME.begin(), name_season[0]);

  for (int ii = 1; ii < this->numbers_season; ++ii) {
    ref_year = this->reference_year;
    time_str =
        std::to_string(ref_year) + start_season[ii] + std::string("-00:00:00");
    this->SEASON_LIST.emplace(this->SEASON_LIST.begin() + ii, time_str.c_str(),
                              "%Y%m%d-%H:%M:%S");
    this->SEASON_LIST_NAME.insert(this->SEASON_LIST_NAME.begin() + ii,
                                  name_season[ii]);
  }
}



inline TimeInterval Season::get_season_dates(const int season_num,
                                             std::string &season_name) {
  /*
  Given Season number, return the range of Season dates (start and end)
  and the name of Season.
  */
  assert(season_num < this->numbers_season);

  TimeObject start_date(this->SEASON_LIST[season_num]);
  TimeObject end_date;

  if (season_num + 1 == this->numbers_season) {
    end_date = this->SEASON_LIST[0];
    end_date.increment_year(1);
  } else {
    end_date = this->SEASON_LIST[season_num + 1];
  }

  season_name = this->SEASON_LIST_NAME[season_num];

  // TimeInterval TI(start_date,end_date);
  return {start_date, end_date};
}

inline int Season::find_season(const TimeObject &TO) {
  /*
  Takes a date as input and return the number and name of Season where it is in.
  */
  TimeInterval TI;
  std::string season_name;
  TimeObject aux(TO);

  int delta_year = this->reference_year - std::stoi(aux.as_string("%Y"));

  aux.increment_year(delta_year);
  for (int ii = 0; ii < this->numbers_season; ++ii) {
    TI = this->get_season_dates(ii, season_name);
    if (TI.contains(aux)) return ii;
  }

  aux.decrement_year(1);
  for (int ii = 0; ii < this->numbers_season; ++ii) {
    TI = this->get_season_dates(ii, season_name);
    if (TI.contains(aux)) return ii;
  }
  return -1;
}
}  // namespace OGS::Time

#endif
