/* Copyright (C) 2013, 2014, 2015 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "PISMUnits.hh"

#include "base/util/error_handling.hh"

#include "pism_const.hh"

namespace pism {

namespace units {

class ut_system_deleter {
public:
  void operator()(ut_system* p) const {
    ut_free_system(p);
  }
};

/** Initialize the unit system by reading from an XML unit
 * definition file.
 */
System::System(const std::string &path) {
  ut_system *tmp;

  ut_set_error_message_handler(ut_ignore);

  if (path.empty() == false) {
    tmp = ut_read_xml(path.c_str());
  } else {
    tmp = ut_read_xml(NULL);
  }

  if (tmp == NULL) {
    std::string message = std::string("ut_read_xml(") + path + ") failed";
    throw RuntimeError(message);
  }
  ut_set_error_message_handler(ut_write_to_stderr);

  m_system = PISM_SHARED_PTR(ut_system)(tmp, ut_system_deleter());
}

//! \brief Convert a quantity from unit1 to unit2.
/*!
 * Example: convert(1, "m/year", "m/s").
 *
 * Please avoid using in computationally-intensive code.
 */
double convert(System::Ptr system, double input,
               const std::string &spec1, const std::string &spec2) {
  Converter c(Unit(system, spec1), Unit(system, spec2));

  return c(input);
}

Unit::Unit(System::Ptr system, const std::string &spec)
  : m_unit(NULL), m_system(system) {
  m_unit = ut_parse(m_system->m_system.get(), spec.c_str(), UT_ASCII);
  if (m_unit == NULL) {
    std::string message = "unit specification '" + spec + "' is unknown or invalid";
    throw RuntimeError(message);
  }
  m_unit_string = spec;
}

Unit::Unit(const Unit &other)
  : m_system(other.m_system) {

  m_unit = ut_clone(other.m_unit);
  if (m_unit == NULL) {
    throw RuntimeError("ut_clone failed");
  }

  m_system      = other.m_system;
  m_unit_string = other.m_unit_string;
}

Unit& Unit::operator=(const Unit& other) {
  if (this == &other) {
    return *this;
  }

  reset();

  m_system      = other.m_system;
  m_unit_string = other.m_unit_string;

  m_unit = ut_clone(other.m_unit);
  if (m_unit == NULL) {
    throw RuntimeError("ut_clone failed");
  }

  return *this;
}

Unit::~Unit() {
  reset();
}

std::string Unit::format() const {
  return m_unit_string;
}

void Unit::reset() {
  ut_free(m_unit);
  m_unit = NULL;
}

ut_unit* Unit::get() const {
  return m_unit;
}

System::Ptr Unit::get_system() const {
  return m_system;
}

Converter::Converter() {
  m_converter = cv_get_trivial();
}

Converter::Converter(System::Ptr sys,
                     const std::string &spec1, const std::string &spec2) {

  Unit u1(sys, spec1), u2(sys, spec2);

  if (ut_are_convertible(u1.get(), u2.get()) == 0) {
    throw RuntimeError::formatted("cannot convert %s to %s", spec1.c_str(), spec2.c_str());
  }

  m_converter = ut_get_converter(u1.get(), u2.get());
  if (m_converter == NULL) {
    throw RuntimeError::formatted("cannot create a converter from %s to %s",
                                  spec1.c_str(), spec2.c_str());
  }
}

Converter::Converter(const Unit &u1, const Unit &u2) {
  if (ut_are_convertible(u1.get(), u2.get()) == 0) {
    std::string message = "cannot convert " + u1.format() + " to " + u2.format();
    throw RuntimeError(message);
  }

  m_converter = ut_get_converter(u1.get(), u2.get());
  if (m_converter == NULL) {
    std::string message = "cannot create a converter from " + u1.format() + " to " + u2.format();
    throw RuntimeError(message);
  }
}

bool are_convertible(const Unit &u1, const Unit &u2) {
  return ut_are_convertible(u1.get(), u2.get()) != 0;
}

Converter::~Converter() {
  cv_free(m_converter);
  m_converter = NULL;
}

double Converter::operator()(double input) const {
  return cv_convert_double(m_converter, input);
}

void Converter::convert_doubles(double *data, size_t length) const {
  cv_convert_doubles(m_converter, data, length, data);
}

} // end of namespace units

} // end of namespace pism
