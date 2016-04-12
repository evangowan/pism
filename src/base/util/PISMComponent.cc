// Copyright (C) 2008-2016 Ed Bueler and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include <gsl/gsl_math.h>
#include <cassert>

#include "PISMComponent.hh"
#include "base/util/io/PIO.hh"
#include "IceGrid.hh"
#include "pism_const.hh"
#include "pism_utilities.hh"
#include "VariableMetadata.hh"
#include "iceModelVec.hh"
#include "pism_options.hh"
#include "error_handling.hh"
#include "PISMConfigInterface.hh"
#include "MaxTimestep.hh"

namespace pism {

Component::Component(IceGrid::ConstPtr g)
  : m_grid(g), m_config(g->ctx()->config()), m_sys(g->ctx()->unit_system()),
    m_log(g->ctx()->log()) {
  // empty
}

Component::~Component() {
  // empty
}

void Component::define_variables(const std::set<std::string> &vars, const PIO &nc,
                                 IO_Type nctype) {
  this->define_variables_impl(vars, nc, nctype);
}

void Component::add_vars_to_output(const std::string &keyword,
                                   std::set<std::string> &result) {
  this->add_vars_to_output_impl(keyword, result);
}

void Component::write_variables(const std::set<std::string> &vars, const PIO& nc) {
  this->write_variables_impl(vars, nc);
}

void Component::get_diagnostics(std::map<std::string, Diagnostic::Ptr> &dict,
                                std::map<std::string, TSDiagnostic::Ptr> &ts_dict) {
  this->get_diagnostics_impl(dict, ts_dict);
}

void Component::get_diagnostics_impl(std::map<std::string, Diagnostic::Ptr> &dict,
                                     std::map<std::string, TSDiagnostic::Ptr> &ts_dict) {
  (void)dict;
  (void)ts_dict;
}

IceGrid::ConstPtr Component::grid() const {
  return m_grid;
}

//! Finds PISM's input file by reading options `-i` and `-bootstrap`.
/*! This might be useful since coupling fields are usually in the file
  IceModel uses to initialize from.

  Returns `true` if an input file is used and `false` otherwise.
*/
bool Component::find_pism_input(std::string &filename, bool &do_regrid, int &start) {

  // read file name options:
  options::String input_file("-i", "input file name");
  bool bootstrap = options::Bool("-bootstrap", "enable bootstrapping heuristics");

  if (input_file.is_set()) {

    filename = input_file;

    PIO nc(m_grid->com, "netcdf3");      // OK to use netcdf3
    unsigned int last_record;
    nc.open(filename, PISM_READONLY);
    last_record = nc.inq_nrecords() - 1;
    nc.close();

    if (bootstrap) {
      do_regrid = true;
      start     = 0;
    } else {
      do_regrid = false;
      start     = last_record;
    }

    return true;
  } else {
    filename.clear();
    do_regrid = false;
    start     = -1;
    return false;
  }
}

/**
 * Regrid a variable by processing -regrid_file and -regrid_vars.
 *
 * @param[in] module_name Module name, used to annotate options when run with -help.
 *
 * @param[out] variable pointer to an IceModelVec; @c variable has to
 *             have metadata set for this to work.
 *
 * @param[in] flag Regridding flag. If set to
 *            REGRID_WITHOUT_REGRID_VARS, regrid this variable by
 *            default, if =-regrid_vars= was not set. Otherwise a
 *            variable is only regridded if both =-regrid_file= and
 *            =-regrid_vars= are set *and* the name of the variable is
 *            found in the set of names given with =-regrid_vars=.
 */
void Component::regrid(const std::string &module_name, IceModelVec &variable,
                       RegriddingFlag flag) {

  options::String regrid_file("-regrid_file", "regridding file name");

  options::StringSet regrid_vars("-regrid_vars",
                                 "comma-separated list of regridding variables",
                                 "");

  if (not regrid_file.is_set()) {
    return;
  }

  SpatialVariableMetadata &m = variable.metadata();

  if ((regrid_vars.is_set() and set_contains(regrid_vars, m.get_string("short_name"))) or
      (not regrid_vars.is_set() and flag == REGRID_WITHOUT_REGRID_VARS)) {

    m_log->message(2,
               "  %s: regridding '%s' from file '%s' ...\n",
               module_name.c_str(),
               m.get_string("short_name").c_str(), regrid_file->c_str());

    variable.regrid(regrid_file, CRITICAL);
  }
}

Component_TS::Component_TS(IceGrid::ConstPtr g)
  : Component(g) {
  m_t = m_dt = GSL_NAN;
}

Component_TS::~Component_TS() {
  // empty
}

MaxTimestep Component_TS::max_timestep(double t) {
  return this->max_timestep_impl(t);
}

void Component_TS::update(double t, double dt) {
  this->update_impl(t, dt);
}

} // end of namespace pism
