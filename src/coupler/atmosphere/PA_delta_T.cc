// Copyright (C) 2011, 2012, 2013, 2014, 2015 PISM Authors
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

#include "PA_delta_T.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/io/io_helpers.hh"

namespace pism {
namespace atmosphere {

/// delta_T forcing of near-surface air temperatures

Delta_T::Delta_T(IceGrid::ConstPtr g, AtmosphereModel* in)
  : PScalarForcing<AtmosphereModel,PAModifier>(g, in),
    air_temp(m_sys, "air_temp"),
    precipitation(m_sys, "precipitation")
{
  offset = NULL;
  option_prefix = "-atmosphere_delta_T";
  offset_name   = "delta_T";

  offset = new Timeseries(*m_grid, offset_name, m_config->get_string("time_dimension_name"));
  offset->metadata().set_string("units", "Kelvin");
  offset->metadata().set_string("long_name", "near-surface air temperature offsets");
  offset->dimension_metadata().set_string("units", m_grid->ctx()->time()->units_string());
  
  air_temp.set_string("pism_intent", "diagnostic");
  air_temp.set_string("long_name", "near-surface air temperature");
  air_temp.set_string("units", "K");

  precipitation.set_string("pism_intent", "diagnostic");
  precipitation.set_string("long_name", "precipitation, units of ice-equivalent thickness per time");
  precipitation.set_string("units", "m / s");
  precipitation.set_string("glaciological_units", "m / year");
}

void Delta_T::init() {

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  input_model->init();

  m_log->message(2,
             "* Initializing near-surface air temperature forcing using scalar offsets...\n");

  init_internal();
}

MaxTimestep Delta_T::max_timestep_impl(double t) {
  (void) t;
  return MaxTimestep();
}

void Delta_T::init_timeseries(const std::vector<double> &ts) {

  PAModifier::init_timeseries(ts);

  m_offset_values.resize(m_ts_times.size());
  for (unsigned int k = 0; k < m_ts_times.size(); ++k) {
    m_offset_values[k] = (*offset)(m_ts_times[k]);
  }
}

void Delta_T::mean_annual_temp(IceModelVec2S &result) {

  input_model->mean_annual_temp(result);
  offset_data(result);
}

void Delta_T::temp_time_series(int i, int j, std::vector<double> &result) {
  input_model->temp_time_series(i, j, result);

  for (unsigned int k = 0; k < m_ts_times.size(); ++k) {
    result[k] += m_offset_values[k];
  }
}

void Delta_T::temp_snapshot(IceModelVec2S &result) {
  input_model->temp_snapshot(result);
  offset_data(result);
}

void Delta_T::add_vars_to_output_impl(const std::string &keyword, std::set<std::string> &result) {
  input_model->add_vars_to_output(keyword, result);

  if (keyword == "medium" || keyword == "big") {
    result.insert("air_temp");
    result.insert("precipitation");
  }
}


void Delta_T::define_variables_impl(const std::set<std::string> &vars_input, const PIO &nc,
                                            IO_Type nctype) {
  std::set<std::string> vars = vars_input;
  std::string order = m_grid->ctx()->config()->get_string("output_variable_order");

  if (set_contains(vars, "air_temp")) {
    io::define_spatial_variable(air_temp, *m_grid, nc, nctype, order, false);
    vars.erase("air_temp");
  }

  if (set_contains(vars, "precipitation")) {
    io::define_spatial_variable(precipitation, *m_grid, nc, nctype, order, true);
    vars.erase("precipitation");
  }

  input_model->define_variables(vars, nc, nctype);
}


void Delta_T::write_variables_impl(const std::set<std::string> &vars_input, const PIO &nc) {
  std::set<std::string> vars = vars_input;

  if (set_contains(vars, "air_temp")) {
    IceModelVec2S tmp;
    tmp.create(m_grid, "air_temp", WITHOUT_GHOSTS);
    tmp.metadata() = air_temp;

    mean_annual_temp(tmp);

    tmp.write(nc);

    vars.erase("air_temp");
  }

  if (set_contains(vars, "precipitation")) {
    IceModelVec2S tmp;
    tmp.create(m_grid, "precipitation", WITHOUT_GHOSTS);
    tmp.metadata() = precipitation;

    mean_precipitation(tmp);

    tmp.write_in_glaciological_units = true;
    tmp.write(nc);

    vars.erase("precipitation");
  }

  input_model->write_variables(vars, nc);
}

} // end of namespace atmosphere
} // end of namespace pism
