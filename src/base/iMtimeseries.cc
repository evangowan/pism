// Copyright (C) 2009-2015 Constantine Khroulev
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

#include <sstream>
#include <algorithm>
#include <gsl/gsl_interp.h>

#include "iceModel.hh"

#include "base/stressbalance/PISMStressBalance.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/PISMDiagnostic.hh"
#include "base/util/PISMTime.hh"
#include "base/util/error_handling.hh"
#include "base/util/io/PIO.hh"
#include "base/util/pism_options.hh"
#include "base/util/PISMVars.hh"
#include "base/util/io/io_helpers.hh"
#include "base/util/MaxTimestep.hh"
#include "base/util/Profiling.hh"

namespace pism {

//! Initializes the code writing scalar time-series.
void IceModel::init_timeseries() {

  options::String ts_file("-ts_file", "Specifies the time-series output file name");
  ts_filename = ts_file;

  options::String times("-ts_times", "Specifies a MATLAB-style range or a list of requested times");

  options::StringSet vars("-ts_vars", "Specifies a comma-separated list of veriables to save",
                          "");

  // default behavior is to move the file aside if it exists already; option allows appending
  bool append = options::Bool("-ts_append", "append scalar time-series");


  IO_Mode mode = PISM_READWRITE;
  if (not append) {
    mode = PISM_READWRITE_MOVE;
  }

  if (ts_file.is_set() ^ times.is_set()) {
    throw RuntimeError("you need to specity both -ts_file and -ts_times to save diagnostic time-series.");
  }

  // If neither -ts_file nor -ts_times is set, we're done.
  if (not ts_file.is_set() && not times.is_set()) {
    save_ts = false;
    return;
  }

  save_ts = true;

  try {
    m_time->parse_times(times, ts_times);
  } catch (RuntimeError &e) {
    e.add_context("parsing the -ts_times argument %s", times->c_str());
    throw;
  }

  if (times->empty()) {
    throw RuntimeError("no argument for -ts_times option.");
  }

  m_log->message(2, "saving scalar time-series to '%s'; ",
             ts_file->c_str());

  m_log->message(2, "times requested: %s\n", times->c_str());

  current_ts = 0;

  if (vars.is_set()) {
    m_log->message(2, "variables requested: %s\n", vars.to_string().c_str());
    ts_vars = vars;
  } else {
    std::map<std::string,TSDiagnostic*>::iterator j = ts_diagnostics.begin();
    while (j != ts_diagnostics.end()) {
      ts_vars.insert(j->first);
      ++j;
    }
  }

  PIO nc(m_grid->com, "netcdf3");      // Use NetCDF-3 to write time-series.
  nc.open(ts_file, mode);

  if (append == true) {
    double time_max;
    std::string time_name = m_config->get_string("time_dimension_name");
    bool time_exists = false;

    time_exists = nc.inq_var(time_name);
    if (time_exists == true) {
      nc.inq_dim_limits(time_name, NULL, &time_max);

      while (current_ts < ts_times.size() && ts_times[current_ts] < time_max) {
        current_ts++;
      }

      if (current_ts > 0) {
        m_log->message(2,
                   "skipping times before the last record in %s (at %s)\n",
                   ts_file->c_str(), m_time->date(time_max).c_str());
      }
    }
  }

  write_metadata(nc, false, false);

  nc.close();


  // set the output file:
  std::map<std::string,TSDiagnostic*>::iterator j = ts_diagnostics.begin();
  while (j != ts_diagnostics.end()) {
    (j->second)->init(ts_file);
    ++j;
  }

  // ignore times before (and including) the beginning of the run:
  while (current_ts < ts_times.size() && ts_times[current_ts] < m_time->start()) {
    current_ts++;
  }

  if (ts_times.size() == current_ts) {
    save_ts = false;
    return;
  }

  // discard requested times before the beginning of the run
  std::vector<double> tmp(ts_times.size() - current_ts);
  for (unsigned int k = 0; k < tmp.size(); ++k) {
    tmp[k] = ts_times[current_ts + k];
  }

  ts_times = tmp;
  current_ts = 0;
}

//! Write time-series.
void IceModel::write_timeseries() {

  // return if no time-series requested
  if (!save_ts) {
     return;
  }

  // return if wrote all the records already
  if (current_ts == ts_times.size()) {
    return;
  }

  // return if did not yet reach the time we need to save at
  if (ts_times[current_ts] > m_time->current()) {
    return;
  }

  for (std::set<std::string>::iterator j = ts_vars.begin(); j != ts_vars.end(); ++j) {
    TSDiagnostic *diag = ts_diagnostics[*j];

    if (diag != NULL) {
      diag->update(m_time->current() - dt, m_time->current());
    }
  }


  // Interpolate to put them on requested times:
  for (; current_ts < ts_times.size() && ts_times[current_ts] <= m_time->current(); current_ts++) {

    // the very first time (current_ts == 0) defines the left endpoint of the
    // first time interval; we don't write a report at that time
    if (current_ts == 0) {
      continue;
    }

    for (std::set<std::string>::iterator j = ts_vars.begin(); j != ts_vars.end(); ++j) {
      TSDiagnostic *diag = ts_diagnostics[*j];

      if (diag != NULL) {
        diag->save(ts_times[current_ts - 1], ts_times[current_ts]);
      }
    }
  }
}


//! Initialize the code saving spatially-variable diagnostic quantities.
void IceModel::init_extras() {

  last_extra = 0;               // will be set in write_extras()
  next_extra = 0;

  options::String extra_file("-extra_file", "Specifies the output file");
  extra_filename = extra_file;

  options::String times("-extra_times", "Specifies times to save at");

  options::StringSet vars("-extra_vars",
                          "Specifies a comma-separated list of variables to save", "");

  bool split  = options::Bool("-extra_split", "Specifies whether to save to separate files");
  bool append = options::Bool("-extra_append", "append spatial diagnostics");

  if (extra_file.is_set() ^ times.is_set()) {
    throw RuntimeError("you need to specify both -extra_file and -extra_times to save spatial time-series.");
  }

  if (!extra_file.is_set() && !times.is_set()) {
    save_extra = false;
    return;
  }

  try {
    m_time->parse_times(times, extra_times);
  } catch (RuntimeError &e) {
    e.add_context("parsing the -extra_times argument %s", times->c_str());
    throw;
  }

  if (extra_times.size() == 0) {
    throw RuntimeError("no argument for -extra_times option.");
  }

  if (append && split) {
    throw RuntimeError("both -extra_split and -extra_append are set.");
  }

  if (append) {
    PIO nc(m_grid->com, m_config->get_string("output_format"));
    std::string time_name = m_config->get_string("time_dimension_name");
    bool time_exists;

    nc.open(extra_filename, PISM_READONLY);
    time_exists = nc.inq_var(time_name);

    if (time_exists == true) {
      double time_max;
      nc.inq_dim_limits(time_name, NULL, &time_max);

      while (next_extra + 1 < extra_times.size() && extra_times[next_extra + 1] < time_max) {
        next_extra++;
      }

      if (next_extra > 0) {
        m_log->message(2,
                   "skipping times before the last record in %s (at %s)\n",
                   extra_filename.c_str(), m_time->date(time_max).c_str());
      }

      // discard requested times before the beginning of the run
      std::vector<double> tmp(extra_times.size() - next_extra);
      for (unsigned int k = 0; k < tmp.size(); ++k) {
        tmp[k] = extra_times[next_extra + k];
      }

      extra_times = tmp;
      next_extra = 0;
    }
    nc.close();
  }

  save_extra          = true;
  extra_file_is_ready = false;
  split_extra         = false;

  if (split) {
    split_extra = true;
    m_log->message(2, "saving spatial time-series to '%s+year.nc'; ",
               extra_filename.c_str());
  } else {
    if (!ends_with(extra_filename, ".nc")) {
      m_log->message(2,
                 "PISM WARNING: spatial time-series file name '%s' does not have the '.nc' suffix!\n",
                 extra_filename.c_str());
    }
    m_log->message(2, "saving spatial time-series to '%s'; ",
               extra_filename.c_str());
  }

  m_log->message(2, "times requested: %s\n", times->c_str());

  if (extra_times.size() > 500) {
    m_log->message(2,
               "PISM WARNING: more than 500 times requested. This might fill your hard-drive!\n");
  }

  if (vars.is_set()) {
    m_log->message(2, "variables requested: %s\n", vars.to_string().c_str());
    extra_vars = vars;
  } else {
    m_log->message(2, "PISM WARNING: -extra_vars was not set."
               " Writing model_state, mapping and climate_steady variables...\n");

    std::set<std::string> vars_set = m_grid->variables().keys();

    std::set<std::string>::iterator i;
    for (i = vars_set.begin(); i != vars_set.end(); ++i) {
      const SpatialVariableMetadata &m = m_grid->variables().get(*i)->metadata();

      std::string intent = m.get_string("pism_intent");

      if (intent == "model_state" ||
          intent == "mapping"     ||
          intent == "climate_steady") {
        extra_vars.insert(*i);
      }
    }

    std::set<std::string> list;
    if (stress_balance) {
      stress_balance->add_vars_to_output("small", extra_vars);
    }

  } // end of the else clause after "if (extra_vars_set)"

  if (extra_vars.size() == 0) {
    m_log->message(2,
               "PISM WARNING: no variables list after -extra_vars ... writing empty file ...\n");
  }
}

//! Write spatially-variable diagnostic quantities.
void IceModel::write_extras() {
  double saving_after = -1.0e30; // initialize to avoid compiler warning; this
                                 // value is never used, because saving_after
                                 // is only used if save_now == true, and in
                                 // this case saving_after is guaranteed to be
                                 // initialized. See the code below.
  char filename[PETSC_MAX_PATH_LEN];
  unsigned int current_extra;
  // determine if the user set the -save_at and -save_to options
  if (!save_extra) {
    return;
  }

  // do we need to save *now*?
  if (next_extra < extra_times.size() &&
      (m_time->current() >= extra_times[next_extra] ||
       fabs(m_time->current() - extra_times[next_extra]) < 1.0)) {
    // the condition above is "true" if we passed a requested time or got to
    // within 1 second from it

    current_extra = next_extra;

    // update next_extra
    while (next_extra < extra_times.size() &&
           (extra_times[next_extra] <= m_time->current() ||
            fabs(m_time->current() - extra_times[next_extra]) < 1.0)) {
      next_extra++;
    }

    saving_after = extra_times[current_extra];
  } else {
    return;
  }

  if (current_extra == 0) {
    // The first time defines the left end-point of the first reporting
    // interval; we don't write a report at this time, but we still need to
    // store cumulative quantities that may be needed to compute rates of
    // change.

    std::set<std::string>::iterator j = extra_vars.begin();
    while(j != extra_vars.end()) {
      Diagnostic *diag = diagnostics[*j];

      if (diag != NULL) {
        diag->update_cumulative();
      }
      ++j;
    }

    // This line re-initializes last_extra (the correct value is not known at
    // the time init_extras() is calles).
    last_extra = m_time->current();

    return;
  }

  if (saving_after < m_time->start()) {
    // Suppose a user tells PISM to write data at times 0:1000:10000. Suppose
    // also that PISM writes a backup file at year 2500 and gets stopped.
    //
    // When restarted, PISM will decide that it's time to write data for time
    // 2000, but
    // * that record was written already and
    // * PISM will end up writing at year 2500, producing a file containing one
    //   more record than necessary.
    //
    // This check makes sure that this never happens.
    return;
  }

  const Profiling &profiling = m_ctx->profiling();

  profiling.begin("extra_file reporting");

  if (split_extra) {
    extra_file_is_ready = false;        // each time-series record is written to a separate file
    snprintf(filename, PETSC_MAX_PATH_LEN, "%s-%s.nc",
             extra_filename.c_str(), m_time->date().c_str());
  } else {
    strncpy(filename, extra_filename.c_str(), PETSC_MAX_PATH_LEN);
  }

  m_log->message(3,
             "\nsaving spatial time-series to %s at %s\n\n",
             filename, m_time->date().c_str());

  // find out how much time passed since the beginning of the run
  double wall_clock_hours = pism::wall_clock_hours(m_grid->com, start_time);

  PIO nc(m_grid->com, m_config->get_string("output_format"));

  if (extra_file_is_ready == false) {
    // default behavior is to move the file aside if it exists already; option allows appending
    bool append = options::Bool("-extra_append", "append -extra_file output");

    IO_Mode mode = PISM_READWRITE;
    if (not append) {
      mode = PISM_READWRITE_MOVE;
    }

    // Prepare the file:
    nc.open(filename, mode);
    io::define_time(nc, m_config->get_string("time_dimension_name"),
                    m_time->calendar(),
                    m_time->CF_units_string(),
                    m_sys);
    nc.put_att_text(m_config->get_string("time_dimension_name"),
                    "bounds", "time_bounds");

    write_metadata(nc, true, false);

    extra_file_is_ready = true;
  } else {
    // In this case the extra file should be present.
    nc.open(filename, PISM_READWRITE);
  }

  double      current_time = m_time->current();
  std::string time_name    = m_config->get_string("time_dimension_name");


  unsigned int time_length = nc.inq_dimlen(time_name);
  size_t time_start = static_cast<size_t>(time_length);

  // This call will extend the time dimension, but that will not
  // happen until nc.enddef() is called. (We don't want to switch to
  // "data mode" before we're done defining all variables, including
  // time bounds). This is why time_start = time_length above (and not
  // time_length - 1).
  io::append_time(nc, time_name, current_time);

  std::vector<double> data(2);
  data[0] = last_extra;
  data[1] = current_time;
  io::write_time_bounds(nc, extra_bounds, time_start, data);

  io::write_timeseries(nc, timestamp, time_start, wall_clock_hours);

  write_variables(nc, extra_vars, PISM_FLOAT);

  nc.close();

  // flush time-series buffers
  flush_timeseries();

  last_extra = current_time;

  profiling.end("extra_file reporting");
}

static MaxTimestep reporting_max_timestep(const std::vector<double> &times, double t) {

  const size_t N = times.size();
  if (t >= times.back()) {
    return MaxTimestep();
  }

  size_t j = 0;
  double dt = 0.0;
  if (t < times[0]) {
    j = -1;
  } else {
    j = gsl_interp_bsearch(&times[0], t, 0, N - 1);
  }

  dt = times[j + 1] - t;

  // now make sure that we don't end up taking a time-step of less than 1
  // second long
  if (dt < 1.0) {
    if (j + 2 < N) {
      return MaxTimestep(times[j + 2] - t);
    } else {
      return MaxTimestep();
    }
  } else {
    return MaxTimestep(dt);
  }
}

//! Computes the maximum time-step we can take and still hit all `-extra_times`.
MaxTimestep IceModel::extras_max_timestep(double my_t) {

  if ((not save_extra) or
      (not m_config->get_boolean("extras_force_output_times"))) {
    return MaxTimestep();
  }

  return reporting_max_timestep(extra_times, my_t);
}

//! Computes the maximum time-step we can take and still hit all `-ts_times`.
MaxTimestep IceModel::ts_max_timestep(double my_t) {

  if ((not save_ts) or
      (not m_config->get_boolean("ts_force_output_times"))) {
    return MaxTimestep();
  }

  return reporting_max_timestep(ts_times, my_t);
}

//! Flush scalar time-series.
void IceModel::flush_timeseries() {
  // flush all the time-series buffers:
  for (std::set<std::string>::iterator j = ts_vars.begin(); j != ts_vars.end(); ++j) {
    TSDiagnostic *diag = ts_diagnostics[*j];

    if (diag != NULL) {
      diag->flush();
    }
  }
}

} // end of namespace pism
