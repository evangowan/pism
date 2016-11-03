/* Copyright (C) 2015, 2016 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
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

#include <vector>
#include <cassert>

#include <gsl/gsl_interp.h>

#include "base/util/pism_options.hh"
#include "base/util/error_handling.hh"
#include "base/util/IceGrid.hh"
#include "base/util/io/PIO.hh"

namespace pism {

static void validate_range(const std::string &axis,
                           const std::vector<double> &x,
                           double x_min, double x_max) {
  if (x_min >= x_max) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "invalid -%s_range: %s_min >= %s_max (%f >= %f).",
                                  axis.c_str(), axis.c_str(), axis.c_str(),
                                  x_min, x_max);
  }

  if (x_min >= x.back()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "invalid -%s_range: %s_min >= max(%s) (%f >= %f).",
                                  axis.c_str(), axis.c_str(), axis.c_str(),
                                  x_min, x.back());
  }

  if (x_max <= x.front()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "invalid -%s-range: %s_max <= min(%s) (%f <= %f).",
                                  axis.c_str(), axis.c_str(), axis.c_str(),
                                  x_max, x.front());
  }
}

static void subset_extent(const std::string& axis,
                          const std::vector<double> &x,
                          double x_min, double x_max,
                          double &x0, double &Lx, unsigned int &Mx) {

  validate_range(axis, x, x_min, x_max);

  size_t x_start = gsl_interp_bsearch(&x[0], x_min, 0, x.size() - 1);
  // include one more point if we can
  if (x_start > 0) {
    x_start -= 1;
  }

  size_t x_end = gsl_interp_bsearch(&x[0], x_max, 0, x.size() - 1);
  // include one more point if we can
  x_end = std::min(x.size() - 1, x_end + 1);

  Lx = (x[x_end] - x[x_start]) / 2.0;

  x0 = (x[x_start] + x[x_end]) / 2.0;

  assert(x_end + 1 >= x_start);
  Mx = x_end + 1 - x_start;
}

//! Create a grid using command-line options and (possibly) an input file.
/** Processes options -i, -bootstrap, -Mx, -My, -Mz, -Lx, -Ly, -Lz, -x_range, -y_range.
 */
IceGrid::Ptr regional_grid_from_options(Context::Ptr ctx) {

  const Periodicity p = string_to_periodicity(ctx->config()->get_string("grid.periodicity"));

  const options::String input_file("-i", "Specifies a PISM input file");
  const bool bootstrap = options::Bool("-bootstrap", "enable bootstrapping heuristics");
  const options::RealList x_range("-x_range",
                                  "range of X coordinates in the selected subset");
  const options::RealList y_range("-y_range",
                                  "range of Y coordinates in the selected subset");

  if (input_file.is_set() and bootstrap and x_range.is_set() and y_range.is_set()) {
    // bootstrapping; get domain size defaults from an input file, allow overriding all grid
    // parameters using command-line options

    if (x_range->size() != 2) {
      throw RuntimeError(PISM_ERROR_LOCATION, "invalid -x_range argument: need 2 numbers.");
    }

    if (y_range->size() != 2) {
      throw RuntimeError(PISM_ERROR_LOCATION, "invalid -y_range argument: need 2 numbers.");
    }

    GridParameters input_grid(ctx->config());

    std::vector<std::string> names;
    names.push_back("enthalpy");
    names.push_back("temp");
    names.push_back("land_ice_thickness");
    names.push_back("bedrock_altitude");
    names.push_back("thk");
    names.push_back("topg");
    bool grid_info_found = false;

    PIO file(ctx->com(), "netcdf3", input_file, PISM_READONLY);
    for (unsigned int i = 0; i < names.size(); ++i) {

      grid_info_found = file.inq_var(names[i]);
      if (not grid_info_found) {
        std::string dummy1;
        bool dummy2;
        // Failed to find using a short name. Try using names[i] as a
        // standard name...
        file.inq_var("dummy", names[i], grid_info_found, dummy1, dummy2);
      }

      if (grid_info_found) {
        input_grid = GridParameters(ctx, file, names[i], p);

        grid_info full = grid_info(file, names[i], ctx->unit_system(), p);

        // x direction
        subset_extent("x", full.x, x_range[0], x_range[1],
                      input_grid.x0, input_grid.Lx, input_grid.Mx);
        // y direction
        subset_extent("y", full.y, y_range[0], y_range[1],
                      input_grid.y0, input_grid.Ly, input_grid.My);

        // Set periodicity to "NONE" to that IceGrid computes coordinates correctly.
        input_grid.periodicity = NONE;

        break;
      }
    }

    if (not grid_info_found) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "no geometry information found in '%s'",
                                    input_file->c_str());
    }

    // ignore -Lx, -Ly, -Mx, -My
    options::ignored(*ctx->log(), "-Mx");
    options::ignored(*ctx->log(), "-My");
    options::ignored(*ctx->log(), "-Lx");
    options::ignored(*ctx->log(), "-Ly");

    // process options controlling vertical grid parameters, overriding values read from a file
    input_grid.vertical_grid_from_options(ctx->config());

    // process options controlling ownership ranges
    input_grid.ownership_ranges_from_options(ctx->size());

    return IceGrid::Ptr(new IceGrid(ctx, input_grid));
  } else if (x_range.is_set() ^ y_range.is_set()) {
    throw RuntimeError(PISM_ERROR_LOCATION, "Please set both -x_range and -y_range.");
  } else {
    return IceGrid::FromOptions(ctx);
  }
}


} // end of namespace pism
