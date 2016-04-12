// Copyright (C) 2009--2016 Ed Bueler, Constantine Khroulev, and David Maxwell
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

#include "SSATestCase.hh"
#include "SSAFD.hh"
#include "SSAFEM.hh"
#include "base/util/Mask.hh"
#include "base/util/PISMTime.hh"
#include "base/util/io/PIO.hh"
#include "base/util/pism_options.hh"
#include "base/util/io/io_helpers.hh"
#include "base/util/pism_utilities.hh"

namespace pism {
namespace stressbalance {

//! Initialize the storage for the various coefficients used as input to the SSA
//! (ice elevation, thickness, etc.)
void SSATestCase::buildSSACoefficients()
{

  const unsigned int WIDE_STENCIL = m_config->get_double("grid_max_stencil_width");

  // ice surface elevation
  m_surface.create(m_grid, "usurf", WITH_GHOSTS, WIDE_STENCIL);
  m_surface.set_attrs("diagnostic", "ice upper surface elevation", "m",
                      "surface_altitude");
  m_grid->variables().add(m_surface);

  // land ice thickness
  m_thickness.create(m_grid, "thk", WITH_GHOSTS, WIDE_STENCIL);
  m_thickness.set_attrs("model_state", "land ice thickness", "m",
                        "land_ice_thickness");
  m_thickness.metadata().set_double("valid_min", 0.0);
  m_grid->variables().add(m_thickness);

  // bedrock surface elevation
  m_bed.create(m_grid, "topg", WITH_GHOSTS, WIDE_STENCIL);
  m_bed.set_attrs("model_state", "bedrock surface elevation", "m",
                  "bedrock_altitude");
  m_grid->variables().add(m_bed);

  // yield stress for basal till (plastic or pseudo-plastic model)
  m_tauc.create(m_grid, "tauc", WITH_GHOSTS, WIDE_STENCIL);
  m_tauc.set_attrs("diagnostic",
                   "yield stress for basal till (plastic or pseudo-plastic model)", "Pa", "");
  m_grid->variables().add(m_tauc);

  // enthalpy
  m_ice_enthalpy.create(m_grid, "enthalpy", WITH_GHOSTS, WIDE_STENCIL);
  m_ice_enthalpy.set_attrs("model_state",
                       "ice enthalpy (includes sensible heat, latent heat, pressure)",
                       "J kg-1", "");
  m_grid->variables().add(m_ice_enthalpy);


  // dirichlet boundary condition (FIXME: perhaps unused!)
  m_bc_values.create(m_grid, "_bc", WITH_GHOSTS, WIDE_STENCIL); // u_bc and v_bc
  m_bc_values.set_attrs("intent",
                     "X-component of the SSA velocity boundary conditions",
                     "m s-1", "", 0);
  m_bc_values.set_attrs("intent",
                     "Y-component of the SSA velocity boundary conditions",
                     "m s-1", "", 1);

  Config::ConstPtr config = m_grid->ctx()->config();
  units::System::Ptr sys = m_grid->ctx()->unit_system();
  double fill_value = units::convert(sys, config->get_double("fill_value"), "m year-1", "m second-1");

  m_bc_values.metadata(0).set_string("glaciological_units", "m year-1");
  m_bc_values.metadata(0).set_double("valid_min", units::convert(m_sys, -1e6, "m year-1", "m second-1"));
  m_bc_values.metadata(0).set_double("valid_max", units::convert(m_sys,  1e6, "m year-1", "m second-1"));
  m_bc_values.metadata(0).set_double("_FillValue", fill_value);

  m_bc_values.metadata(1).set_string("glaciological_units", "m year-1");
  m_bc_values.metadata(1).set_double("valid_min", units::convert(m_sys, -1e6, "m year-1", "m second-1"));
  m_bc_values.metadata(1).set_double("valid_max", units::convert(m_sys,  1e6, "m year-1", "m second-1"));
  m_bc_values.metadata(1).set_double("_FillValue", fill_value);

  m_bc_values.write_in_glaciological_units = true;
  m_bc_values.set(fill_value);

  // grounded_dragging_floating integer mask
  m_ice_mask.create(m_grid, "mask", WITH_GHOSTS, WIDE_STENCIL);
  m_ice_mask.set_attrs("model_state",
                       "grounded_dragging_floating integer mask", "", "");
  std::vector<double> mask_values(4);
  mask_values[0] = MASK_ICE_FREE_BEDROCK;
  mask_values[1] = MASK_GROUNDED;
  mask_values[2] = MASK_FLOATING;
  mask_values[3] = MASK_ICE_FREE_OCEAN;
  m_ice_mask.metadata().set_doubles("flag_values", mask_values);
  m_ice_mask.metadata().set_string("flag_meanings",
                                   "ice_free_bedrock grounded_ice floating_ice ice_free_ocean");
  m_grid->variables().add(m_ice_mask);

  m_ice_mask.set(MASK_GROUNDED);

  // Dirichlet B.C. mask
  m_bc_mask.create(m_grid, "bc_mask", WITH_GHOSTS, WIDE_STENCIL);
  m_bc_mask.set_attrs("model_state",
                      "grounded_dragging_floating integer mask", "", "");
  mask_values.resize(2);
  mask_values[0] = 0;
  mask_values[1] = 1;
  m_bc_mask.metadata().set_doubles("flag_values", mask_values);
  m_bc_mask.metadata().set_string("flag_meanings",
                                  "no_data ssa_dirichlet_bc_location");
  m_grid->variables().add(m_bc_mask);

  m_melange_back_pressure.create(m_grid, "melange_back_pressure_fraction",
                                 WITH_GHOSTS, WIDE_STENCIL);
  m_melange_back_pressure.set_attrs("boundary_condition",
                                    "melange back pressure fraction", "", "");
  m_melange_back_pressure.set(0.0);
}

SSATestCase::SSATestCase(Context::Ptr ctx)
  : m_com(ctx->com()), m_config(ctx->config()), m_ctx(ctx),
    m_sys(ctx->unit_system()), m_ssa(NULL) {
  // empty
}

SSATestCase::~SSATestCase()
{
  delete m_ssa;
}

//! Initialize the test case at the start of a run
void SSATestCase::init(int Mx, int My, SSAFactory ssafactory) {
  // Subclass builds grid->
  initializeGrid(Mx, My);

  // Subclass builds ice flow law, basal resistance, etc.
  initializeSSAModel();

  // We setup storage for the coefficients.
  buildSSACoefficients();

  // Allocate the actual SSA solver.
  m_ssa = ssafactory(m_grid, m_enthalpyconverter);
  m_ssa->init(); // vars was setup preivouisly with buildSSACoefficients

  // Allow the subclass to setup the coefficients.
  initializeSSACoefficients();
}

//! Solve the SSA
void SSATestCase::run() {
  // Solve (fast==true means "no update"):
  m_ctx->log()->message(2, "* Solving the SSA stress balance ...\n");

  bool fast = false;
  double sea_level = 0.0;
  m_ssa->update(fast, sea_level, m_melange_back_pressure);
}

//! Report on the generated solution
void SSATestCase::report(const std::string &testname) {

  m_ctx->log()->message(3, m_ssa->stdout_report());

  double
    maxvecerr  = 0.0,
    avvecerr   = 0.0,
    avuerr     = 0.0,
    avverr     = 0.0,
    maxuerr    = 0.0,
    maxverr    = 0.0;
  double
    gmaxvecerr = 0.0,
    gavvecerr  = 0.0,
    gavuerr    = 0.0,
    gavverr    = 0.0,
    gmaxuerr   = 0.0,
    gmaxverr   = 0.0;

  if (m_config->get_boolean("do_pseudo_plastic_till") &&
      m_config->get_double("pseudo_plastic_q") != 1.0) {
    m_ctx->log()->message(1,
                          "WARNING: numerical errors not valid for pseudo-plastic till\n");
  }
  m_ctx->log()->message(1,
                        "NUMERICAL ERRORS in velocity relative to exact solution:\n");

  const IceModelVec2V &vel_ssa = m_ssa->velocity();

  IceModelVec::AccessList list;
  list.add(vel_ssa);

  double exactvelmax = 0, gexactvelmax = 0;
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double uexact, vexact;
    double myx = m_grid->x(i), myy = m_grid->y(j);

    exactSolution(i,j,myx,myy,&uexact,&vexact);

    double exactnormsq=sqrt(uexact*uexact+vexact*vexact);
    exactvelmax = std::max(exactnormsq,exactvelmax);

    // compute maximum errors
    const double uerr = fabs(vel_ssa(i,j).u - uexact);
    const double verr = fabs(vel_ssa(i,j).v - vexact);
    avuerr = avuerr + uerr;
    avverr = avverr + verr;
    maxuerr = std::max(maxuerr,uerr);
    maxverr = std::max(maxverr,verr);
    const double vecerr = sqrt(uerr * uerr + verr * verr);
    maxvecerr = std::max(maxvecerr,vecerr);
    avvecerr = avvecerr + vecerr;
  }

  unsigned int N = (m_grid->Mx()*m_grid->My());

  gexactvelmax = GlobalMax(m_grid->com, exactvelmax);
  gmaxuerr     = GlobalMax(m_grid->com, maxuerr);
  gmaxverr     = GlobalMax(m_grid->com, maxverr);
  gavuerr      = GlobalSum(m_grid->com, avuerr);
  gavuerr      = gavuerr / N;
  gavverr      = GlobalSum(m_grid->com, avverr);
  gavverr      = gavverr / N;
  gmaxvecerr   = GlobalMax(m_grid->com, maxvecerr);
  gavvecerr    = GlobalSum(m_grid->com, avvecerr);
  gavvecerr    = gavvecerr / N;

  using pism::units::convert;

  m_ctx->log()->message(1,
                        "velocity  :  maxvector   prcntavvec      maxu      maxv       avu       avv\n");
  m_ctx->log()->message(1,
                        "           %11.4f%13.5f%10.4f%10.4f%10.4f%10.4f\n",
                        convert(m_sys, gmaxvecerr, "m second-1", "m year-1"),
                        (gavvecerr/gexactvelmax)*100.0,
                        convert(m_sys, gmaxuerr, "m second-1", "m year-1"),
                        convert(m_sys, gmaxverr, "m second-1", "m year-1"),
                        convert(m_sys, gavuerr, "m second-1", "m year-1"),
                        convert(m_sys, gavverr, "m second-1", "m year-1"));

  m_ctx->log()->message(1, "NUM ERRORS DONE\n");

  report_netcdf(testname,
                convert(m_sys, gmaxvecerr, "m second-1", "m year-1"),
                (gavvecerr/gexactvelmax)*100.0,
                convert(m_sys, gmaxuerr, "m second-1", "m year-1"),
                convert(m_sys, gmaxverr, "m second-1", "m year-1"),
                convert(m_sys, gavuerr, "m second-1", "m year-1"),
                convert(m_sys, gavverr, "m second-1", "m year-1"));
}

void SSATestCase::report_netcdf(const std::string &testname,
                                double max_vector,
                                double rel_vector,
                                double max_u,
                                double max_v,
                                double avg_u,
                                double avg_v) {
  TimeseriesMetadata err("N", "N", m_grid->ctx()->unit_system());
  unsigned int start;
  VariableMetadata global_attributes("PISM_GLOBAL", m_grid->ctx()->unit_system());

  options::String filename("-report_file", "NetCDF error report file");

  if (not filename.is_set()) {
    return;
  }

  err.set_string("units", "1");

  m_ctx->log()->message(2, "Also writing errors to '%s'...\n", filename->c_str());

  bool append = options::Bool("-append", "Append the NetCDF error report");

  IO_Mode mode = PISM_READWRITE;
  if (not append) {
    mode = PISM_READWRITE_MOVE;
  }

  global_attributes.set_string("source", std::string("PISM ") + PISM_Revision);

  // Find the number of records in this file:
  PIO nc(m_grid->com, "netcdf3");      // OK to use NetCDF3.
  nc.open(filename, mode);
  start = nc.inq_dimlen("N");

  io::write_global_attributes(nc, global_attributes);

  // Write the dimension variable:
  io::write_timeseries(nc, err, (size_t)start, (double)(start + 1), PISM_INT);

  // Always write grid parameters:
  err.set_name("dx");
  err.set_string("units", "meters");
  io::write_timeseries(nc, err, (size_t)start, m_grid->dx());
  err.set_name("dy");
  io::write_timeseries(nc, err, (size_t)start, m_grid->dy());

  // Always write the test name:
  err.clear_all_strings(); err.clear_all_doubles(); err.set_string("units", "1");
  err.set_name("test");
  io::write_timeseries(nc, err, (size_t)start, testname[0], PISM_BYTE);

  err.clear_all_strings(); err.clear_all_doubles(); err.set_string("units", "1");
  err.set_name("max_velocity");
  err.set_string("units", "m year-1");
  err.set_string("long_name", "maximum ice velocity magnitude error");
  io::write_timeseries(nc, err, (size_t)start, max_vector);

  err.clear_all_strings(); err.clear_all_doubles(); err.set_string("units", "1");
  err.set_name("relative_velocity");
  err.set_string("units", "percent");
  err.set_string("long_name", "relative ice velocity magnitude error");
  io::write_timeseries(nc, err, (size_t)start, rel_vector);

  err.clear_all_strings(); err.clear_all_doubles(); err.set_string("units", "1");
  err.set_name("maximum_u");
  err.set_string("units", "m year-1");
  err.set_string("long_name", "maximum error in the X-component of the ice velocity");
  io::write_timeseries(nc, err, (size_t)start, max_u);

  err.clear_all_strings(); err.clear_all_doubles(); err.set_string("units", "1");
  err.set_name("maximum_v");
  err.set_string("units", "m year-1");
  err.set_string("long_name", "maximum error in the Y-component of the ice velocity");
  io::write_timeseries(nc, err, (size_t)start, max_v);

  err.clear_all_strings(); err.clear_all_doubles(); err.set_string("units", "1");
  err.set_name("average_u");
  err.set_string("units", "m year-1");
  err.set_string("long_name", "average error in the X-component of the ice velocity");
  io::write_timeseries(nc, err, (size_t)start, avg_u);

  err.clear_all_strings(); err.clear_all_doubles(); err.set_string("units", "1");
  err.set_name("average_v");
  err.set_string("units", "m year-1");
  err.set_string("long_name", "average error in the Y-component of the ice velocity");
  io::write_timeseries(nc, err, (size_t)start, avg_v);

  nc.close();
}

void SSATestCase::exactSolution(int /*i*/, int /*j*/,
                                double /*x*/, double /*y*/,
                                double *u, double *v) {
  *u=0; *v=0;
}

//! Save the computation and data to a file.
void SSATestCase::write(const std::string &filename) {

  // Write results to an output file:
  PIO pio(m_grid->com, m_grid->ctx()->config()->get_string("output_format"));
  pio.open(filename, PISM_READWRITE_MOVE);
  io::define_time(pio, m_config->get_string("time_dimension_name"),
                  m_grid->ctx()->time()->calendar(),
                  m_grid->ctx()->time()->CF_units_string(),
                  m_grid->ctx()->unit_system());
  io::append_time(pio, m_config->get_string("time_dimension_name"), 0.0);

  m_surface.write(pio);
  m_thickness.write(pio);
  m_bc_mask.write(pio);
  m_tauc.write(pio);
  m_bed.write(pio);
  m_ice_enthalpy.write(pio);
  m_bc_values.write(pio);

  const IceModelVec2V &vel_ssa = m_ssa->velocity();
  vel_ssa.write(pio);

  IceModelVec2V exact;
  exact.create(m_grid, "_exact", WITHOUT_GHOSTS);
  exact.set_attrs("diagnostic",
                  "X-component of the SSA exact solution",
                  "m s-1", "", 0);
  exact.set_attrs("diagnostic",
                  "Y-component of the SSA exact solution",
                  "m s-1", "", 1);
  exact.metadata(0).set_string("glaciological_units", "m year-1");
  exact.metadata(1).set_string("glaciological_units", "m year-1");
  exact.write_in_glaciological_units = true;

  IceModelVec::AccessList list(exact);
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    exactSolution(i, j, m_grid->x(i), m_grid->y(j),
                  &(exact(i,j).u), &(exact(i,j).v));
  }
  exact.write(pio);

  pio.close();
}

} // end of namespace stressbalance
} // end of namespace pism
