/* Copyright (C) 2016 PISM Authors
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

#include "EnergyModel.hh"
#include "base/util/MaxTimestep.hh"
#include "base/stressbalance/PISMStressBalance.hh"
#include "base/util/io/PIO.hh"
#include "base/util/PISMVars.hh"
#include "utilities.hh"
#include "base/enthalpyConverter.hh"
#include "bootstrapping.hh"

namespace pism {
namespace energy {

static void check_input(const IceModelVec *ptr, const char *name) {
  if (ptr == NULL) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "energy balance model input %s was not provided", name);
  }
}

EnergyModelInputs::EnergyModelInputs() {
  basal_frictional_heating = NULL;
  basal_heat_flux          = NULL;
  cell_type                = NULL;
  ice_thickness            = NULL;
  surface_liquid_fraction  = NULL;
  shelf_base_temp          = NULL;
  surface_temp             = NULL;
  till_water_thickness     = NULL;

  strain_heating3          = NULL;
  u3                       = NULL;
  v3                       = NULL;
  w3                       = NULL;
}

void EnergyModelInputs::check() const {
  check_input(cell_type,                "cell_type");
  check_input(basal_frictional_heating, "basal_frictional_heating");
  check_input(basal_heat_flux,          "basal_heat_flux");
  check_input(ice_thickness,            "ice_thickness");
  check_input(surface_liquid_fraction,  "surface_liquid_fraction");
  check_input(shelf_base_temp,          "shelf_base_temp");
  check_input(surface_temp,             "surface_temp");
  check_input(till_water_thickness,     "till_water_thickness");

  check_input(strain_heating3, "strain_heating3");
  check_input(u3, "u3");
  check_input(v3, "v3");
  check_input(w3, "w3");
}

EnergyModelStats::EnergyModelStats() {
  bulge_counter = 0;
  reduced_accuracy_counter = 0;
  low_temperature_counter = 0;
  liquified_ice_volume = 0.0;
}

EnergyModel::EnergyModel(IceGrid::ConstPtr grid,
                         stressbalance::StressBalance *stress_balance)
  : Component_TS(grid), m_stress_balance(stress_balance) {

  const unsigned int WIDE_STENCIL = m_config->get_double("grid.max_stencil_width");

  {
    m_ice_enthalpy.create(m_grid, "enthalpy", WITH_GHOSTS, WIDE_STENCIL);
    // POSSIBLE standard name = land_ice_enthalpy
    m_ice_enthalpy.set_attrs("model_state",
                             "ice enthalpy (includes sensible heat, latent heat, pressure)",
                             "J kg-1", "");
  }

  {
    m_basal_melt_rate.create(m_grid, "bmelt", WITHOUT_GHOSTS);
    // ghosted to allow the "redundant" computation of tauc
    m_basal_melt_rate.set_attrs("model_state",
                                "ice basal melt rate from energy conservation and subshelf melt, in ice thickness per time",
                                "m s-1", "land_ice_basal_melt_rate");
    m_basal_melt_rate.metadata().set_string("glaciological_units", "m year-1");
    m_basal_melt_rate.write_in_glaciological_units = true;
    m_basal_melt_rate.metadata().set_string("comment", "positive basal melt rate corresponds to ice loss");
  }

  // a 3d work vector
  {
    m_work.create(m_grid, "work_vector", WITHOUT_GHOSTS);
    m_work.set_attrs("internal",
                     "usually new values of temperature or enthalpy during time step",
                     "", "");
  }
}

void EnergyModel::init_enthalpy(const PIO &input_file, bool do_regrid, int record) {

  const IceModelVec2S & ice_thickness = *m_grid->variables().get_2d_scalar("land_ice_thickness");

  const bool
    enthalpy_exists = input_file.inq_var("enthalpy"),
    temp_exists     = input_file.inq_var("temp"),
    liqfrac_exists  = input_file.inq_var("liqfrac");

  if (enthalpy_exists) {
    if (do_regrid) {
      m_ice_enthalpy.regrid(input_file, CRITICAL);
    } else {
      m_ice_enthalpy.read(input_file, record);
    }
  } else if (temp_exists) {
    IceModelVec3
      &temp    = m_work,
      &liqfrac = m_ice_enthalpy;

    SpatialVariableMetadata enthalpy_metadata = m_ice_enthalpy.metadata();
    temp.set_name("temp");
    temp.metadata(0).set_name("temp");
    temp.set_attrs("temporary", "ice temperature", "Kelvin",
                   "land_ice_temperature");

    if (do_regrid) {
      temp.regrid(input_file, CRITICAL);
    } else {
      temp.read(input_file, record);
    }

    if (liqfrac_exists and not m_config->get_boolean("energy.temperature_based")) {
      liqfrac.set_name("liqfrac");
      liqfrac.metadata(0).set_name("liqfrac");
      liqfrac.set_attrs("temporary", "ice liquid water fraction",
                        "1", "");

      if (do_regrid) {
        liqfrac.regrid(input_file, CRITICAL);
      } else {
        liqfrac.read(input_file, record);
      }

      m_log->message(2,
                     "* Computing enthalpy using ice temperature,"
                     "  liquid water fraction and thickness...\n");
      compute_enthalpy(temp, liqfrac, ice_thickness, m_ice_enthalpy);
    } else {
      m_log->message(2,
                     "* Computing enthalpy using ice temperature and thickness "
                     "and assuming zero liquid water fraction...\n");
      compute_enthalpy_cold(temp, ice_thickness, m_ice_enthalpy);
    }

    m_ice_enthalpy.metadata() = enthalpy_metadata;
  } else {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "neither enthalpy nor temperature was found in '%s'.\n",
                                  input_file.inq_filename().c_str());
  }
}


void EnergyModel::init(const InputOptions &opts) {
  this->init_impl(opts);
}

void EnergyModel::update(double t, double dt, const EnergyModelInputs &inputs) {
  this->update_impl(t, dt, inputs);
}

void EnergyModel::update_impl(double t, double dt) {
  (void) t;
  (void) dt;
  throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                "EnergyModel::update_impl(t, dt) is not implemented");
}

MaxTimestep EnergyModel::max_timestep_impl(double t) const {
  // fix a compiler warning
  (void) t;

  if (m_stress_balance == NULL) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "EnergyModel: no stress balance provided."
                                  " Cannot compute max. time step.");
  }

  return m_stress_balance->max_timestep_cfl_3d().dt_max;
}

const EnergyModelStats& EnergyModel::stats() const {
  return m_stats;
}

const IceModelVec3 & EnergyModel::enthalpy() const {
  return m_ice_enthalpy;
}

const IceModelVec2S & EnergyModel::basal_melt_rate() const {
  return m_basal_melt_rate;
}

void EnergyModel::define_model_state_impl(const PIO &output) const {
  m_ice_enthalpy.define(output);
  m_basal_melt_rate.define(output);
}

void EnergyModel::write_model_state_impl(const PIO &output) const {
  m_ice_enthalpy.write(output);
  m_basal_melt_rate.write(output);
}

} // end of namespace energy

} // end of namespace pism