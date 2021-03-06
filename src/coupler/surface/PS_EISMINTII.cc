/* Copyright (C) 2014, 2015 PISM Authors
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

#include "PS_EISMINTII.hh"
#include "coupler/PISMAtmosphere.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/pism_options.hh"
#include "base/util/error_handling.hh"
#include "base/util/IceGrid.hh"
#include "base/util/MaxTimestep.hh"

namespace pism {
namespace surface {

EISMINTII::EISMINTII(IceGrid::ConstPtr g, int experiment)
  : PSFormulas(g), m_experiment(experiment) {
  // empty
}

EISMINTII::~EISMINTII() {
  // empty
}

void EISMINTII::init_impl() {

  m_log->message(2,
             "setting parameters for surface mass balance"
             " and temperature in EISMINT II experiment %c ... \n",
             m_experiment);

  // EISMINT II specified values for parameters
  m_S_b = units::convert(m_sys, 1.0e-2 * 1e-3, "1/year", "1/s"); // Grad of accum rate change
  m_S_T = 1.67e-2 * 1e-3;       // K/m  Temp gradient

  // these are for A,E,G,H,I,K:
  m_M_max = units::convert(m_sys, 0.5, "m/year", "m/s"); // Max accumulation
  m_R_el  = 450.0e3;            // Distance to equil line (SMB=0)
  m_T_min = 238.15;

  switch (m_experiment) {
  case 'B':                     // supposed to start from end of experiment A and:
    m_T_min = 243.15;
    break;
  case 'C':
  case 'J':
  case 'L':                     // supposed to start from end of experiment A (for C;
    //   resp I and K for J and L) and:
    m_M_max = units::convert(m_sys, 0.25, "m/year", "m/s");
    m_R_el  = 425.0e3;
    break;
  case 'D':                     // supposed to start from end of experiment A and:
    m_R_el  = 425.0e3;
    break;
  case 'F':                     // start with zero ice and:
    m_T_min = 223.15;
    break;
  }

  // if user specifies Tmin, Tmax, Mmax, Sb, ST, Rel, then use that (override above)
  m_T_min = options::Real("-Tmin", "T min, Kelvin", m_T_min);

  options::Real Mmax("-Mmax", "Maximum accumulation, m/year",
                     units::convert(m_sys, m_M_max, "m/s", "m/year"));
  if (Mmax.is_set()) {
    m_M_max = units::convert(m_sys, Mmax, "m/year", "m/second");
  }

  options::Real Sb("-Sb", "radial gradient of accumulation rate, (m/year)/km",
                   units::convert(m_sys, m_S_b,   "m/second/m", "m/year/km"));
  if (Sb.is_set()) {
    m_S_b = units::convert(m_sys, Sb, "m/year/km", "m/second/m");
  }

  options::Real ST("-ST", "radial gradient of surface temperature, K/km",
                   units::convert(m_sys, m_S_T, "K/m", "K/km"));
  if (ST.is_set()) {
    m_S_T = units::convert(m_sys, ST, "K/km", "K/m");
  }

  options::Real Rel("-Rel", "radial distance to equilibrium line, km",
                    units::convert(m_sys, m_R_el, "m", "km"));
  if (Rel.is_set()) {
    m_R_el = units::convert(m_sys, Rel, "km", "m");
  }

  initialize_using_formulas();
}

MaxTimestep EISMINTII::max_timestep_impl(double t) {
  (void) t;
  return MaxTimestep();
}

void EISMINTII::initialize_using_formulas() {

  PetscScalar cx = m_grid->Lx(), cy = m_grid->Ly();
  if (m_experiment == 'E') {
    // shift center
    cx += 100.0e3;
    cy += 100.0e3;
  }

  // now fill in accum and surface temp
  IceModelVec::AccessList list;
  list.add(m_ice_surface_temp);
  list.add(m_climatic_mass_balance);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // r is distance from center of grid; if E then center is shifted (above)
    const double r = sqrt(PetscSqr(-cx + m_grid->dx()*i)
                          + PetscSqr(-cy + m_grid->dy()*j));
    // set accumulation from formula (7) in (Payne et al 2000)
    m_climatic_mass_balance(i,j) = std::min(m_M_max, m_S_b * (m_R_el-r));
    // set surface temperature
    m_ice_surface_temp(i,j) = m_T_min + m_S_T * r;  // formula (8) in (Payne et al 2000)
  }

  // convert from [m/s] to [kg m-2 s-1]
  m_climatic_mass_balance.scale(m_config->get_double("ice_density"));
}

void EISMINTII::update_impl(PetscReal t, PetscReal dt) {
  (void) t;
  (void) dt;

  // do nothing (but an implementation is required)
}

} // end of namespace surface
} // end of namespace pism
