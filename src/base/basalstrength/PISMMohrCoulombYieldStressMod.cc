// Copyright (C) 2004--2016 PISM Authors
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

#include <cmath>
#include <cassert>
#include <gsl/gsl_math.h>

#include "PISMMohrCoulombYieldStress.hh"
#include "PISMMohrCoulombYieldStressMod.hh"
#include "base/hydrology/PISMHydrology.hh"
#include "base/hydrology/PISMHydrologyMod.hh"
#include "base/util/IceGrid.hh"
#include "base/util/Mask.hh"
#include "base/util/PISMVars.hh"
#include "base/util/error_handling.hh"
#include "base/util/io/PIO.hh"
#include "base/util/pism_options.hh"
#include "base/util/MaxTimestep.hh"
#include "base/util/pism_utilities.hh"
#include "base/util/IceModelVec2CellType.hh"

namespace pism {

MohrCoulombYieldStressMod::MohrCoulombYieldStressMod(IceGrid::ConstPtr g,
                                               hydrology::Hydrology *hydro)
  : MohrCoulombYieldStress(g, hydro), m_hydrology(dynamic_cast<hydrology::HydrologyMod*>(hydro)) {


// in addition to the variables created from the base MohrCoulombYieldStress, there needs to be variables for:

// the precentage of area that has channels
// the percentage of area that is till covered

//  unsigned int stencil_width = m_config->get_double("grid_max_stencil_width");





/*
  const unsigned int stencil_width = m_config->get_double("grid.max_stencil_width");

  m_till_phi.create(m_grid, "tillphi", WITH_GHOSTS, stencil_width);
  m_till_phi.set_attrs("model_state",
                       "friction angle for till under grounded ice sheet",
                       "degrees", "");
  m_till_phi.set_time_independent(true);
  // in this model; need not be time-independent in general

  // internal working space; stencil width needed because redundant computation
  // on overlaps
  m_tillwat.create(m_grid, "tillwat_for_MohrCoulomb", WITH_GHOSTS, stencil_width);
  m_tillwat.set_attrs("internal",
                      "copy of till water thickness held by MohrCoulombYieldStress",
                      "m", "");

  m_Po.create(m_grid, "overburden_pressure_for_MohrCoulomb",
              WITH_GHOSTS, stencil_width);
  m_Po.set_attrs("internal",
                 "copy of overburden pressure held by MohrCoulombYieldStress",
                 "Pa", "");

  if (m_config->get_boolean("basal_yield_stress.add_transportable_water")) {
    m_bwat.create(m_grid, "bwat_for_MohrCoulomb", WITHOUT_GHOSTS);
    m_bwat.set_attrs("internal",
                     "copy of transportable water thickness held by MohrCoulombYieldStress",
                     "m", "");
  }
*/

  m_fraction_till_internal.create(m_grid, "fraction_till_for_MohrCoulomb",
              WITHOUT_GHOSTS);
  m_fraction_till_internal.set_attrs("internal",
                 "copy of fraction till held by MohrCoulombYieldStressMod",
                 "1", "");

  m_fraction_channel_internal.create(m_grid, "fraction_channel_for_MohrCoulomb",
              WITHOUT_GHOSTS);
  m_fraction_channel_internal.set_attrs("internal",
                 "copy of fraction channel held by MohrCoulombYieldStressMod",
                 "1", "");
}

MohrCoulombYieldStressMod::~MohrCoulombYieldStressMod() {
  // empty
}

void MohrCoulombYieldStressMod::init_impl() {

  // call the based MohrCoulombYieldStress to initialize the shear friction angle
  MohrCoulombYieldStress::init_impl();


/* ******************


problem here

sorry for the mess


******************* */
    hydrology::HydrologyMod *hydrology_Mod = dynamic_cast<hydrology::HydrologyMod*>(m_hydrology);
    if (hydrology_Mod == NULL) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "The the Mohr-Coulomb yield stress model needs a hydrology::HydrologyMod\n"
                                    "object. The current Hydrology instance is not suitable.  Set flag\n"
                                    "hydrology_model to 'hydrologymod'.");
    }

  // initialize the sediment distrbution

  

}

void MohrCoulombYieldStressMod::define_model_state_impl(const PIO &output) const {
  MohrCoulombYieldStress::define_model_state_impl(output);
//  m_fraction_till.define(output);
}

void MohrCoulombYieldStressMod::write_model_state_impl(const PIO &output) const {
  MohrCoulombYieldStress::write_model_state_impl(output);
//  m_fraction_till.write(output);
}



// override the MohrCoulombYieldStress updating function so that we can include a crude coupling 
// with basal hydrology and changes in the distribution of till
void MohrCoulombYieldStressMod::update_impl() {

  bool slippery_grounding_lines = m_config->get_boolean("basal_yield_stress.slippery_grounding_lines"),
       add_transportable_water  = m_config->get_boolean("basal_yield_stress.add_transportable_water");

  hydrology::HydrologyMod* hydroMod = dynamic_cast<hydrology::HydrologyMod*>(m_hydrology);
 // hydrology::Routing* hydrowithtransport = dynamic_cast<hydrology::Routing*>(m_hydrology);
//  if (hydroMod) {
//    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
//                                  "Cannot use mohr_coulomb_mod without also using hydrologymod");
//  }

  const double high_tauc   = m_config->get_double("basal_yield_stress.ice_free_bedrock"),
               tillwat_max = m_config->get_double("hydrology.tillwat_max"),
               c0          = m_config->get_double("basal_yield_stress.mohr_coulomb.till_cohesion"),
               N0          = m_config->get_double("basal_yield_stress.mohr_coulomb.till_reference_effective_pressure"),
               e0overCc    = (m_config->get_double("basal_yield_stress.mohr_coulomb.till_reference_void_ratio") /
                              m_config->get_double("basal_yield_stress.mohr_coulomb.till_compressibility_coefficient")),
               delta       = m_config->get_double("basal_yield_stress.mohr_coulomb.till_effective_fraction_overburden"),
               tlftw       = m_config->get_double("basal_yield_stress.mohr_coulomb.till_log_factor_transportable_water");

  if (m_hydrology) {
    m_hydrology->till_water_thickness(m_tillwat);
    m_hydrology->overburden_pressure(m_Po);

    m_hydrology -> fraction_till(m_fraction_till_internal);
    m_hydrology -> fraction_channel(m_fraction_channel_internal);

    if (add_transportable_water) {
      hydroMod->subglacial_water_thickness(m_bwat);
    }
  }

  const IceModelVec2CellType &mask           = *m_grid->variables().get_2d_cell_type("mask");
  const IceModelVec2S        &bed_topography = *m_grid->variables().get_2d_scalar("bedrock_altitude");

  IceModelVec::AccessList list{&m_tillwat, &m_till_phi, &m_basal_yield_stress, &mask,
      &bed_topography, &m_Po, &m_fraction_till_internal, &m_fraction_channel_internal};
  if (add_transportable_water) {
    list.add(m_bwat);
  }

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (mask.ocean(i, j)) {
      m_basal_yield_stress(i, j) = 0.0;
    } else if (mask.ice_free(i, j)) {
      m_basal_yield_stress(i, j) = high_tauc;  // large yield stress if grounded and ice-free
    } else { // grounded and there is some ice
      // user can ask that marine grounding lines get special treatment
      const double sea_level = 0.0; // FIXME: get sea-level from correct PISM source
      double water = m_tillwat(i,j); // usual case
      if (slippery_grounding_lines &&
          bed_topography(i,j) <= sea_level &&
          (mask.next_to_floating_ice(i,j) || mask.next_to_ice_free_ocean(i,j))) {
        water = tillwat_max;
      } else if (add_transportable_water) {
      //  water = m_tillwat(i,j) + tlftw * log(1.0 + m_bwat(i,j) / tlftw); // already done in PISMHydrologyMod
      }
      double
        s    = water / tillwat_max,
        Ntil = N0 * pow(delta * m_Po(i,j) / N0, s) * pow(10.0, e0overCc * (1.0 - s));
      Ntil = std::min(m_Po(i,j), Ntil);

      m_basal_yield_stress(i, j) = (high_tauc * (1.0 - m_fraction_till_internal(i,j)) + c0 + Ntil * tan((M_PI/180.0) * m_till_phi(i, j)) * m_fraction_till_internal(i,j)) * (1.0 - m_fraction_channel_internal(i,j));
    }
  }

  m_basal_yield_stress.update_ghosts();
}

} // end namespace pism

