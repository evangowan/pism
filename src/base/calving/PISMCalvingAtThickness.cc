/* Copyright (C) 2013, 2014, 2015, 2016 PISM Authors
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

#include "PISMCalvingAtThickness.hh"
#include "base/util/Mask.hh"
#include "base/util/error_handling.hh"
#include "base/util/IceGrid.hh"
#include "base/util/pism_const.hh"

namespace pism {

//! @brief Calving and iceberg removal code.
namespace calving {

CalvingAtThickness::CalvingAtThickness(IceGrid::ConstPtr g)
  : Component(g) {
  m_calving_threshold = m_config->get_double("thickness_calving_threshold");

  m_old_mask.create(m_grid, "old_mask", WITH_GHOSTS, 1);
}

CalvingAtThickness::~CalvingAtThickness() {
  // empty
}


void CalvingAtThickness::init() {
  m_log->message(2,
             "* Initializing the 'calving at a threshold thickness' mechanism...\n"
             "  thickness threshold: %3.3f meters\n", m_calving_threshold);
}

/**
 * Updates ice cover mask and the ice thickness according to the
 * calving rule removing ice at the shelf front that is thinner than a
 * given threshold.
 *
 * @param[in,out] pism_mask ice cover mask
 * @param[in,out] ice_thickness ice thickness
 *
 * @return 0 on success
 */
void CalvingAtThickness::update(IceModelVec2CellType &pism_mask,
                                IceModelVec2S &ice_thickness) {

  // this call fills ghosts of m_old_mask
  m_old_mask.copy_from(pism_mask);

  IceModelVec::AccessList list;
  list.add(pism_mask);
  list.add(ice_thickness);
  list.add(m_old_mask);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (m_old_mask.floating_ice(i, j)           &&
        m_old_mask.next_to_ice_free_ocean(i, j) &&
        ice_thickness(i, j) < m_calving_threshold) {
      ice_thickness(i, j) = 0.0;
      pism_mask(i, j)     = MASK_ICE_FREE_OCEAN;
    }
  }

  pism_mask.update_ghosts();
  ice_thickness.update_ghosts();
}


void CalvingAtThickness::add_vars_to_output_impl(const std::string &/*keyword*/,
                                                std::set<std::string> &/*result*/) {
  // empty
}

void CalvingAtThickness::define_variables_impl(const std::set<std::string> &/*vars*/,
                                                        const PIO &/*nc*/,
                                                        IO_Type /*nctype*/) {
  // empty
}

void CalvingAtThickness::write_variables_impl(const std::set<std::string> &/*vars*/,
                                                       const PIO& /*nc*/) {
  // empty
}


} // end of namespace calving
} // end of namespace pism
