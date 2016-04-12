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

#include "PISMEigenCalving.hh"
#include "base/stressbalance/PISMStressBalance.hh"
#include "base/util/IceGrid.hh"
#include "base/util/Mask.hh"
#include "base/util/PISMVars.hh"
#include "base/util/error_handling.hh"
#include "base/util/pism_const.hh"
#include "base/util/MaxTimestep.hh"
#include "base/util/pism_utilities.hh"
#include "base/util/IceModelVec2CellType.hh"

namespace pism {
namespace calving {

EigenCalving::EigenCalving(IceGrid::ConstPtr g,
                           stressbalance::StressBalance *stress_balance)
  : Component(g), m_stencil_width(2),
    m_stress_balance(stress_balance) {
  m_strain_rates.create(m_grid, "strain_rates", WITH_GHOSTS,
                        m_stencil_width,
                        2);

  m_thk_loss.create(m_grid, "temporary_storage", WITH_GHOSTS, 1);

  m_strain_rates.metadata(0).set_name("eigen1");
  m_strain_rates.set_attrs("internal",
                           "major principal component of horizontal strain-rate",
                           "second-1", "", 0);

  m_strain_rates.metadata(1).set_name("eigen2");
  m_strain_rates.set_attrs("internal",
                           "minor principal component of horizontal strain-rate",
                           "second-1", "", 1);

  m_K = m_config->get_double("eigen_calving_K");
  m_restrict_timestep = m_config->get_boolean("cfl_eigen_calving");
}

EigenCalving::~EigenCalving() {
  // empty
}

void EigenCalving::init() {

  m_log->message(2,
             "* Initializing the 'eigen-calving' mechanism...\n");

  if (fabs(m_grid->dx() - m_grid->dy()) / std::min(m_grid->dx(), m_grid->dy()) > 1e-2) {
    throw RuntimeError::formatted("-calving eigen_calving using a non-square grid cell is not implemented (yet);\n"
                                  "dx = %f, dy = %f, relative difference = %f",
                                  m_grid->dx(), m_grid->dy(),
                                  fabs(m_grid->dx() - m_grid->dy()) / std::max(m_grid->dx(), m_grid->dy()));
  }

  m_strain_rates.set(0.0);

}

//! \brief Uses principal strain rates to apply "eigencalving" with constant K.
/*!
  See equation (26) in [\ref Winkelmannetal2011].
*/
void EigenCalving::update(double dt,
                          IceModelVec2CellType &mask,
                          IceModelVec2S &Href,
                          IceModelVec2S &ice_thickness) {

  // Distance (grid cells) from calving front where strain rate is evaluated
  int offset = m_stencil_width;

  m_thk_loss.set(0.0);

  update_strain_rates();

  IceModelVec::AccessList list;
  list.add(ice_thickness);
  list.add(mask);
  list.add(Href);
  list.add(m_strain_rates);
  list.add(m_thk_loss);

  for (Points pt(*m_grid); pt; pt.next()) {
    const int i = pt.i(), j = pt.j();
    // Average of strain-rate eigenvalues in adjacent floating grid
    // cells to be used for eigen-calving:
    double
      eigen1    = 0.0,
      eigen2    = 0.0;
    // Neighbor-averaged ice thickness:
    double
      H_average = 0.0;

    int N_floating_neighbors = 0, M = 0;
    // M is the number of cells used to compute averages of strain
    // rate components.

    // Find partially filled or empty grid boxes on the icefree ocean, which
    // have floating ice neighbors after the mass continuity step
    if (mask.ice_free_ocean(i, j) &&
        mask.next_to_floating_ice(i, j)) {

      if (mask.floating_ice(i + 1, j)) {
        N_floating_neighbors += 1;
        H_average += ice_thickness(i + 1, j);
      }

      if (mask.floating_ice(i - 1, j)) {
        N_floating_neighbors += 1;
        H_average += ice_thickness(i - 1, j);
      }

      if (mask.floating_ice(i, j + 1)) {
        N_floating_neighbors += 1;
        H_average += ice_thickness(i, j + 1);
      }

      if (mask.floating_ice(i, j - 1)) {
        N_floating_neighbors += 1;
        H_average += ice_thickness(i, j - 1);
      }

      if (N_floating_neighbors > 0) {
        H_average /= N_floating_neighbors;
      }

      for (int p = -1; p < 2; p += 2) {
        int i_offset = p * offset;
        if (mask.floating_ice(i + i_offset, j) &&
            not mask.ice_margin(i + i_offset, j)) {
          eigen1 += m_strain_rates(i + i_offset, j, 0);
          eigen2 += m_strain_rates(i + i_offset, j, 1);
          M += 1;
        }
      }

      for (int q = -1; q < 2; q += 2) {
        int j_offset = q * offset;
        if (mask.floating_ice(i, j + j_offset) &&
            not mask.ice_margin(i , j + j_offset)) {
          eigen1 += m_strain_rates(i, j + j_offset, 0);
          eigen2 += m_strain_rates(i, j + j_offset, 1);
          M += 1;
        }
      }

      if (M > 0) {
        eigen1 /= M;
        eigen2 /= M;
      }

      double
        calving_rate_horizontal = 0.0,
        eigenCalvOffset         = 0.0;
      // eigenCalvOffset allows adjusting the transition from
      // compressive to extensive flow regime

      // Calving law
      //
      // eigen1 * eigen2 has units [s^-2] and calving_rate_horizontal
      // [m*s^1] hence, eigen_calving_K has units [m*s]
      if (eigen2 > eigenCalvOffset and eigen1 > 0.0) {
        // spreading in all directions
        calving_rate_horizontal = m_K * eigen1 * (eigen2 - eigenCalvOffset);
      }

      // calculate mass loss with respect to the associated ice thickness and the grid size:
      double calving_rate = calving_rate_horizontal * H_average / m_grid->dx(); // in m/s

      // apply calving rate at partially filled or empty grid cells
      if (calving_rate > 0.0) {
        Href(i, j) -= calving_rate * dt; // in m

        if (Href(i, j) < 0.0) {
          // Partially filled grid cell became ice-free

          m_thk_loss(i, j) = -Href(i, j); // in m, corresponds to additional ice loss
          Href(i, j)       = 0.0;

          // additional mass loss will be distributed among
          // N_floating_neighbors:
          if (N_floating_neighbors > 0) {
            m_thk_loss(i, j) /= N_floating_neighbors;
          }
        }
      }

    } // end of "if (ice_free_ocean && next_to_floating)"
  }

  m_thk_loss.update_ghosts();

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    double thk_loss_ij = 0.0;

    if (mask.floating_ice(i, j) &&
        (m_thk_loss(i + 1, j) > 0.0 || m_thk_loss(i - 1, j) > 0.0 ||
         m_thk_loss(i, j + 1) > 0.0 || m_thk_loss(i, j - 1) > 0.0)) {

      thk_loss_ij = (m_thk_loss(i + 1, j) + m_thk_loss(i - 1, j) +
                     m_thk_loss(i, j + 1) + m_thk_loss(i, j - 1));     // in m/s

      // Note std::max: we do not account for further calving
      // ice-inwards! Alternatively CFL criterion for time stepping
      // could be adjusted to maximum of calving rate
      Href(i, j) = std::max(ice_thickness(i, j) - thk_loss_ij, 0.0); // in m

      ice_thickness(i, j) = 0.0;
      mask(i, j) = MASK_ICE_FREE_OCEAN;
    }
  }

  mask.update_ghosts();

  remove_narrow_tongues(mask, ice_thickness);

  ice_thickness.update_ghosts();
  mask.update_ghosts();
}


/**
 * @brief Compute the maximum time-step length allowed by the CFL
 * condition applied to the calving rate.
 *
 * Note: this code uses the mask variable obtained from the Vars
 * dictionary. This is not the same mask that is used in the update()
 * call, since max_timestep() is called *before* the mass-continuity
 * time-step.
 *
 * @return 0 on success
 */
MaxTimestep EigenCalving::max_timestep() {

  if (not m_restrict_timestep) {
    return MaxTimestep();
  }

  // About 9 hours which corresponds to 10000 km year-1 on a 10 km grid
  double dt_min = units::convert(m_sys, 0.001, "years", "seconds");

  // Distance (grid cells) from calving front where strain rate is evaluated
  int offset = m_stencil_width;
  double
    my_calving_rate_max     = 0.0,
    my_calving_rate_mean    = 0.0,
    my_calving_rate_counter = 0.0;

  const IceModelVec2CellType &mask = *m_grid->variables().get_2d_cell_type("mask");

  update_strain_rates();

  IceModelVec::AccessList list;
  list.add(mask);
  list.add(m_strain_rates);

  for (Points pt(*m_grid); pt; pt.next()) {
    const int i = pt.i(), j = pt.j();
    // Average of strain-rate eigenvalues in adjacent floating grid cells to
    // be used for eigencalving
    double eigen1 = 0.0, eigen2 = 0.0;
    // Number of cells used in computing eigen1 and eigen2:
    int M = 0;

    // find partially filled or empty grid boxes on the ice-free
    // ocean which have floating ice neighbors
    if ((mask.ice_free_ocean(i, j) and not mask.next_to_floating_ice(i, j))) {
      continue;
    }

    double
      calving_rate_horizontal = 0.0,
      eigenCalvOffset = 0.0;

    for (int p = -1; p < 2; p += 2) {
      int i_offset = p * offset;
      if (mask.floating_ice(i + i_offset, j) and not mask.ice_margin(i + i_offset, j)) {
        eigen1 += m_strain_rates(i + i_offset, j, 0);
        eigen2 += m_strain_rates(i + i_offset, j, 1);
        M += 1;
      }
    }

    for (int q = -1; q < 2; q += 2) {
      int j_offset = q * offset;
      if (mask.floating_ice(i, j + j_offset) and not mask.ice_margin(i,   j + j_offset)) {
        eigen1 += m_strain_rates(i, j + j_offset, 0);
        eigen2 += m_strain_rates(i, j + j_offset, 1);
        M += 1;
      }
    }

    if (M > 0) {
      eigen1 /= M;
      eigen2 /= M;
    }

    // calving law
    if (eigen2 > eigenCalvOffset && eigen1 > 0.0) { // if spreading in all directions
      calving_rate_horizontal = m_K * eigen1 * (eigen2 - eigenCalvOffset);
      my_calving_rate_counter += 1.0;
      my_calving_rate_mean += calving_rate_horizontal;
      my_calving_rate_max = std::max(my_calving_rate_max, calving_rate_horizontal);
    }
  }

  double calving_rate_max = 0.0, calving_rate_mean = 0.0, calving_rate_counter = 0.0;
  calving_rate_mean    = GlobalSum(m_grid->com, my_calving_rate_mean);
  calving_rate_counter = GlobalSum(m_grid->com, my_calving_rate_counter);
  calving_rate_max     = GlobalMax(m_grid->com, my_calving_rate_max);

  if (calving_rate_counter > 0.0) {
    calving_rate_mean /= calving_rate_counter;
  } else {
    calving_rate_mean = 0.0;
  }

  double denom = calving_rate_max / m_grid->dx();
  const double epsilon = units::convert(m_sys, 0.001 / (m_grid->dx() + m_grid->dy()), "seconds", "years");

  double dt = 1.0 / (denom + epsilon);

  m_log->message(3,
             "  eigencalving: max c_rate = %.2f m/a ... gives dt=%.5f a; mean c_rate = %.2f m/a over %d cells\n",
             units::convert(m_sys, calving_rate_max, "m second-1", "m year-1"),
             units::convert(m_sys, dt, "seconds", "years"),
             units::convert(m_sys, calving_rate_mean, "m second-1", "m year-1"),
             (int)calving_rate_counter);

  return MaxTimestep(std::max(dt, dt_min));
}

void EigenCalving::add_vars_to_output_impl(const std::string &/*keyword*/, std::set<std::string> &/*result*/) {
  // empty
}

void EigenCalving::define_variables_impl(const std::set<std::string> &/*vars*/, const PIO &/*nc*/,
                                                  IO_Type /*nctype*/) {
  // empty
}

void EigenCalving::write_variables_impl(const std::set<std::string> &/*vars*/, const PIO& /*nc*/) {
  // empty
}

/**
 * Update the strain rates field.
 *
 * Note: this code uses the mask obtained from the Vars
 * dictionary, because the velocity field used to compute it need not
 * extend past the ice margin corresponding to the *beginning* of the
 * time-step.
 *
 * @return 0 on success
 */
void EigenCalving::update_strain_rates() {
  const IceModelVec2V        &ssa_velocity = m_stress_balance->advective_velocity();
  const IceModelVec2CellType &mask         = *m_grid->variables().get_2d_cell_type("mask");
  m_stress_balance->compute_2D_principal_strain_rates(ssa_velocity, mask, m_strain_rates);
}

/** Remove tips of one-cell-wide ice tongues ("noses").  Changes ice thickness.
 *
 * The center icy cell in ice tongues like this one (and equivalent)
 *
 * @code
   O O ?
   X X O
   O O ?
   @endcode
 *
 * where "O" is ice-free and "?" is any mask value, are removed.
 * Ice tongues like this one
 *
 * @code
   # O ?
   X X O
   # O ?
   @endcode
 * where one or two of the "#" cells are ice-filled, are not removed.
 *
 * See the code for the precise rule, which uses `ice_free_ocean()` for the "O"
 * cells if the center cell has grounded ice, and uses `ice_free()` if the
 * center cell has floating ice.
 *
 * @note We use `pism_mask` (and not ice_thickness) to make decisions.
 * This means that we can update `ice_thickness` in place without
 * introducing a dependence on the grid traversal order.
 *
 * @param[in,out] pism_mask cell type mask
 * @param[in,out] ice_thickness modeled ice thickness
 *
 * @return 0 on success
 */
void EigenCalving::remove_narrow_tongues(IceModelVec2CellType &mask,
                                         IceModelVec2S &ice_thickness) {

  IceModelVec::AccessList list;
  list.add(mask);
  list.add(ice_thickness);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    if (mask.ice_free(i, j)) {
      // FIXME: it might be better to have access to bedrock elevation b(i,j)
      // and sea level SL so that the predicate can be
      //   mask.ice_free(i,j) || (mask.grounded_ice(i,j) && (b(i,j) >= SL)))
      continue;
    }

    bool ice_free_N = false,  ice_free_E = false,
      ice_free_S = false, ice_free_W = false,
      ice_free_NE = false, ice_free_NW = false,
      ice_free_SE = false, ice_free_SW = false;

    if (mask.grounded_ice(i,j)) {
      // if (i,j) is grounded ice then we will remove it if it has
      // exclusively ice-free ocean neighbors
      ice_free_N  = mask.ice_free_ocean(i, j + 1);
      ice_free_E  = mask.ice_free_ocean(i + 1, j);
      ice_free_S  = mask.ice_free_ocean(i, j - 1);
      ice_free_W  = mask.ice_free_ocean(i - 1, j);
      ice_free_NE = mask.ice_free_ocean(i + 1, j + 1);
      ice_free_NW = mask.ice_free_ocean(i - 1, j + 1);
      ice_free_SE = mask.ice_free_ocean(i + 1, j - 1);
      ice_free_SW = mask.ice_free_ocean(i - 1, j - 1);
    } else if (mask.floating_ice(i,j)) {
      // if (i,j) is floating then we will remove it if its neighbors are
      // ice-free, whether ice-free ocean or ice-free ground
      ice_free_N  = mask.ice_free(i, j + 1);
      ice_free_E  = mask.ice_free(i + 1, j);
      ice_free_S  = mask.ice_free(i, j - 1);
      ice_free_W  = mask.ice_free(i - 1, j);
      ice_free_NE = mask.ice_free(i + 1, j + 1);
      ice_free_NW = mask.ice_free(i - 1, j + 1);
      ice_free_SE = mask.ice_free(i + 1, j - 1);
      ice_free_SW = mask.ice_free(i - 1, j - 1);
    }

    if ((not ice_free_W &&
         ice_free_NW         &&
         ice_free_SW         &&
         ice_free_N          &&
         ice_free_S          &&
         ice_free_E)         ||
        (not ice_free_N &&
         ice_free_NW         &&
         ice_free_NE         &&
         ice_free_W          &&
         ice_free_E          &&
         ice_free_S)         ||
        (not ice_free_E &&
         ice_free_NE         &&
         ice_free_SE         &&
         ice_free_W          &&
         ice_free_S          &&
         ice_free_N)         ||
        (not ice_free_S &&
         ice_free_SW         &&
         ice_free_SE         &&
         ice_free_W          &&
         ice_free_E          &&
         ice_free_N)) {
      mask(i, j)          = MASK_ICE_FREE_OCEAN;
      ice_thickness(i, j) = 0.0;
    }
  }
}

} // end of namespace calving
} // end of namespace pism
