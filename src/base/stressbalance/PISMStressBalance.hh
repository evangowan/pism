// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2016 Constantine Khroulev and Ed Bueler
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

#ifndef _PISMSTRESSBALANCE_H_
#define _PISMSTRESSBALANCE_H_

#include "base/util/PISMComponent.hh"     // derives from Component
#include "base/util/iceModelVec.hh"
#include "base/timestepping.hh"

namespace pism {

class IceModelVec2CellType;

//! Stress balance models and related diagnostics.
namespace stressbalance {

class StressBalanceInputs {
public:
  StressBalanceInputs();
  double sea_level;

  // geometry
  const IceModelVec2S        *ice_thickness;
  const IceModelVec2S        *bed_elevation;
  const IceModelVec2S        *surface_elevation;
  const IceModelVec2CellType *cell_type;
  const IceModelVec2S        *grounded_cell_fraction;

  // energy
  const IceModelVec3 *ice_enthalpy;

  // velocity
  const IceModelVec2S *basal_melt_rate;

  // boundary conditions
  const IceModelVec2S *melange_back_pressure;
  const IceModelVec2S *basal_yield_stress;

  const IceModelVec2S *fracture_density;

  // Support for direct specification of driving stress to the FEM SSA solver. This helps with
  // certain test cases where the grid is periodic but the driving stress cannot be the gradient of
  // a periodic function. (See commit ffb4be16.)
  const IceModelVec2S *driving_stress_x;
  const IceModelVec2S *driving_stress_y;
};


class ShallowStressBalance;
class SSB_Modifier;

//! The class defining PISM's interface to the shallow stress balance code.
/*!
  Generally all the nontrivial fields are updated by a call to update().  The rest
  of the methods generally provide access to precomputed results.  The following
  diagram shows where these results are generally used in the rest of PISM.  (It 
  does not show the call graph, as would doxygen.)

  \image html stressbalance-out.png "\b Methods of StressBalance, and the uses of their results.  Dotted edges show scalars and dashed edges show fields.  Dashed boxes inside the StressBalance object are important methods which may be present in shallow cases.  The age time step has inputs which are a strict subset of the inputs of the energy time step."

  this command fails: \dotfile stressbalance-out.dot
*/
class StressBalance : public Component
{
public:
  StressBalance(IceGrid::ConstPtr g, ShallowStressBalance *sb, SSB_Modifier *ssb_mod);
  virtual ~StressBalance();

  //! \brief Initialize the StressBalance object.
  void init();

  //! \brief Update all the fields if (not fast), only update diffusive flux
  //! and max. diffusivity otherwise.
  void update(bool fast, const StressBalanceInputs &inputs);

  //! \brief Get the thickness-advective (SSA) 2D velocity.
  const IceModelVec2V& advective_velocity() const;

  //! \brief Get the diffusive (SIA) vertically-averaged flux on the staggered grid.
  const IceModelVec2Stag& diffusive_flux() const;

  //! \brief Get the max diffusivity (for the adaptive time-stepping).
  double max_diffusivity() const;

  CFLData max_timestep_cfl_2d() const;
  CFLData max_timestep_cfl_3d() const;

  // for the energy/age time step:

  //! \brief Get components of the the 3D velocity field.
  const IceModelVec3& velocity_u() const;
  const IceModelVec3& velocity_v() const;
  const IceModelVec3& velocity_w() const;

  //! \brief Get the basal frictional heating.
  const IceModelVec2S& basal_frictional_heating() const;

  const IceModelVec3& volumetric_strain_heating() const;

  // for the calving, etc.:

  //! \brief Get the components of the 2D deviatoric stress tensor.
  void compute_2D_stresses(const IceModelVec2V &velocity,
                           const IceModelVec2CellType &mask,
                           IceModelVec2 &result) const;

  //! \brief Produce a report string for the standard output.
  std::string stdout_report() const;

  //! \brief Returns a pointer to a shallow stress balance solver implementation.
  const ShallowStressBalance* shallow() const;

  //! \brief Returns a pointer to a stress balance modifier implementation.
  const SSB_Modifier* modifier() const;
protected:
  virtual void get_diagnostics_impl(std::map<std::string, Diagnostic::Ptr> &dict,
                                    std::map<std::string, TSDiagnostic::Ptr> &ts_dict) const;

  virtual void define_model_state_impl(const PIO &output) const;
  virtual void write_model_state_impl(const PIO &output) const;

  virtual void compute_vertical_velocity(const IceModelVec3 &u,
                                         const IceModelVec3 &v,
                                         const IceModelVec2CellType &mask,
                                         const IceModelVec2S *bmr,
                                         IceModelVec3 &result);
  virtual void compute_volumetric_strain_heating(const IceModelVec2S &thickness,
                                                 const IceModelVec3 &enthalpy,
                                                 const IceModelVec2CellType &mask);

  CFLData compute_cfl_2d(const IceModelVec2S &ice_thickness,
                         const IceModelVec2CellType &cell_type);
  CFLData compute_cfl_3d(const IceModelVec2S &ice_thickness,
                         const IceModelVec2CellType &cell_type);

  CFLData m_cfl_2d, m_cfl_3d;

  IceModelVec3 m_w, m_strain_heating;

  ShallowStressBalance *m_shallow_stress_balance;
  SSB_Modifier *m_modifier;
};

void compute_2D_principal_strain_rates(const IceModelVec2V &velocity,
                                       const IceModelVec2CellType &mask,
                                       IceModelVec2 &result);

} // end of namespace stressbalance
} // end of namespace pism

#endif /* _PISMSTRESSBALANCE_H_ */

