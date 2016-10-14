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

#ifndef _SHALLOWSTRESSBALANCE_H_
#define _SHALLOWSTRESSBALANCE_H_

#include "base/util/PISMComponent.hh"
#include "base/util/iceModelVec.hh"
#include "base/enthalpyConverter.hh"

namespace pism {
namespace rheology {
class FlowLaw;
}

class IceGrid;
class IceBasalResistancePlasticLaw;
class IceModelVec2CellType;

namespace stressbalance {

class ShallowStressBalanceInputs {
public:
  ShallowStressBalanceInputs();
  double sea_level;

  const IceModelVec2V *bc_values;
  const IceModelVec2Int *bc_mask;

  const IceModelVec2S *melange_back_pressure;
  const IceModelVec2S *ice_thickness;
  const IceModelVec2S *bed_elevation;
  const IceModelVec2S *surface_elevation;
  const IceModelVec3  *ice_enthalpy;

  // Support for direct specification of driving stress to the FEM SSA solver. This helps
  // with certain test cases where the grid is periodic but the driving stress cannot be the
  // gradient of a periodic function. (See commit ffb4be16.)
  const IceModelVec2S *driving_stress_x;
  const IceModelVec2S *driving_stress_y;
};

//! Shallow stress balance (such as the SSA).
class ShallowStressBalance : public Component {
public:
  ShallowStressBalance(IceGrid::ConstPtr g);
  virtual ~ShallowStressBalance();

  //  initialization and I/O:

  void init();
  void set_boundary_conditions(const IceModelVec2Int &locations,
                               const IceModelVec2V &velocities);

  virtual void update(bool fast, const ShallowStressBalanceInputs &inputs) = 0;

  //! \brief Get the thickness-advective 2D velocity.
  const IceModelVec2V& velocity() const;

  //! \brief Get the basal frictional heating (for the adaptive energy time-stepping).
  const IceModelVec2S& basal_frictional_heating();

  void compute_2D_stresses(const IceModelVec2V &velocity,
                           const IceModelVec2CellType &mask,
                           IceModelVec2 &result) const;

  void compute_basal_frictional_heating(const IceModelVec2V &velocity,
                                        const IceModelVec2S &tauc,
                                        const IceModelVec2CellType &mask,
                                        IceModelVec2S &result) const;
  // helpers:

  //! \brief Produce a report string for the standard output.
  virtual std::string stdout_report() const;

  const rheology::FlowLaw* flow_law() const;

  EnthalpyConverter::Ptr enthalpy_converter() const;

  const IceBasalResistancePlasticLaw* sliding_law() const;
protected:
  virtual void init_impl();
  
  virtual void get_diagnostics_impl(std::map<std::string, Diagnostic::Ptr> &dict,
                                    std::map<std::string, TSDiagnostic::Ptr> &ts_dict) const;

  double m_sea_level;
  IceBasalResistancePlasticLaw *m_basal_sliding_law;
  rheology::FlowLaw *m_flow_law;
  EnthalpyConverter::Ptr m_EC;

  IceModelVec2V m_velocity;
  const IceModelVec2V *m_bc_values;
  const IceModelVec2Int *m_bc_mask;
  IceModelVec2S m_basal_frictional_heating;
};

//! Returns zero velocity field, zero friction heating, and zero for D^2.
/*!
  This derived class is used in the non-sliding SIA approximation. This
  implementation ignores any basal resistance fields (e.g. yield stress from
  the IceModel or other user of this class).
*/
class ZeroSliding : public ShallowStressBalance {
public:
  ZeroSliding(IceGrid::ConstPtr g);
  virtual ~ZeroSliding();
  
  virtual void update(bool fast, const ShallowStressBalanceInputs &inputs);

protected:
};

class PrescribedSliding : public ZeroSliding {
public:
  PrescribedSliding(IceGrid::ConstPtr g);
  virtual ~PrescribedSliding();
  virtual void update(bool fast, const ShallowStressBalanceInputs &inputs);
  virtual void init();
};

} // end of namespace stressbalance
} // end of namespace pism

#endif /* _SHALLOWSTRESSBALANCE_H_ */
