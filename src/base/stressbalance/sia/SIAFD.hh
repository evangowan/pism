// Copyright (C) 2004--2016 Jed Brown, Ed Bueler and Constantine Khroulev
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

#ifndef _SIAFD_H_
#define _SIAFD_H_

#include "base/stressbalance/SSB_Modifier.hh"      // derivesfrom SSB_Modifier

namespace pism {

namespace stressbalance {

class BedSmoother;

/** Implements the shallow ice approximation stress balance.
 *
 * Inputs:
 *
 * - ice geometry (thickness, bed elevation, surface elevation, cell
 *   type mask)
 * - ice enthalpy
 * - ice age (could be used to compute the grain size)
 * - sliding velocity
 *
 * Outputs:
 *
 * - horizontal velocity (3D fields)
 * - diffusive ice flux (for use in the geometry update)
 * - maximum diffusivity (used to determine the maximum allowed time
 *   step length)
 * - volumetric strain heating
 */
class SIAFD : public SSB_Modifier
{
  friend class SIAFD_schoofs_theta;
  friend class SIAFD_topgsmooth;
  friend class SIAFD_thksmooth;
  friend class SIAFD_diffusivity;
  friend class SIAFD_diffusivity_staggered;
  friend class SIAFD_h_x;
  friend class SIAFD_h_y;
public:
  SIAFD(IceGrid::ConstPtr g, EnthalpyConverter::Ptr e);

  virtual ~SIAFD();

  virtual void init();

  virtual void update(const IceModelVec2V &vel_input, bool fast);

protected:
  virtual void get_diagnostics_impl(std::map<std::string, Diagnostic::Ptr> &dict,
                                    std::map<std::string, TSDiagnostic::Ptr> &ts_dict);
  virtual void write_variables_impl(const std::set<std::string> &vars, const PIO &nc);
  virtual void add_vars_to_output_impl(const std::string &keyword,
                                  std::set<std::string> &result);
  virtual void define_variables_impl(const std::set<std::string> &vars, const PIO &nc,
                                     IO_Type nctype);

  virtual void compute_surface_gradient(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y);

  virtual void surface_gradient_eta(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y);
  virtual void surface_gradient_haseloff(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y);
  virtual void surface_gradient_mahaffy(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y);

  virtual void compute_diffusive_flux(const IceModelVec2Stag &h_x, const IceModelVec2Stag &h_y,
                                      IceModelVec2Stag &result, bool fast);

  virtual void compute_3d_horizontal_velocity(const IceModelVec2Stag &h_x,
                                              const IceModelVec2Stag &h_y,
                                              const IceModelVec2V &vel_input,
                                              IceModelVec3 &u_out, IceModelVec3 &v_out);

  virtual void compute_I();

  virtual double grainSizeVostok(double age) const;

  virtual void compute_diffusivity(IceModelVec2S &result);
  virtual void compute_diffusivity_staggered(IceModelVec2Stag &result);

  bool interglacial(double accumulation_time);

  //! temporary storage for eta, theta and the smoothed thickness
  IceModelVec2S m_work_2d[2];
  //! temporary storage for the surface gradient
  IceModelVec2Stag m_work_2d_stag[2];
  //! temporary storage for delta on the staggered grid
  IceModelVec3 m_delta[2];
  //! temporary storage used to store I and strain_heating on the staggered grid
  IceModelVec3 m_work_3d[2];

  BedSmoother *m_bed_smoother;
  int m_bed_state_counter;

  // profiling
  int m_event_sia;

  // unit conversion
  double m_second_to_kiloyear;
  // enhancement factor-age coupling parameters
  double m_holocene_start;
  double m_eemian_start;
  double m_eemian_end;
};

} // end of namespace stressbalance
} // end of namespace pism

#endif /* _SIAFD_H_ */
