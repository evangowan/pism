// Copyright (C) 2012--2016 Ricarda Winkelmann, Ronja Reese and Torsten Albrecht
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
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

#ifndef _POOCEANBOXMODEL_H_
#define _POOCEANBOXMODEL_H_

#include "PGivenClimate.hh"
#include "POModifier.hh"
#include "Timeseries.hh"

namespace pism {
namespace ocean {

//! A class defining the interface of a PISM ocean model modifier.

//! \brief A class implementing an ocean model.
//! Computes the subshelf melt/refreezing rate based on a simple ocean box model
//by Olbers & Hellmer (2010).

class BoxModel : public PGivenClimate<POModifier, PISMOceanModel> {
public:
  BoxModel(IceGrid::ConstPtr g);
  virtual ~BoxModel();

  virtual void add_vars_to_output(std::string keyword,
                                  std::set<std::string> &result);

  virtual void define_variables(std::set<std::string> vars, const PIO &nc,
                                PISM_IO_Type nctype);

  virtual void write_variables(std::set<std::string> vars, const PIO &nc);

  class POBMConstants {
  public:
    POBMConstants(const PISMConfig &config);

    double earth_grav, rhoi, rhow, rho_star, nu, latentHeat, c_p_ocean, lambda,
        a, b, c, alpha, beta;

    double gamma_T, value_C, T_dummy, S_dummy;

    double gamma_T_o, meltFactor, meltSalinity, b2;
    double continental_shelf_depth;

    int numberOfBasins;
  };

protected:
  virtual void shelf_base_mass_flux_impl(IceModelVec2S &result);
  virtual void melange_back_pressure_fraction_impl(IceModelVec2S &result);
  virtual void shelf_base_temperature_impl(IceModelVec2S &result);
  virtual void sea_level_elevation_impl(double &result);
  virtual void update_impl(double t, double dt);
  virtual void init_impl();

private:
  IceModelVec2S shelfbtemp, shelfbmassflux;

  IceModelVec2T *theta_ocean, *salinity_ocean;

  IceModelVec2S *ice_thickness, *topg, *basins; // not owned by this class

  IceModelVec2Int *mask; // not owned by this class

  virtual void initBasinsOptions(const POBMConstants &constants);
  virtual void roundBasins();
  virtual void identifyMASK(IceModelVec2S &inputmask, std::string masktype);
  virtual void computeOCEANMEANS();
  virtual void extentOfIceShelves();
  virtual void identifyBOXMODELmask();
  virtual void extendGLBox();
  virtual void extendIFBox();
  virtual void oceanTemperature(const POBMConstants &constants);
  virtual void basalMeltRateForGroundingLineBox(const POBMConstants &constants);
  virtual void basalMeltRateForIceFrontBox(const POBMConstants &constants);
  virtual void basalMeltRateForOtherShelves(const POBMConstants &constants);

  static const int box_unidentified, box_noshelf, box_GL, box_neighboring,
      box_IF, box_other,
      maskfloating, maskocean, maskgrounded,
      imask_inner, imask_outer, imask_exclude, imask_unidentified;

  double counter_box_unidentified;

  std::vector<double> Toc_base_vec, Soc_base_vec, gamma_T_star_vec, C_vec,

      mean_salinity_GLbox_vector, mean_meltrate_GLbox_vector,
      mean_overturning_GLbox_vector,

      counter, counter_GLbox, counter_CFbox;

  IceModelVec2S ICERISESmask, BOXMODELmask,
      OCEANMEANmask, // FIXME delete OCEANMEANmask
      CHECKmask, // FIXME delete CHECKmask
      Soc, Soc_base, Toc, Toc_base, Toc_inCelsius, T_star, Toc_anomaly,
      overturning, heatflux, basalmeltrate_shelf;

  double gamma_T, value_C, T_dummy, S_dummy, continental_shelf_depth;

  int numberOfBasins;

protected:
  std::vector<IceModelVec*> m_variables;
  Timeseries *delta_T;
  double delta_T_factor;
  double temp_anomaly;

  bool ocean_oceanboxmodel_deltaT_set, exicerises_set,
      continental_shelf_depth_set;
};

} // end of namespace ocean
} // end of namespace pism

#endif /* _POOCEANBOXMODEL_H_ */
