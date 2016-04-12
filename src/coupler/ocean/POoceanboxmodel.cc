// Copyright (C) 2012-2016 Ricarda Winkelmann, Ronja Reese, Torsten Albrecht and Constantine Khroulev
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

#include "POoceanboxmodel.hh"
#include "base/util/IceGrid.hh"
#include "base/util/PISMConfig.hh"
#include "base/util/PISMVars.hh"

namespace pism {
namespace ocean {

/*!
  The aim of these routines is to simulate the ocean circulation underneath the
  ice shelves and compute the basal melt/refreezing rates according to the ocean
  box model described in olbers_hellmer10.
*/

BoxModel::POBMConstants::POBMConstants(const PISMConfig &config) {

  numberOfBasins = 20;

  continental_shelf_depth = -800;

  T_dummy = -1.5; // FIXME why these values?
  S_dummy = 34.5; //

  earth_grav = config.get("standard_gravity");
  rhoi       = config.get("ice_density");
  rhow       = config.get("sea_water_density");
  rho_star   = 1033;            // kg/m^3
  nu         = rhoi / rho_star; // no unit

  latentHeat = config.get("water_latent_heat_fusion");
  c_p_ocean  = 3974.0;                 // J/(Kelvin*kg), specific heat capacity of ocean mixed layer
  lambda     = latentHeat / c_p_ocean; // Celsius, NOTE Kelvin vs Celsius

  a = -0.057;  // Celsius/psu
  b = 0.0832;  // Celsius
  c = 7.64e-4; // Celsius/dbar

  alpha = 7.5e-5; // 1/Celsius, NOTE Kelvin vs Celsius
  beta  = 7.7e-4; // 1/psu

  gamma_T = 1e-6;
  value_C = 5e6;

  // other ice shelves
  gamma_T_o    = 1.0e-4;
  meltFactor   = 0.002; // FIXME!!!! (Wrong value!) FIXME config
  meltSalinity = 35.0;
  b2           = 0.0939;
}

const int BoxModel::box_unidentified = -99; // This should never show up in the .nc-files.
const int BoxModel::box_neighboring  = -1;  // This should never show up in the .nc-files.
const int BoxModel::box_noshelf      = 0;
const int BoxModel::box_GL           = 1; // ocean box covering the grounding line region
const int BoxModel::box_IF           = 2; // ocean box covering the rest of the ice shelf
const int BoxModel::box_other        = 3; // ice_shelf but there is no GL_box in the corresponding basin

const int BoxModel::maskfloating = MASK_FLOATING;
const int BoxModel::maskocean    = MASK_ICE_FREE_OCEAN;
const int BoxModel::maskgrounded = MASK_GROUNDED;

const int BoxModel::imask_inner        = 2;
const int BoxModel::imask_outer        = 0;
const int BoxModel::imask_exclude      = 1;
const int BoxModel::imask_unidentified = -1;

BoxModel::BoxModel(IceGrid &g, const PISMConfig &conf) : PGivenClimate<POModifier, OceanModel>(g, conf, NULL) {
  option_prefix = "-ocean_oceanboxmodel";

  // will be de-allocated by the parent's destructor
  theta_ocean    = new IceModelVec2T;
  salinity_ocean = new IceModelVec2T;

  m_fields["theta_ocean"]    = theta_ocean;
  m_fields["salinity_ocean"] = salinity_ocean;

  process_options();

  std::map<std::string, std::string> standard_names;
  set_vec_parameters(standard_names);

  theta_ocean->create(grid, "theta_ocean");
  theta_ocean->set_attrs("climate_forcing", "absolute potential temperature of the adjacent ocean", "Kelvin", "");

  salinity_ocean->create(grid, "salinity_ocean");
  salinity_ocean->set_attrs("climate_forcing", "salinity of the adjacent ocean", "g/kg", "");

  shelfbtemp.create(grid, "shelfbtemp", WITHOUT_GHOSTS);
  shelfbtemp.set_attrs("climate_forcing", "absolute temperature at ice shelf base", "Kelvin", "");
  m_variables.push_back(&shelfbtemp);

  shelfbmassflux.create(grid, "shelfbmassflux", WITHOUT_GHOSTS);
  shelfbmassflux.set_attrs(
      "climate_forcing", "ice mass flux from ice shelf base (positive flux is loss from ice shelf)", "kg m-2 s-1", "");
  shelfbmassflux.set_glaciological_units("kg m-2 year-1");
  shelfbmassflux.write_in_glaciological_units = true;
  m_variables.push_back(&shelfbmassflux);
  //////////////////////////////////////////////////////////////////////////

  // mask to identify the ocean boxes
  BOXMODELmask.create(grid, "BOXMODELmask", WITH_GHOSTS);
  BOXMODELmask.set_attrs("model_state", "mask displaying ocean box model grid", "", "");
  m_variables.push_back(&BOXMODELmask);

  // mask to identify the grounded ice rises
  ICERISESmask.create(grid, "ICERISESmask", WITH_GHOSTS);
  ICERISESmask.set_attrs("model_state", "mask displaying ice rises", "", "");
  m_variables.push_back(&ICERISESmask);

  exicerises_set = options::Bool("-exclude_icerises", "FIXME: add description");

  // mask displaying continental shelf - region where mean salinity and ocean
  // temperature is calculated
  OCEANMEANmask.create(grid, "OCEANMEANmask", WITH_GHOSTS);
  OCEANMEANmask.set_attrs("model_state", "mask displaying ocean region for parameter input", "", "");
  m_variables.push_back(&OCEANMEANmask);

  // salinity
  Soc.create(grid, "Soc", WITHOUT_GHOSTS);
  Soc.set_attrs("model_state", "ocean salinity field", "", "ocean salinity field"); // NOTE unit=psu
  m_variables.push_back(&Soc);

  Soc_base.create(grid, "Soc_base", WITHOUT_GHOSTS);
  Soc_base.set_attrs("model_state", "ocean base salinity field", "", "ocean base salinity field"); // NOTE unit=psu
  m_variables.push_back(&Soc_base);

  // temperature
  Toc.create(grid, "Toc", WITHOUT_GHOSTS);
  Toc.set_attrs("model_state", "ocean temperature field", "Kelvin", "ocean temperature field");
  m_variables.push_back(&Toc);

  Toc_base.create(grid, "Toc_base", WITHOUT_GHOSTS);
  Toc_base.set_attrs("model_state", "ocean base temperature", "Kelvin", "ocean base temperature");
  m_variables.push_back(&Toc_base);

  Toc_inCelsius.create(grid, "Toc_inCelsius", WITHOUT_GHOSTS);
  Toc_inCelsius.set_attrs("model_state", "ocean box model temperature field", "degree C",
                          "ocean box model temperature field");
  m_variables.push_back(&Toc_inCelsius);

  T_star.create(grid, "T_star", WITHOUT_GHOSTS);
  T_star.set_attrs("model_state", "T_star field", "degree C", "T_star field");
  m_variables.push_back(&T_star);

  Toc_anomaly.create(grid, "Toc_anomaly", WITHOUT_GHOSTS);
  Toc_anomaly.set_attrs("model_state", "ocean temperature anomaly", "Kelvin", "ocean temperature anomaly");
  m_variables.push_back(&Toc_anomaly);

  // overturning rate
  overturning.create(grid, "overturning", WITHOUT_GHOSTS);
  overturning.set_attrs("model_state", "cavity overturning", "m^3 s-1", "cavity overturning");
  m_variables.push_back(&overturning);

  // heat flux
  heatflux.create(grid, "ocean heat flux", WITHOUT_GHOSTS);
  heatflux.set_attrs("climate_state", "ocean heat flux", "W/m^2", "");
  m_variables.push_back(&heatflux);

  // basal melt rate
  basalmeltrate_shelf.create(grid, "basal melt rate from ocean box model", WITHOUT_GHOSTS);
  basalmeltrate_shelf.set_attrs("climate_state", "basal melt rate from ocean box model", "m/s", "");
  basalmeltrate_shelf.set_glaciological_units("m year-1");
  basalmeltrate_shelf.write_in_glaciological_units = true;
  m_variables.push_back(&basalmeltrate_shelf);

  ///////// forcing
  ////////////////////////////////////////////////////////////////////////////////////

  // option for scalar forcing of ocean temperature
  ocean_oceanboxmodel_deltaT_set = options::Bool("-ocean_obm_deltaT", "FIXME: add description");

  if (ocean_oceanboxmodel_deltaT_set) {
    std::string delta_T_file;

    options::String delta_T_file("-ocean_obm_deltaT",
                                 "Specifies the ocean temperature offsets file to use with -ocean_obm_deltaT");

    m_log->message(2, "  reading delta_T data from forcing file %s for -ocean_obm_deltaT actions ...\n",
                   delta_T_file.c_str());

    delta_T = new Timeseries(&grid, "delta_T", grid.config.get_string("time_dimension_name"));
    delta_T->set_units("Kelvin", "");
    delta_T->set_dimension_units(grid.time->units_string(), "");
    delta_T->set_attr("long_name", "ocean temperature offsets");
    // delta_T->read(delta_T_file, grid.time->use_reference_date());

    {
      PIO nc(grid.com, "netcdf3", grid.get_unit_system());
      nc.open(delta_T_file, PISM_NOWRITE);
      delta_T->read(nc, grid.time);
      nc.close();
    }

    delta_T_factor = options::Real("-ocean_obm_factor", "ocean_obm_factor set", 1.0);
  }
}

BoxModel::~BoxModel() {
  // empty
}

void BoxModel::init_impl() {

  m_log->message(2, "* Initializing the ocean box model (based on Olbers & Hellmer (2010)...\n");

  m_t = m_dt = GSL_NAN; // every re-init restarts the clock

  ice_thickness = vars.get_2d_scalar("land_ice_thickness");

  topg = vars.get_2d_scalar("bedrock_altitude");

  mask = vars.get_2d_mask("mask");

  basins = vars.get_2d_scalar("drainage_basins"); // if option drainageBasins
  // set

  bool omeans_set = options::Bool("-ocean_means", "read mean salinity and temperatures")

                    // FIXME: not necessary when -ocean_means set
                    theta_ocean->init(filename, bc_period, bc_reference_time);
  salinity_ocean->init(filename, bc_period, bc_reference_time);

  // read time-independent data right away:
  if (theta_ocean->get_n_records() == 1 and salinity_ocean->get_n_records() == 1) {
    update(grid.time->current(), 0); // dt is irrelevant
  }

  POBMConstants cc(config);
  initBasinsOptions(cc);
}

void BoxModel::add_vars_to_output(std::string keyword, std::set<std::string> &result) {
  PGivenClimate<POModifier, OceanModel>::add_vars_to_output(keyword, result);

  if (keyword != "none" and keyword != "small") {
    result.insert("shelfbtemp");
    result.insert("shelfbmassflux");
  }
}

void BoxModel::initBasinsOptions(const POBMConstants &cc) {
  m_log->message(4, "0b : set number of Basins\n");

  // set number of basins per option
  numberOfBasins = options::Integer("-number_of_basins", "Number of Drainage Basins", cc.numberOfBasins);

  Toc_base_vec.resize(numberOfBasins);
  Soc_base_vec.resize(numberOfBasins);
  gamma_T_star_vec.resize(numberOfBasins);
  C_vec.resize(numberOfBasins);

  counter.resize(numberOfBasins);
  counter_GLbox.resize(numberOfBasins);
  counter_CFbox.resize(numberOfBasins);
  // k_basins.resize(numberOfBasins);;

  mean_salinity_GLbox_vector.resize(numberOfBasins);
  mean_meltrate_GLbox_vector.resize(numberOfBasins);
  mean_overturning_GLbox_vector.resize(numberOfBasins);

  // set gamma_T and value_C per option
  gamma_T = options::Real("-gamma_T", "FIXME: add description", cc.gamma_T);

  value_C = options::Real("-value_C", "FIXME: add description", cc.value_C);

  ///////////////////////////////////////////////////////////////////////////////////
  // data have been calculated previously for the 18 Rignot basins

  // const double Toc_base_schmidtko[18] = {0.0, 271.56415203, 271.63356482,
  // 271.42074233, 271.46720524, 272.30038929, 271.52821139, 271.5440751,
  // 271.58201494, 272.90159695, 273.61058862, 274.19203524, 274.32083917,
  // 272.55938554, 271.35349906, 271.39337366, 271.49926019, 271.49473924};
  // //Schmidtko
  // const double Soc_base_schmidtko[18] = {0.0, 34.49308909, 34.50472554,
  // 34.70187911, 34.65306507, 34.7137078, 34.74817136, 34.89206844,
  // 34.78056731, 34.60528314, 34.72521443, 34.86210624, 34.83836297,
  // 34.73392011, 34.83617088, 34.82137147, 34.69477334, 34.48145265};
  // //Schmidtko

  // data have been calculated previously for the 20 Zwally basins
  const double Toc_base_schmidtko[20] = {
    0.0,          271.39431005, 271.49081157, 271.49922596, 271.56714804, 271.63507013, 271.42228667,
    271.46720524, 272.42253843, 271.53779093, 271.84942002, 271.31676801, 271.56846696, 272.79372542,
    273.61694268, 274.19168456, 274.31958227, 273.38372579, 271.91951514, 271.35349906
  }; // Schmidtko

  const double Soc_base_schmidtko[20] = {
    0.0,         34.82193374, 34.69721226, 34.47641407, 34.48950162, 34.50258917, 34.70101507,
    34.65306507, 34.73295029, 34.74859586, 34.8368573,  34.9529016,  34.79486795, 34.58380953,
    34.7260615,  34.86198383, 34.8374212,  34.70418016, 34.75598208, 34.83617088
  }; // Schmidtko

  // const double Toc_base_woa[18] = {0.0, 272.28351693, 272.10101401,
  // 271.65965597, 271.50766979, 273.02732277, 272.12473624, 271.79505722,
  // 271.93548261, 273.37866926, 272.98126135, 273.73564726, 273.95971315,
  // 273.02383769, 272.56732024, 271.75152607, 271.93962932, 272.46601985};
  // //World Ocean Atlas
  // const double Soc_base_woa[18] = {0.0, 34.48230812, 34.46499742,
  // 34.51939747, 34.40695979, 34.62947893, 34.59932424, 34.77118004,
  // 34.71666183, 34.54603987, 34.44824601, 34.62416923, 34.57034648,
  // 34.60459029, 34.68356516, 34.67159838, 34.62308218, 34.49961882}; //World
  // Ocean Atlas

  const double Toc_base_woa[20] = {
    272.99816667, 271.27814004, 272.1840257,  272.04435251, 272.20415662, 272.36396072, 271.48763831,
    271.99695864, 272.06504052, 272.27114732, 272.66657018, 271.18920729, 271.74067699, 273.01811291,
    272.15295572, 273.08542047, 272.74584469, 273.14263356, 272.58496563, 272.45217911
  }; // World Ocean Atlas
  const double Soc_base_woa[20] = {
    34.6810522,  34.78161073, 34.67151084, 34.66538478, 34.67127468, 34.67716458, 34.75327377,
    34.69213327, 34.72086382, 34.70670158, 34.71210592, 34.80229468, 34.76588022, 34.69745763,
    34.7090778,  34.68690903, 34.66379606, 34.64572337, 34.6574402,  34.65813983
  }; // World Ocean Atlas

  options::String ocean_means("-ocean_means", "Input data name");

  /////////////////////////////////////////////////////////////////////////////////////

  for (int k = 0; k < numberOfBasins; k++) {
    if (ocean_means) {
      if (ocean_means == "schmidtko") {
        Toc_base_vec[k] = Toc_base_schmidtko[k] - 273.15;
        Soc_base_vec[k] = Soc_base_schmidtko[k];
      } else if (ocean_means == "woa") {
        Toc_base_vec[k] = Toc_base_woa[k] - 273.15;
        Soc_base_vec[k] = Soc_base_woa[k];
      } else {
        Toc_base_vec[k] = cc.T_dummy; // dummy, FIXME why these values?
        Soc_base_vec[k] = cc.S_dummy; // dummy
      }
    } else {
      Toc_base_vec[k] = cc.T_dummy; // dummy, FIXME why these values?
      Soc_base_vec[k] = cc.S_dummy; // dummy
    }

    gamma_T_star_vec[k] = gamma_T;
    C_vec[k]            = value_C;
  }

  m_log->message(5, "     Using %d drainage basins and default values: \n"
                    "     gamma_T_star= %.2e, C = %.2e... \n",
                 numberOfBasins, gamma_T, value_C);

  if (not ocean_means) {
    m_log->message(5, "  calculate Soc and Toc from thetao and salinity... \n");

    // set continental shelf depth // FIXME -800 might be too high for
    // Antarctica
    continental_shelf_depth = cc.continental_shelf_depth;
    options::Real shelf_depth("-continental_shelf_depth", "-continental_shelf_depth", continental_shelf_depth);
    if (shelf_depth) {
      continental_shelf_depth = shelf_depth;
      m_log->message(5, "  Depth of continental shelf for computation of temperature and salinity input\n"
                        "  is set for whole domain to continental_shelf_depth=%.0f meter\n",
                     shelf_depth);
    }
  }
}

void BoxModel::update(double my_t, double my_dt) {

  // Make sure that sea water salinity and sea water potential
  // temperature fields are up to date:
  update_internal(my_t, my_dt);

  theta_ocean->average(m_t, m_dt);
  salinity_ocean->average(m_t, m_dt);

  if ((delta_T != NULL) and ocean_oceanboxmodel_deltaT_set) {
    temp_anomaly = (*delta_T)(my_t + 0.5 * my_dt);
    m_log->message(4, "0a : set global temperature anomaly = %.3f\n", temp_anomaly);
  }

  bool omeans_set = options::Bool("-ocean_means", "read mean salinity and temperatures");

  POBMConstants cc(config);

  initBasinsOptions(cc);

  roundBasins();
  if (omeans_set) {
    m_log->message(4, "0c : reading mean salinity and temperatures\n");
  } else {
    m_log->message(4, "0c : calculating mean salinity and temperatures\n");
    identifyMASK(OCEANMEANmask, "ocean");
    computeOCEANMEANS();
  }

  // geometry of ice shelves and temperatures
  m_log->message(4, "A  : calculating shelf_base_temperature\n");
  if (exicerises_set) {
    identifyMASK(ICERISESmask, "icerises");
  }
  extentOfIceShelves();
  identifyBOXMODELmask();
  oceanTemperature(cc);
  Toc.copy_to(shelfbtemp);

  // basal melt rates underneath ice shelves
  m_log->message(4, "B  : calculating shelf_base_mass_flux\n");
  basalMeltRateForGroundingLineBox(cc);
  basalMeltRateForIceFrontBox(cc);
  // TODO Diese Routinen woanders aufrufen (um Dopplung zu vermeiden)
  basalMeltRateForOtherShelves(cc); // Assumes that mass flux is proportional to
                                    // the shelf-base heat flux.
  // const double secpera=31556926.0;
  basalmeltrate_shelf.scale(cc.rhoi);
  basalmeltrate_shelf.copy_to(shelfbmassflux);
  // TODO Check if scaling with ice density
}

//! Round basin mask non integer values to an integral value of the next neighbor
void BoxModel::roundBasins() {
  // FIXME: THIS routine should be applied once in init, and roundbasins should
  // be stored as a field

  bool round_basins = options::Bool("-round_basins", "FIXME: add description");

  list.add(*basins);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double
      id_fractional = (*basins)(i, j),
      id_fr_ne      = (*basins)(i + 1, j + 1),
      id_fr_nw      = (*basins)(i - 1, j + 1),
      id_fr_sw      = (*basins)(i - 1, j - 1),
      id_fr_se      = (*basins)(i + 1, j - 1);

    int
      id_rounded = static_cast<int>(round(id_fractional)),
      id_ro_ne   = static_cast<int>(round(id_fr_ne)),
      id_ro_nw   = static_cast<int>(round(id_fr_nw)),
      id_ro_sw   = static_cast<int>(round(id_fr_sw)),
      id_ro_se   = static_cast<int>(round(id_fr_se)),
      id         = -1;

    if (round_basins) {

      if (PetscAbs(id_fractional - static_cast<float>(id_rounded)) > 0.0) {
        // if id_fractional differs from integer value

        if (id_fr_sw == static_cast<float>(id_ro_sw) and id_fr_sw != 0) {
          id = id_ro_sw;
        } else if (id_fr_se == static_cast<float>(id_ro_se) and id_fr_se != 0) {
          id = id_ro_se;
        } else if (id_fr_nw == static_cast<float>(id_ro_nw) and id_fr_nw != 0) {
          id = id_ro_nw;
        } else if (id_fr_ne == static_cast<float>(id_ro_ne) and id_fr_ne != 0) {
          id = id_ro_ne;
        } else {
          // if no neighbor has an integer id
          id = id_rounded;
        }
      } else {
        // if id_rounded == id_fractional
        id = id_rounded;
      }
    } else {
      // if -round_basins not set
      id = id_rounded;
    }
    (*basins)(i, j) = id;
  } // end of loop over grid points
}

//! Identify ocean: identify ocean up to continental shelf without detached submarine islands
//! regions icerises: identify grounded regions without detached ice rises

void BoxModel::identifyMASK(IceModelVec2S &inputmask, std::string masktype) {

  m_log->message(4, "0b1: in identifyMASK routine\n");

  int seed_x = (grid.Mx - 1) / 2, seed_y = (grid.My - 1) / 2;

  double linner_identified = 0.0, all_inner_identified = 1.0, previous_step_identified = 0.0;
  ;

  list.add(inputmask);
  inputmask.set(imask_unidentified);
  if ((seed_x >= grid.xs) and (seed_x < grid.xs + grid.xm) and (seed_y >= grid.ys) and (seed_y <= grid.ys + grid.ym)) {
    inputmask(seed_x, seed_y) = imask_inner;
  }

  int iteration_round = 0;
  // find inner region first
  while (all_inner_identified > previous_step_identified) {

    iteration_round += 1;
    previous_step_identified = all_inner_identified;

    list.add(inputmask);
    list.add(*mask);
    list.add(*topg);

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      bool masktype_condition = false;
      if (masktype == "ocean") {
        masktype_condition = ((*mask)(i, j) != maskocean or (*topg)(i, j) >= continental_shelf_depth);
      } else if (masktype == "icerises") {
        masktype_condition = ((*mask)(i, j) == maskgrounded);
      }

      if (masktype_condition and inputmask(i, j) == imask_unidentified and
          (inputmask(i, j + 1) == imask_inner or inputmask(i, j - 1) == imask_inner or
           inputmask(i + 1, j) == imask_inner or inputmask(i - 1, j) == imask_inner)) {
        inputmask(i, j) = imask_inner;
        linner_identified += 1;
      } else if (masktype_condition == false) {
        inputmask(i, j) = imask_outer;
      }
    } // end of loop over grid points

    inputmask.update_ghosts();

    all_inner_identified = GlobalSum(grid.com, linner_identified);
  }

  // set value for excluded areas (ice rises or submarine islands)
  list.add(inputmask);
  list.add(*mask);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (inputmask(i, j) == imask_unidentified) {
      inputmask(i, j) = imask_exclude;
    }

    if (masktype == "ocean") {
      // exclude ice covered parts
      if ((*mask)(i, j) != maskocean and inputmask(i, j) == imask_inner) {
        inputmask(i, j) = imask_outer;
      }
    }
  } // end of loop over grid points
}

//! When ocean_given is set compute mean salinity and temperature in each basin.
void BoxModel::computeOCEANMEANS() {
  // FIXME currently the mean is also calculated over submarine islands which
  // are higher than continental_shelf_depth

  m_log->message(4, "0b2: in computeOCEANMEANS routine \n");

  // count cells to compute mean over for each basin
  std::vector<double> basin_size(numberOfBasins, 0.0);
  // add temperature for each basin
  std::vector<double> basin_temperature(numberOfBasins, 0.0);
  // add salinity for each basin
  std::vector<double> basin_salinity(numberOfBasins, 0.0);

  list.add(OCEANMEANmask);
  list.add(*theta_ocean);
  list.add(*salinity_ocean); // salinity
  list.add(*basins);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (OCEANMEANmask(i, j) == imask_inner) {
      int shelf_id = (*basins)(i, j);

      basin_size[shelf_id]        += 1;
      basin_salinity[shelf_id]    += (*salinity_ocean)(i, j);
      basin_temperature[shelf_id] += (*theta_ocean)(i, j);
    }
  } // end of loop over grid points

  for (int k = 0; k < numberOfBasins; k++) {

    basin_size[k]        = GlobalSum(grid.com, basin_size[k]);
    basin_salinity[k]    = GlobalSum(grid.com, basin_salinity[k]);
    basin_temperature[k] = GlobalSum(grid.com, basin_temperature[k]);

    if (k > 0 and basin_size[k] == 0) {
      // if basin is not dummy basin 0 or there are no ocean cells in this basin
      // to take the mean over.
      m_log->message(2, "PISM_WARNING: basin %d contains no ocean mean cells,\n"
                        "              no mean salinity or temperature values are computed!\n",
                     k);
    } else {
      basin_salinity[k]    = basin_salinity[k] / basin_size[k];
      basin_temperature[k] = basin_temperature[k] / basin_size[k];

      Toc_base_vec[k] = basin_temperature[k] - 273.15;
      Soc_base_vec[k] = basin_salinity[k];
      m_log->message(4, "  %d: temp =%.3f, salinity=%.3f\n", k, Toc_base_vec[k], Soc_base_vec[k]);
    }
  }
}

//! Compute the extent of the ice shelves of each basin/region (i.e. counter)
//and find grounding line and calving front in BOXMODELmask

void BoxModel::extentOfIceShelves() {

  m_log->message(4, "A1b: in extent of ice shelves routine\n");

  double lcounter_box_unidentified = 0; // count the total amount of unidentified shelf boxes.
  double lcounter[numberOfBasins];
  double lcounter_CFbox[numberOfBasins];
  double lcounter_GLbox[numberOfBasins];

  for (int k = 0; k < numberOfBasins; k++) {
    lcounter[k]       = 0.0;
    lcounter_CFbox[k] = 0.0;
    lcounter_GLbox[k] = 0.0;
  }

  list.add(*mask);
  list.add(*basins);
  list.add(BOXMODELmask);
  if (exicerises_set) {
    list.add(ICERISESmask);
  }

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if ((*mask)(i, j) == maskfloating) {
      // if this is a ice shelf cell
      int shelf_id = (*basins)(i, j);
      lcounter[shelf_id]++;

      bool neighbor_to_land;
      if (exicerises_set) {
        neighbor_to_land = (ICERISESmask(i, j + 1) == imask_inner or ICERISESmask(i, j - 1) == imask_inner or
                            ICERISESmask(i + 1, j) == imask_inner or ICERISESmask(i - 1, j) == imask_inner or
                            ICERISESmask(i + 1, j + 1) == imask_inner or ICERISESmask(i + 1, j - 1) == imask_inner or
                            ICERISESmask(i - 1, j + 1) == imask_inner or ICERISESmask(i - 1, j - 1) == imask_inner);
      } else {
        neighbor_to_land = ((*mask)(i, j + 1) < maskfloating or (*mask)(i, j - 1) < maskfloating or
                            (*mask)(i + 1, j) < maskfloating or (*mask)(i - 1, j) < maskfloating or
                            (*mask)(i + 1, j + 1) < maskfloating or (*mask)(i + 1, j - 1) < maskfloating or
                            (*mask)(i - 1, j + 1) < maskfloating or (*mask)(i - 1, j - 1) < maskfloating);
      }

      if (neighbor_to_land) {
        // i.e. there is a grounded neighboring cell (which is not ice rise!)
        BOXMODELmask(i, j) = box_GL;
        lcounter_GLbox[shelf_id]++;

      } else if ((*mask)(i, j + 1) == maskocean or (*mask)(i, j - 1) == maskocean or (*mask)(i + 1, j) == maskocean or
                 (*mask)(i - 1, j) == maskocean or (*mask)(i + 1, j + 1) == maskocean or
                 (*mask)(i + 1, j - 1) == maskocean or (*mask)(i - 1, j + 1) == maskocean or
                 (*mask)(i - 1, j - 1) == maskocean) {
        // i.e. there is an ocean neighboring cell
        BOXMODELmask(i, j) = box_IF;
        lcounter_CFbox[shelf_id]++;
      } else {
        // i.e., all other floating boxes
        BOXMODELmask(i, j) = box_unidentified;
        lcounter_box_unidentified++;
      }

    } else {
      // i.e., not floating
      BOXMODELmask(i, j) = box_noshelf;
    }
  } // end of loop over grid points

  BOXMODELmask.update_ghosts();

  counter_box_unidentified = GlobalSum(grid.com, lcounter_box_unidentified);
  for (int k = 0; k < numberOfBasins; k++) {
    counter[k]       = GlobalSum(grid.com, lcounter[k]);
    counter_CFbox[k] = GlobalSum(grid.com, lcounter_CFbox[k]);
    counter_GLbox[k] = GlobalSum(grid.com, lcounter_GLbox[k]);

    m_log->message(5, "  %d: cnt[k]=%.0f, cnt_CFbox=%.0f, cnt_GLbox=%.0f\n", k, counter[k], counter_CFbox[k],
                   counter_GLbox[k]);
  }
}

//! Compute the BOXMODELmask, calculate the extent of each box in each region
void BoxModel::identifyBOXMODELmask() {

  m_log->message(4, "A1c: in identify box model mask routine\n");

  double lcounter_box_unidentified = counter_box_unidentified + 1.0;

  while ((counter_box_unidentified > 0.0) and (lcounter_box_unidentified != counter_box_unidentified)) {

    lcounter_box_unidentified = counter_box_unidentified;
    m_log->message(5, "     cnt_box_unidentified=%.0f\n", lcounter_box_unidentified);

    for (int l = 0; l < 3; l++) {
      // FIXME size depends on how often this routine is called
      extendIFBox();
    }
    extendGLBox();
  }

  // FIXME How to handle shelf-lakes (shelf which lies in the sheet), at the
  // moment: GL_Box, better: Beckmann-Goose or no melting at all?
  if (counter_box_unidentified > 0.0) {
    // if there are still floating cells which were not labels before, i.e. a
    // shelf connected to an island (in the ex_icerises case)

    list.add(BOXMODELmask);
    list.add(*basins);
    double lcounter_GLbox[numberOfBasins];
    double all_counter_GLbox[numberOfBasins];
    for (int k = 0; k < numberOfBasins; k++) {
      lcounter_GLbox[k]    = 0.0;
      all_counter_GLbox[k] = 0.0;
    }

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (BOXMODELmask(i, j) == box_unidentified) {
        // i.e. this is an unidentified shelf cell with a neighbor that is in the IFbox
        m_log->message(5, "   left over i=%d, j=%d \n", i, j);
        BOXMODELmask(i, j) = box_GL;
        int shelf_id = (*basins)(i, j);
        lcounter_GLbox[shelf_id]++;
      }
    } // end of loop over grid points


    for (int k = 0; k < numberOfBasins; k++) {
      all_counter_GLbox[k] = GlobalSum(grid.com, lcounter_GLbox[k]);
      counter_GLbox[k] += all_counter_GLbox[k];
    }
  }

  for (int k = 0; k < numberOfBasins; k++) {
    m_log->message(
        5, "  %d: cnt[i] = %.0f, cnt_CFbox = %.0f, cnt_GLbox = %.0f, ratio_CF_box = %.3f, ratio_GL_box = %.3f\n", k,
        counter[k], counter_CFbox[k], counter_GLbox[k], counter_CFbox[k] / counter[k], counter_GLbox[k] / counter[k]);
  }
}

//! extent the grounding line box with neighboring unidentified shelf cells
void BoxModel::extendGLBox() {

  double lcounter_box_unidentified = 0.0, all_counter_box_unidentified = 0.0;
  double lcounter_GLbox[numberOfBasins], all_counter_GLbox[numberOfBasins];

  for (int k = 0; k < numberOfBasins; k++) {
    lcounter_GLbox[k]    = 0.0;
    all_counter_GLbox[k] = 0.0;
  }

  list.add(*mask);
  list.add(*basins);
  list.add(BOXMODELmask);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (BOXMODELmask(i, j) == box_unidentified and
        (BOXMODELmask(i, j + 1) == box_GL or BOXMODELmask(i, j - 1) == box_GL or BOXMODELmask(i + 1, j) == box_GL or
         BOXMODELmask(i - 1, j) == box_GL)) {
      // i.e. this is an unidentified shelf cell with a neighbor that is in the GLbox
      BOXMODELmask(i, j) = box_neighboring;
    }
  } // end of loop over grid points

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (BOXMODELmask(i, j) == box_neighboring) {
      BOXMODELmask(i, j) = box_GL;
      lcounter_box_unidentified++;
      int shelf_id = (*basins)(i, j);
      lcounter_GLbox[shelf_id]++;
    }
  } // end of loop over grid points


  BOXMODELmask.update_ghosts();

  all_counter_box_unidentified = GlobalSum(grid.com, lcounter_box_unidentified);
  counter_box_unidentified -= all_counter_box_unidentified;

  for (int k = 0; k < numberOfBasins; k++) {
    all_counter_GLbox[k] = GlobalSum(grid.com, lcounter_GLbox[k]);
    counter_GLbox[k] += all_counter_GLbox[k];
  }
}

//! extend the ice_front box with neighboring unidentified shelf cells.
void BoxModel::extendIFBox() {

  double lcounter_box_unidentified    = 0.0;
  double all_counter_box_unidentified = 0.0;
  double lcounter_CFbox[numberOfBasins];
  double all_counter_CFbox[numberOfBasins];
  for (int k = 0; k < numberOfBasins; k++) {
    lcounter_CFbox[k]    = 0.0;
    all_counter_CFbox[k] = 0.0;
  }

  list.add(*mask);
  list.add(*basins);
  list.add(BOXMODELmask);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (BOXMODELmask(i, j) == box_unidentified and
        (BOXMODELmask(i, j + 1) == box_IF or BOXMODELmask(i, j - 1) == box_IF or BOXMODELmask(i + 1, j) == box_IF or
         BOXMODELmask(i - 1, j) == box_IF)) {
      // i.e. this is an unidentified shelf cell with a neighbor that is in the IFbox
      BOXMODELmask(i, j) = box_neighboring;
    }
  } // end of loop over grid points

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (BOXMODELmask(i, j) == box_neighboring) {
      BOXMODELmask(i, j) = box_IF;
      lcounter_box_unidentified++;
      int shelf_id = (*basins)(i, j);
      lcounter_CFbox[shelf_id]++;
    }
  } // end of loop over grid points

  BOXMODELmask.update_ghosts();

  all_counter_box_unidentified = GlobalSum(grid.com, lcounter_box_unidentified);
  counter_box_unidentified -= all_counter_box_unidentified;

  for (int k = 0; k < numberOfBasins; k++) {
    all_counter_CFbox[k] = GlobalSum(grid.com, lcounter_CFbox[k]);
    counter_CFbox[k] += all_counter_CFbox[k];
  }
}

/*!
  Compute ocean temperature outside of the ice shelf cavities.
*/

void BoxModel::oceanTemperature(const POBMConstants &cc) {

  m_log->message(4, "A2 : in ocean temp routine\n");

  list.add(*mask);
  list.add(*basins);
  list.add(*ice_thickness);

  list.add(Soc_base);
  list.add(Toc_base);
  list.add(Toc_anomaly);
  list.add(Toc);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // make sure all temperatures are zero at the beginning of each time step
    Toc(i, j)         = 273.15; // in Kelvin
    Toc_base(i, j)    = 273.15; // in Kelvin
    Toc_anomaly(i, j) = 0.0;    // in Kelvin or Celsius
    Soc_base(i, j)    = 0.0;    // in psu

    if ((*mask)(i, j) == maskfloating) {
      int shelf_id = (*basins)(i, j);
      Toc_base(i, j) = 273.15 + Toc_base_vec[shelf_id];
      Soc_base(i, j) = Soc_base_vec[shelf_id];

      //! salinity and temperature for grounding line box
      if (Soc_base(i, j) == 0.0 or Toc_base_vec[shelf_id] == 0.0) {
        m_log->message(2, "PISM_ERROR: Missing Soc_base and Toc_base for %d, %d, basin %d \n   Aborting... \n", i, j,
                       shelf_id);
        PISMEnd();
      }

      // Add temperature anomalies from given nc-file  // FIXME different
      // nc-files for each basin!
      if ((delta_T != NULL) and ocean_oceanboxmodel_deltaT_set) {
        // Toc_anomaly(i, j) = delta_T_factor * (*delta_T)(m_t + 0.5*m_dt);
        Toc_anomaly(i, j) = delta_T_factor * temp_anomaly;

      } else {

        Toc_anomaly(i, j) = 0.0;
      }

      ////////////////
      // prevent ocean temp from being below pressure melting temperature

      // const double shelfbaseelev = - (cc.rhoi / cc.rhow) *
      // (*ice_thickness)(i, j);

      const double pressure = cc.rhoi * cc.earth_grav * (*ice_thickness)(i, j) * 1e-4;
      // MUST be in dbar  // NOTE 1dbar = 10000 Pa = 1e4 kg m-1 s-2,
      // FIXME need to include atmospheric pressure?
      const double T_pmt = cc.a * Soc_base(i, j) + cc.b - cc.c * pressure;

      Toc_anomaly(i, j) = PetscMax(T_pmt + 273.15 - Toc_base(i, j), Toc_anomaly(i, j));
      ////////////////

      Toc(i, j) = Toc_base(i, j) + Toc_anomaly(i, j); // in Kelvin

    } // end if herefloating
  }   // end of loop over grid points
}

// NOTE Mean Gl_box melt rate is needed for basalMeltRateForIceFrontBox(). Here,
// mean is taken over all shelves in one drainage basin!

//! Compute the basal melt / refreezing rates for each shelf cell bordering the
//grounding line box
void BoxModel::basalMeltRateForGroundingLineBox(const POBMConstants &cc) {
  m_log->message(4, "B1 : in basal melt rate gl routine\n");

  double lcounter_edge_of_GLbox_vector[numberOfBasins], lmean_salinity_GLbox_vector[numberOfBasins],
      lmean_meltrate_GLbox_vector[numberOfBasins], lmean_overturning_GLbox_vector[numberOfBasins];

  for (int k = 0; k < numberOfBasins; k++) {
    lcounter_edge_of_GLbox_vector[k]  = 0.0;
    lmean_salinity_GLbox_vector[k]    = 0.0;
    lmean_meltrate_GLbox_vector[k]    = 0.0;
    lmean_overturning_GLbox_vector[k] = 0.0;
  }

  list.add(*ice_thickness);
  list.add(*basins);
  list.add(BOXMODELmask);
  list.add(T_star);
  list.add(Toc_base);
  list.add(Toc_anomaly);
  list.add(Toc_inCelsius);
  list.add(Toc);
  list.add(Soc_base);
  list.add(Soc);
  list.add(overturning);
  list.add(basalmeltrate_shelf);

  double countHelpterm = 0, lcountHelpterm = 0;

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    int shelf_id = (*basins)(i, j);

    // Make sure everything is at default values at the beginning of each
    // time step
    T_star(i, j)        = 0.0; // in Celsius
    Toc_inCelsius(i, j) = 0.0; // in Celsius
    Soc(i, j)           = 0.0; // in psu

    basalmeltrate_shelf(i, j) = 0.0;
    overturning(i, j)         = 0.0;

    if (BOXMODELmask(i, j) == box_GL and shelf_id > 0.0) {

      const double pressure = cc.rhoi * cc.earth_grav * (*ice_thickness)(i, j) * 1e-4; // MUST be in dbar
      // NOTE 1dbar = 10000 Pa = 1e4 kg m-1 s-2
      // FIXME need to include atmospheric pressure?
      T_star(i, j) =
        cc.a * Soc_base(i, j) + cc.b - cc.c * pressure - (Toc_base(i, j) - 273.15 + Toc_anomaly(i, j)); // in Celsius

      double gamma_T_star = gamma_T_star_vec[shelf_id], C1 = C_vec[shelf_id],
        g1 = ((counter_GLbox[shelf_id] * grid.dx * grid.dy) * gamma_T_star / (C1 * cc.rho_star));

      //! temperature for grounding line box

      double helpterm1 = (g1 / (cc.beta * (Soc_base(i, j) / (cc.nu * cc.lambda)) - cc.alpha)); // in 1 / (1/Celsius) = Celsius
      double helpterm2 = ((g1 * T_star(i, j)) /
                          (cc.beta * (Soc_base(i, j) / (cc.nu * cc.lambda)) - cc.alpha)); // in Celsius / (1/Celsius) = Celsius^2

      if ((0.25 * PetscSqr(helpterm1) - helpterm2) < 0.0) {
        helpterm2 = 0.25 * PetscSqr(helpterm1);
        // FIXME: This might be wrong!
        lcountHelpterm += 1;
      }

      // NOTE Careful, Toc_base(i, j) is in Kelvin, Toc_inCelsius(i, j) NEEDS to be
      // in Celsius!
      Toc_inCelsius(i, j) = ((Toc_base(i, j) - 273.15 + Toc_anomaly(i, j)) -
                             (-0.5 * helpterm1 + sqrt(0.25 * PetscSqr(helpterm1) - helpterm2)));

      //! salinity for grounding line box
      Soc(i, j) = (Soc_base(i, j) -
                   (Soc_base(i, j) / (cc.nu * cc.lambda)) *
                   ((Toc_base(i, j) - 273.15 + Toc_anomaly(i, j)) - Toc_inCelsius(i, j))); // in psu

      //! basal melt rate for grounding line box
      basalmeltrate_shelf(i, j) = ((-gamma_T_star / (cc.nu * cc.lambda)) *
                                   (cc.a * Soc(i, j) + cc.b - cc.c * pressure - Toc_inCelsius(i, j))); // in m/s

      //! overturning
      //
      // NOTE Actually, there is of course no overturning-FIELD, it is only a scalar for each
      // shelf.
      //
      // Here, I compute overturning as
      //
      // MEAN[C1*cc.rho_star* (cc.beta*(Soc_base(i, j)-Soc(i, j)) - cc.alpha*((Toc_base(i,
      // j)-273.15+Toc_anomaly(i, j))-Toc_inCelsius(i, j)))]
      //
      // while in fact it should be
      //
      // C1*cc.rho_star* (cc.beta*(Soc_base-MEAN[Soc(i, j)]) -
      // cc.alpha*((Toc_base-273.15+Toc_anomaly)-MEAN[Toc_inCelsius(i, j)]))
      //
      // which is the SAME since Soc_base, Toc_base and Toc_anomaly are the
      // same FOR ALL i, j CONSIDERED, so this is just nomenclature!
      overturning(i, j) =
        C1 * cc.rho_star *
        (cc.beta * (Soc_base(i, j) - Soc(i, j)) -
         cc.alpha * ((Toc_base(i, j) - 273.15 + Toc_anomaly(i, j)) - Toc_inCelsius(i, j))); // in m^3/s

      if (BOXMODELmask(i - 1, j) == box_IF or BOXMODELmask(i + 1, j) == box_IF or BOXMODELmask(i, j - 1) == box_IF or
          BOXMODELmask(i, j + 1) == box_IF) {
        // i.e., if this cell is from the GL box and one of the neighbors is from the CF box - It
        // is important to only take the border of the grounding line box to the calving front box
        // into account, because the following mean value will be used to compute the value for
        // the calving front box. I.e., this helps avoiding discontinuities!

        lcounter_edge_of_GLbox_vector[shelf_id] += 1;
        lmean_salinity_GLbox_vector[shelf_id] += Soc(i, j);
        lmean_meltrate_GLbox_vector[shelf_id] += basalmeltrate_shelf(i, j);
        lmean_overturning_GLbox_vector[shelf_id] += overturning(i, j);
      }
      // no else-case necessary since all variables are set to zero at the beginning of this
      // routine

    } else { // i.e., not GL_box
      basalmeltrate_shelf(i, j) = 0.0;
    }
  } // end of loop over grid points

  for (int k = 0; k < numberOfBasins; k++) {
    double counter_edge_of_GLbox_vector = 0.0;
    counter_edge_of_GLbox_vector        = GlobalSum(grid.com, lcounter_edge_of_GLbox_vector[k]);
    mean_meltrate_GLbox_vector[k]       = GlobalSum(grid.com, lmean_meltrate_GLbox_vector[k]);
    mean_salinity_GLbox_vector[k]       = GlobalSum(grid.com, lmean_salinity_GLbox_vector[k]);
    mean_overturning_GLbox_vector[k]    = GlobalSum(grid.com, lmean_overturning_GLbox_vector[k]);

    if (counter_edge_of_GLbox_vector > 0.0) {
      mean_salinity_GLbox_vector[k]    = mean_salinity_GLbox_vector[k] / counter_edge_of_GLbox_vector;
      mean_meltrate_GLbox_vector[k]    = mean_meltrate_GLbox_vector[k] / counter_edge_of_GLbox_vector;
      mean_overturning_GLbox_vector[k] = mean_overturning_GLbox_vector[k] / counter_edge_of_GLbox_vector;
    } else { // This means that there is no [cell from the GLbox neighboring a
             // cell from the CFbox], NOT necessarily that there is no GLbox!
      mean_salinity_GLbox_vector[k]    = 0.0;
      mean_meltrate_GLbox_vector[k]    = 0.0;
      mean_overturning_GLbox_vector[k] = 0.0;
    }

    m_log->message(5, "  %d: cnt=%.0f, sal=%.3f, melt=%.3e, over=%.1e \n", k, counter_edge_of_GLbox_vector,
                   mean_salinity_GLbox_vector[k], mean_meltrate_GLbox_vector[k], mean_overturning_GLbox_vector[k]);
  }

  countHelpterm = GlobalSum(grid.com, lcountHelpterm);

  if (countHelpterm > 0) {
    m_log->message(2, "B1!: PISM_WARNING: square-root has been negative in %.0f cases!\n", countHelpterm);
  }
}

//! Compute the basal melt / refreezing rates for each shelf cell bordering the ice front box
void BoxModel::basalMeltRateForIceFrontBox(const POBMConstants &cc) {

  // FIXME redo verbprintfs!
  m_log->message(4, "B2 : in bm ice front routine\n");

  list.add(*ice_thickness);
  list.add(*basins);
  list.add(BOXMODELmask);

  list.add(T_star);
  list.add(Toc_base);
  list.add(Toc_anomaly);
  list.add(Toc_inCelsius);
  list.add(Toc);
  list.add(Soc_base);
  list.add(Soc);
  list.add(overturning);
  list.add(basalmeltrate_shelf);

  double countk4 = 0, lcountk4 = 0, countGl0 = 0, lcountGl0 = 0, countSqr = 0, lcountSqr = 0, countMean0 = 0,
         lcountMean0 = 0;

  //! The ice front box = BOX I
  for (Points p(*m_grid); p; p.next()) { // FIXME REPAIR
    const int i = p.i(), j = p.j();

    int shelf_id = (*basins)(i, j);

    if (BOXMODELmask(i, j) == box_IF and shelf_id > 0.0) {

      const double pressure = cc.rhoi * cc.earth_grav * (*ice_thickness)(i, j) * 1e-4; // MUST be in dbar
      // NOTE 1dbar = 10000 Pa = 1e4 kg m-1 s-2
      // FIXME need to include atmospheric pressure?
      T_star(i, j) =
        cc.a * Soc_base(i, j) + cc.b - cc.c * pressure - (Toc_base(i, j) - 273.15 + Toc_anomaly(i, j)); // in Celsius

      double gamma_T_star, area_GLbox, area_CFbox, mean_salinity_in_GLbox, mean_meltrate_in_GLbox,
        mean_overturning_in_GLbox;

      gamma_T_star              = gamma_T_star_vec[shelf_id];
      area_CFbox                = (counter_CFbox[shelf_id] * grid.dx * grid.dy);
      area_GLbox                = (counter_GLbox[shelf_id] * grid.dx * grid.dy);
      mean_salinity_in_GLbox    = mean_salinity_GLbox_vector[shelf_id];
      mean_meltrate_in_GLbox    = mean_meltrate_GLbox_vector[shelf_id];
      mean_overturning_in_GLbox = mean_overturning_GLbox_vector[shelf_id];

      if (area_GLbox == 0) {
        // if there is no grounding line box in the current basin, set BOXMODELmask to 3 and
        // compute basal melt rate by Beckmann-Gose
        BOXMODELmask(i, j) = box_other;
        lcountGl0 += 1;

      } else {

        double k1 = (area_CFbox * gamma_T_star) / (cc.nu * cc.lambda);
        // in (m^2*m/s)/(Celsius) = m^3/(s*Celsius)
        double k2 = 1 / (mean_overturning_in_GLbox + area_CFbox * gamma_T_star);
        // in s/m^3
        double k3 = (area_CFbox * gamma_T_star * T_star(i, j) - cc.nu * cc.lambda * area_GLbox * mean_meltrate_in_GLbox);
        // in m^2 * m/s * Celsius = m^3 * Celsius / s
        double k4 = (-k1 * k2 * area_CFbox * gamma_T_star * cc.a + k1 * cc.a);
        // in m^3/(s*Celsius) * s/m^3 * m^2 * m/s * Celsius/psu = m^3/(s*psu)
        double k5 = (mean_overturning_in_GLbox + Soc_base(i, j) * k1 * k2 * area_CFbox * gamma_T_star * cc.a -
                     k1 * Soc_base(i, j) * cc.a - k1 * T_star(i, j) + k1 * k2 * k3);
        // in m^3/s
        double k6 = (k1 * Soc_base(i, j) * T_star(i, j) - k1 * k2 * Soc_base(i, j) * k3 -
                     area_GLbox * mean_meltrate_in_GLbox * mean_salinity_in_GLbox);
        // in psu * m^3/s

        //! salinity for calving front box
        if (k4 == 0.0) {
          // FIXME rewrite this warning? Do not stop but calculate melt rates
          // according to Beckmann-Gose?
          lcountk4 += 1;
          BOXMODELmask(i, j) = box_other;
          continue;
        }

        if ((0.25 * PetscSqr(k5 / k4) - (k6 / k4)) < 0.0) {
          lcountSqr += 1;
          BOXMODELmask(i, j) = box_other;
          continue;
        }

        Soc(i, j) = (Soc_base(i, j) - (-0.5 * (k5 / k4) - sqrt(0.25 * PetscSqr(k5 / k4) - (k6 / k4)))); // in psu

        //! temperature for calving front box
        // NOTE Careful, Toc_base(i, j) is in Kelvin, Toc_inCelsius(i, j) NEEDS to
        // be in Celsius!
        Toc_inCelsius(i, j) = (Toc_base(i, j) - 273.15 + Toc_anomaly(i, j)) -
          (k2 * area_CFbox * gamma_T_star * cc.a * (Soc_base(i, j) - Soc(i, j)) - k2 * k3);

        //! basal melt rate for calving front box
        basalmeltrate_shelf(i, j) = ((-gamma_T_star / (cc.nu * cc.lambda)) *
                                     (cc.a * Soc(i, j) + cc.b - cc.c * pressure - Toc_inCelsius(i, j))); // in m/s

        if (mean_salinity_in_GLbox == 0.0 or mean_meltrate_in_GLbox == 0.0 or mean_overturning_in_GLbox == 0.0) {
          // this must not occur since there must always be a GL_box neighbor
          lcountMean0 += 1;
          BOXMODELmask(i, j) = box_other;
          continue;
        }
      }
    }
    // NOTE NO else-case, since basalMeltRateForGroundingLineBox() and
    // basalMeltRateForOtherShelves() cover all other cases and we would overwrite those results
    // here.
  } // end of loop over grid points

  countk4    = GlobalSum(grid.com, lcountk4);
  countGl0   = GlobalSum(grid.com, lcountGl0);
  countSqr   = GlobalSum(grid.com, lcountSqr);
  countMean0 = GlobalSum(grid.com, lcountMean0);

  if (countk4 > 0) {
    m_log->message(2, "B2!: PISM_WARNING: k4 is zero in %.0f case(s)!\n", countk4);
  }

  if (countGl0 > 0) {
    m_log->message(2, "B2!: PISM_WARNING: no grounding line box in basin in %.0f case(s)!\n", countGl0);
  }

  if (countSqr > 0) {
    m_log->message(2, "B2!: PISM_WARNING: square root is negative in %.0f case(s)!\n", countSqr);
  }

  if (countMean0 > 0) {
    m_log->message(2, "B2!: PISM_WARNING: mean of salinity, melt rate or overturning is zero in %.0f case(s)!\n",
                   countMean0);
  }
}

//! Convert Toc_inCelsius from Celsius to Kelvin and write into Toc for the .nc-file; NOTE
//It is crucial, that Toc_inCelsius is in Celsius for the computation of the basal
//melt rate
//! Compute the melt rate for all other ice shelves.
void BoxModel::basalMeltRateForOtherShelves(const POBMConstants &cc) {
  m_log->message(4, "B3 : in bm others routine\n");

  list.add(*ice_thickness);
  list.add(*basins);
  list.add(BOXMODELmask);
  list.add(Toc_base);
  list.add(Toc_anomaly);
  list.add(Toc_inCelsius);
  list.add(Toc);
  list.add(overturning);
  list.add(basalmeltrate_shelf); // NOTE melt rate has units:   J m-2 s-1 / (J
                                 // kg-1 * kg m-3) = m s-1
  list.add(heatflux);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    int shelf_id = (*basins)(i, j);

    if (shelf_id == 0) { // boundary of computational domain

      basalmeltrate_shelf(i, j) = 0.0;

    } else if (BOXMODELmask(i, j) == box_other) {

      Toc(i, j) = Toc_base(i, j) + Toc_anomaly(i, j);
      // in Kelvin, NOTE: Toc_base is already in Kelvin, so no (+273.15)
      // default: compute the melt rate from the temperature field according
      // to beckmann_goosse03 (see below)

      const double shelfbaseelev = -(cc.rhoi / cc.rhow) * (*ice_thickness)(i, j);

      // FIXME: for consistency reasons there should be constants a, b, c,
      // gamma_T used
      double T_f = 273.15 + (cc.a * cc.meltSalinity + cc.b2 + cc.c * shelfbaseelev); // add 273.15 to get it in
      // Kelvin... 35 is the
      // salinity

      heatflux(i, j) = cc.meltFactor * cc.rhow * cc.c_p_ocean * cc.gamma_T_o * (Toc(i, j) - T_f); // in W/m^2
      basalmeltrate_shelf(i, j) = heatflux(i, j) / (cc.latentHeat * cc.rhoi);                     // in m s-1

    } else if (shelf_id > 0.0) {

      Toc(i, j) = 273.15 + Toc_inCelsius(i, j) + Toc_anomaly(i, j); // in Kelvin
      // FIXME I think Toc should not occur in any of the routines before!

    } else { // This must not happen

      m_log->message(2, "PISM_ERROR: [rank %d] at %d, %d  -- basins(i, j)=%d causes problems.\n   Aborting... \n",
                     grid.rank, i, j, shelf_id);
      PISMEnd();
    }
  } // end of loop over grid points
}

void BoxModel::shelf_base_temperature(IceModelVec2S &result) {
  shelfbtemp.copy_to(result);
}

void BoxModel::shelf_base_mass_flux(IceModelVec2S &result) {
  shelfbmassflux.copy_to(result);
}

void BoxModel::sea_level_elevation(double &result) {
  result = sea_level;
}

void BoxModel::melange_back_pressure_fraction(IceModelVec2S &result) {
  result.set(0.0);
}

void BoxModel::define_variables(std::set<std::string> vars, const PIO &nc, IO_Type nctype) {

  PGivenClimate<POModifier, OceanModel>::define_variables(vars, nc, nctype);

  for (unsigned int k = 0; k < m_variables.size(); ++k) {
    IceModelVec *v = m_variables[k];
    std::string name = v->metadata().get_string("short_name");
    if (set_contains(vars, name)) {
      v->define(nc, nctype);
    }
  }
}

void BoxModel::write_variables(std::set<std::string> vars, const PIO &nc) {

  PGivenClimate<POModifier, OceanModel>::write_variables(vars, nc);

  for (unsigned int k = 0; k < m_variables.size(); ++k) {
    IceModelVec *v = m_variables[k];
    std::string name = v->metadata().get_string("short_name");
    if (set_contains(vars, name)) {
      v->write(nc);
    }
  }
}

} // end of namespace ocean
} // end of namespace pism
