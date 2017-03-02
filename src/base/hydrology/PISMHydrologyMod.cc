// Copyright (C) 2012-2016 PISM Authors
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

#include "PISMHydrologyMod.hh"
#include "base/util/Mask.hh"
#include "base/util/PISMVars.hh"
#include "base/util/error_handling.hh"
#include "base/util/iceModelVec2T.hh"
#include "base/util/io/PIO.hh"
#include "base/util/pism_options.hh"
#include "hydrology_diagnostics.hh"
//#include "hydrology_diagnostics_mod.hh"
#include "base/util/pism_utilities.hh"
#include "base/util/PISMVars.hh"
#include "base/util/MaxTimestep.hh"
#include "base/util/IceModelVec2CellType.hh"
#include <math.h>
#include <stdlib.h>

namespace pism {
namespace hydrology {

class hydrology;

HydrologyMod::HydrologyMod(IceGrid::ConstPtr g)
  : Hydrology(g) {

  m_fraction_till.create(m_grid, "fraction_till", WITHOUT_GHOSTS);
  m_fraction_till.set_attrs("model_state",
                       "fraction of surface that is covered with thick till",
                       "1", "");
  m_fraction_till.metadata().set_double("valid_min", 0.0);
  m_fraction_channel.create(m_grid, "fraction_channel", WITHOUT_GHOSTS);
  m_fraction_channel.set_attrs("model_state",
                       "fraction of surface that is covered with channels",
                       "1", "");
  m_fraction_channel.metadata().set_double("valid_min", 0.0);

  m_till_permeability.create(m_grid, "till_permeability", WITHOUT_GHOSTS);
  m_till_permeability.set_attrs("model_state",
                       "hydraulic permeability of the till",
                       "m s-1", "");
  m_till_permeability.metadata().set_double("valid_min", 0.0);

  m_Pover_ghosts.create(m_grid, "overburden_pressure_internal",WITH_GHOSTS,1);
  m_Pover_ghosts.set_attrs("internal",
                  "overburden pressure",
                  "Pa", "");
  m_Pover_ghosts.metadata().set_double("valid_min", 0.0);


  m_pressure_gradient.create(m_grid, "pressure_gradient",WITH_GHOSTS,1);
  m_pressure_gradient.set_attrs("model_state",
                       "pressure gradient at the base",
                       "Pa m-1", "");


  m_gradient_magnitude.create(m_grid, "gradient_magnitude", WITHOUT_GHOSTS);
  m_gradient_magnitude.set_attrs("model_state",
                       "pressure gradient magnitude at the base",
                       "Pa m-1", "");
  m_gradient_magnitude.metadata().set_double("valid_min", 0.0);

  m_tillwat_flux.create(m_grid, "tillwat_flux", WITH_GHOSTS,1);
  m_tillwat_flux.set_attrs("model_state",
                       "flux of water through the till",
                       "m s-1", "");
  m_tillwat_flux.metadata().set_double("valid_min", 0.0);
//  m_grid->variables().add(m_tillwat_flux);

  m_excess_water.create(m_grid, "excess_water", WITHOUT_GHOSTS); // without ghosts for now, maybe change later
  m_excess_water.set_attrs("model_state",
                       "extra water flux from water not in till",
                       "m s-1", "");
  m_excess_water.metadata().set_double("valid_min", 0.0);
//  m_grid->variables().add(m_excess_water);

  m_theta.create(m_grid, "water_flow_direction", WITH_GHOSTS,1);
  m_theta.set_attrs("internal",
                       "direction of the negative gradient vector, gives direction of water flow",
                       "radians", "");

} // end constructor

HydrologyMod::~HydrologyMod() {
} // end destructor


void HydrologyMod::init() {

  // first initialize the normal Hydrology parameters
  Hydrology::init();


  m_log->message(2,
             "* Including modifications to subglacial hydrology model ...\n");


  // now initialize the fraction_till 


  m_log->message(2, "* Initializing the sediment thickness cover...\n");
  const double default_fraction_till = m_config->get_double("hydrology.default_fraction_till");

  InputOptions opts = process_input_options(m_grid->com);

  if (opts.type == INIT_RESTART) {
    m_fraction_till.read(opts.filename, opts.record);
  } else if (opts.type == INIT_BOOTSTRAP) {
    m_fraction_till.regrid(opts.filename, OPTIONAL, default_fraction_till);
  } else {
    m_fraction_till.set(default_fraction_till);
  }

  // regrid if requested, regardless of how initialized
  regrid("HydrologyMod", m_fraction_till);




  double till_permeability_default = m_config->get_double("hydrology.default_till_permeability");

  switch (opts.type) {
  case INIT_RESTART:
    m_till_permeability.read(opts.filename, opts.record);
    break;
  case INIT_BOOTSTRAP:
    m_till_permeability.regrid(opts.filename, OPTIONAL, till_permeability_default);
    break;
  case INIT_OTHER:
  default:
    m_till_permeability.set(till_permeability_default);
  }

  // whether or not we could initialize from file, we could be asked to regrid from file
  regrid("Hydrology", m_till_permeability);


  // set excess water to zero for now, probably should have options later

  switch (opts.type) {
  case INIT_RESTART:
    m_excess_water.read(opts.filename, opts.record);
    break;
  case INIT_BOOTSTRAP:
    m_excess_water.regrid(opts.filename, OPTIONAL, 0);
    break;
  case INIT_OTHER:
  default:
    m_excess_water.set(0);
  }

  // whether or not we could initialize from file, we could be asked to regrid from file
  regrid("Hydrology", m_excess_water);


} // end init()



//Override the base class get_input_rate to include the last time step water
//
//! Compute the total water input rate into the basal hydrology layer in the ice-covered region, allowing time-varying input from a file.
/*!
The user can specify the total of en- and supra-glacial drainage contributions
to subglacial hydrology in a time-dependent input file using option -hydrology_input_to_bed.
This method includes that possible input along with `bmelt` to get the total water
input into the subglacial hydrology.

This method crops the input rate to the ice-covered region.  It
also uses hydrology_const_bmelt if that is requested.

Call this method using the current \e hydrology time step.  This method
may be called many times per IceModel time step.  See update() method
in derived classes of Hydrology.
 */
void HydrologyMod::get_input_rate(double hydro_t, double hydro_dt,
                               IceModelVec2S &result) {

  m_log->message(2, "* Entering get_input_rate\n");
  bool   use_const   = m_config->get_boolean("hydrology.use_const_bmelt");
  double const_bmelt = m_config->get_double("hydrology.const_bmelt");

  const IceModelVec2S        &bmelt = *m_grid->variables().get_2d_scalar("bmelt");
  const IceModelVec2CellType &mask  = *m_grid->variables().get_2d_cell_type("mask");

  IceModelVec::AccessList list{&m_excess_water};

  if (not m_hold_bmelt) {
    m_bmelt_local.copy_from(bmelt);
  }

  list.add(m_bmelt_local);
  list.add(mask);
  list.add(result);
  if (m_inputtobed != NULL) {
    m_inputtobed->update(hydro_t, hydro_dt);
    m_inputtobed->interp(hydro_t + hydro_dt/2.0);
    list.add(*m_inputtobed);
  }

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (mask.icy(i, j)) {
      result(i,j) = (use_const) ? const_bmelt : m_bmelt_local(i,j);
      if (m_inputtobed != NULL) {
        result(i,j) += (*m_inputtobed)(i,j);
      }
      result(i,j) += m_excess_water(i,j); // added
    } else {
      result(i,j) = 0.0;
    }
  }

  m_log->message(2, "* finished get_input_rate\n");
}


// Update the basal hydrology
// taken from NullTransportHydrology, but modified for my own purposes
void HydrologyMod::update_impl(double t, double dt) {


  m_log->message(2,
             "* calculating hydrology ...\n");


  double m_t = t;
  double m_dt = dt;

  IceModelVec::AccessList list{&m_tillwat_flux, &m_pressure_gradient, &m_total_input,&m_theta,&m_gradient_magnitude, &m_total_input};

  get_input_rate(t, dt, m_total_input); // retrieve the calculated input due to basal melting, etc



  

  pressure_gradient(m_pressure_gradient,m_gradient_magnitude,m_theta); // find the basal pressure gradient

  // calculate till water flux with the given gradient magnitude

  till_drainage(m_tillwat_flux, dt);

  // with the given input rate and till drainage, fill up the till with water

  m_log->message(2, "* something something mask...\n");

  const IceModelVec2CellType &mask = *m_grid->variables().get_2d_cell_type("mask");

  double tillwat_max = m_config->get_double("hydrology.tillwat_max");

  list.add(mask);
  list.add(m_Wtil);
  list.add(m_fraction_till);
  list.add(m_excess_water);

  m_log->message(2, "* starting to calculate the water in the till...\n");
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (mask.ocean(i,j) ) {
      m_log->message(2, "* ocean, should be tillwat_max...\n");
      // the other models say that till in the ocean will be have zero water layer thickness, of course it should be fully saturated, and in my scheme this is a necessary condition
      m_Wtil(i, j) = tillwat_max;
    } else if( mask.ice_free(i,j)) {
      m_log->message(2, "* ice free, should be 0...\n");
      // land without ice will be dry, will need to change it when there are lakes added
      m_Wtil(i,j) = 0.0;
    } else {
      m_log->message(2, "* ice, calculating water in till...\n");
      m_Wtil(i, j) += dt * (m_total_input(i, j) -m_tillwat_flux(i,j)) / m_fraction_till(i,j);
      if(m_Wtil(i,j) < 0) { // should be zero
      m_log->message(2, "* till water is zero...\n");
       m_Wtil(i, j) = 0;
       m_excess_water(i,j) = 0;
      } else if (m_Wtil(i,j)/ m_fraction_till(i,j) > tillwat_max) { // add the extra water to m_excess_water
       m_log->message(2, "* till water is full, put excess in m_excess_water...\n");
       m_excess_water(i,j) = m_Wtil(i,j) - tillwat_max * m_fraction_till(i,j);
       m_Wtil(i,j) = tillwat_max;
      } else { // no excess water and the till is not full
       m_log->message(2, "* till water is not full...\n");
       m_excess_water(i,j) = 0;

      }

    }
  }


  m_log->message(2,
             "* finished calculating hydrology ...\n");


}



void HydrologyMod::define_model_state_impl(const PIO &output) const {
  Hydrology::define_model_state_impl(output);
  m_tillwat_flux.define(output);
  m_excess_water.define(output);
}

void HydrologyMod::write_model_state_impl(const PIO &output) const {
  Hydrology::write_model_state_impl(output);
  m_tillwat_flux.write(output);
  m_excess_water.write(output);
}


/*
std::map<std::string, Diagnostic::Ptr> HydrologyMod::diagnostics_impl() const {
  std::map<std::string, Diagnostic::Ptr> result = {
    {"bwat",       Diagnostic::Ptr(new Hydrology_bwat(this))},
    {"bwp",        Diagnostic::Ptr(new Hydrology_bwp(this))},
    {"bwprel",     Diagnostic::Ptr(new Hydrology_bwprel(this))},
    {"effbwp",     Diagnostic::Ptr(new Hydrology_effbwp(this))},
    {"hydrobmelt", Diagnostic::Ptr(new Hydrology_hydrobmelt(this))},
    {"hydroinput", Diagnostic::Ptr(new Hydrology_hydroinput(this))},
    {"wallmelt",   Diagnostic::Ptr(new Hydrology_wallmelt(this))},
  };
  return result;
}

*/

//! Copies the W variable, the modeled transportable water layer thickness.
void HydrologyMod::subglacial_water_thickness(IceModelVec2S &result) const {
  m_log->message(2,"* entering subglacial_water_thickness ...\n");
  result.set(0.0);
}



// need to add max_timestep_impl()

MaxTimestep HydrologyMod::max_timestep_impl(double t) const {
  (void) t;
  m_log->message(2,"* entering max_timestep_impl ...\n");
  return MaxTimestep("HydrologyMod hydrology");
}






// find basal pressure gradient
void HydrologyMod::pressure_gradient(IceModelVec2V &result, IceModelVec2S &result_mag, IceModelVec2S &result_angle) {

  m_log->message(2,"* entering pressure_gradient ...\n");

  const double
    dx = m_grid->dx(),
    dy = m_grid->dy();



  IceModelVec::AccessList list{&m_Pover_ghosts};

  Hydrology::overburden_pressure(m_Pover_ghosts); // grab the overburden pressure
  m_Pover_ghosts.update_ghosts();

  list.add(result);
  list.add(result_mag);
  list.add(result_angle);


  // third order finite difference method for calculating gradient, should give good results (Skidmore, 1989)
  for (Points p(*m_grid); p; p.next()) {										
    const int i = p.i(), j = p.j();


    result(i,j).u = (m_Pover_ghosts(i+1,j+1) + 2.0 * m_Pover_ghosts(i+1,j) + m_Pover_ghosts(i+1,j-1)) -
			    (m_Pover_ghosts(i-1,j+1) + 2.0 * m_Pover_ghosts(i-1,j) + m_Pover_ghosts(i-1,j-1)) / (8.0 * dx);

    result(i,j).v = (m_Pover_ghosts(i+1,j+1) + 2.0 * m_Pover_ghosts(i,j+1) + m_Pover_ghosts(i-1,j+1)) -
			    (m_Pover_ghosts(i+1,j-1) + 2.0 * m_Pover_ghosts(i,j-1) + m_Pover_ghosts(i-1,j-1)) / (8.0 * dy);


  }

  m_log->message(2,"* calculated gradient, now removing holes ...\n");

  // check to see if there are adjacent grid cells with gradient vector components that oppose each other. If so, the one
  // that has a larger magnitude wins, and the smaller one gets set to zero. This is why the ghosts are reset.
  // Obviously, not exactly mass conserving, basically any water that goes into one of these grid cells will become a black
  // hole for water using my "drain completely" scheme.
  // After, the magnitude of the gradient is calculated. Also, the direction of water flow is calculated.

  result.update_ghosts();

  for (Points p(*m_grid); p; p.next()) {										
    const int i = p.i(), j = p.j();
    // remembering that the gradient will point in the opposite direction of flow
    if( result(i,j).u < 0 && result(i+1,j).u > 0) {
      if(abs(result(i,j).u) < abs(result(i+1,j).u)) {
        result(i,j).u = 0;
      }
    }

    if( result(i,j).u > 0 && result(i-1,j).u < 0) {
      if(abs(result(i,j).u) < abs(result(i+1,j).u)) {
        result(i,j).u = 0;
      }
    }


    if( result(i,j).v < 0 && result(i,j+1).v > 0) {
      if(abs(result(i,j).v) < abs(result(i,j+1).v)) {
        result(i,j).v = 0;
      }
    }

    if( result(i,j).v > 0 && result(i,j-1).v < 0) {
      if(abs(result(i,j).v) < abs(result(i,j+1).v)) {
        result(i,j).v = 0;
      }
    }

    result_mag(i,j) = sqrt(pow(result(i,j).u,2) + (result(i,j).v,2));
    if(result_mag(i,j) > 0) {
      result_angle(i,j) = atan2 (-result(i,j).v,-result(i,j).u);
    }
    else {
      result_angle(i,j) = 0;
    }

  }

  result.update_ghosts();
  result_angle.update_ghosts();

  m_log->message(2,"* finished pressure_gradient ...\n");

} // end pressure_gradient()

//! Returns the (trivial) overburden pressure as the pressure of the non-existent transportable water, because this is the least harmful output if this is misused.
void HydrologyMod::subglacial_water_pressure(IceModelVec2S &result) const {
  m_log->message(2,"* entering subglacial_water_pressure ...\n");
  result.copy_from(m_Pover_ghosts);
}

// calculates the maximum drainage rate of the till based on the sediment type and gradient
void HydrologyMod::till_drainage(IceModelVec2S &result, double dt) {
  m_log->message(2,"* entering till_drainage ...\n");
  double m_dt = dt;
  IceModelVec::AccessList list;
  list.add(m_till_permeability);
  list.add(m_gradient_magnitude);
  list.add(m_Wtil);
  list.add(m_pressure_gradient);
  list.add(m_theta);
  list.add(m_fraction_till);
  
  const double water_viscosity = m_config->get_double("hydrology.water_viscosity");

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // Darcy's Law, combined with fraction sediments
    result(i,j) = - m_till_permeability(i,j) / water_viscosity * m_gradient_magnitude(i,j)* m_fraction_till(i,j);  // flux out of the cell from sediments

   // of course, it is possible that the rate of discharge cannot be so much that it completely drains the till layer before the end of the next time step
   // I'm setting it so that the rate can only be at maximum complete drainage during the time step

   double drainage_thickness = result(i,j) * m_dt;

   if(m_Wtil(i,j)*m_fraction_till(i,j) < drainage_thickness) {
     result(i,j) = m_Wtil(i,j)*m_fraction_till(i,j) / m_dt;
   }


  }  // end for

  m_log->message(2,"* finished calculating basic input, now find how much drains out ...\n");

  result.update_ghosts();

  // next find the amount of water flow into each cell from adjacent cells with the updated ghosts


  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // sediment water flux into the cell due to adjacent sediments
    if(m_pressure_gradient(i+1,j).u > 0) { // water will flow in from the right cell
      result(i,j) -= result(i+1,j)  * cos(m_theta(i+1,j));
    }

    if(m_pressure_gradient(i-1,j).u < 0) { // water will flow in from the left cell
      result(i,j) -= result(i-1,j)* cos(m_theta(i-1,j));
    }

    if(m_pressure_gradient(i,j+1).v > 0) { // water will flow in from the cell above
      result(i,j) -= result(i,j+1)  * sin(m_theta(i,j+1));
    }
    if(m_pressure_gradient(i,j-1).v < 0) { // water will flow in from the cell below
      result(i,j) -= result(i,j-1) * sin(m_theta(i,j-1));
    }

  }
  m_log->message(2,"* finished till_drainage ...\n");
}

/*
Rothlisberger channel flow
finds the fraction of the base that has tunnels, given the excess water
based on the Manning formula (Cuffey and Paterson 2010 - equation 6.17)

The Manning roughness coefficient is a number between 10^-2 and 10^-1 s m^(-1/3), according to Cuffey and Paterson 2010.
According to Nye (1976), this value is controlled by how straight the tunnel is, with values close to 10^-2
for a straight tunnel on a smooth bed, and values of 10^-1 for a meandering tunnel. Nye (1976) calculated a value
of 0.12 for the drainage of subglacial lakes under the Vantajokull in Iceland, which is at the high end of this range.
This value corresponds to the value of "f" in the Arnold and Sharp papers (1992 and 2002) of 700 m^(-8/3) kg.
As such, the value of the roughness should probably be a function of bed properties, but it is not clear
what this should be. For the moment, I am just going to use the value used by Nye, because this seems pretty standard.

*/

void HydrologyMod::tunnels(IceModelVec2S &result) {

  m_log->message(2,"* entering tunnels ...\n");

  const double
    radius = m_config->get_double("hydrology.tunnel_radius"),
    roughness_coefficient = m_config->get_double("hydrology.roughness_coefficient"),
    pi = 3.14159265358979323846;

  double cross_section_tunnel = pi * pow(radius,2.0);
  double Rh = radius/2.0;



  IceModelVec::AccessList list{&m_fraction_channel,&m_excess_water,&result};
  //ParallelSection loop(m_grid->com);

  result.set(0.0);

  m_log->message(2,"* finished tunnels ...\n");

}

} // end namespace hydrology
} // end namespace pism
