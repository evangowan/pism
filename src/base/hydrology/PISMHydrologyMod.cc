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

#include <iostream>
#include <algorithm>    // std::min, max

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
//  m_tillwat_flux.metadata().set_double("valid_min", 0.0);
//  m_grid->variables().add(m_tillwat_flux);

  m_excess_water.create(m_grid, "excess_water", WITH_GHOSTS); 
  m_excess_water.set_attrs("model_state",
                       "extra water flux from water not in till",
                       "m s-1", "");
  m_excess_water.metadata().set_double("valid_min", 0.0);



  m_excess_water_playground.create(m_grid, "excess_water_playground", WITHOUT_GHOSTS); 
  m_excess_water_playground.set_attrs("internal",
                       "extra water flux from water not in till",
                       "m s-1", "");
  m_excess_water_playground.metadata().set_double("valid_min", 0.0);


  m_excess_water_removed.create(m_grid, "excess_water_removed", WITHOUT_GHOSTS); 
  m_excess_water_removed.set_attrs("internal",
                       "extra water flux from water check",
                       "m s-1", "");
  m_excess_water_removed.metadata().set_double("valid_min", 0.0);



  bottom_left.create(m_grid, "translate_bottom_left", WITH_GHOSTS,1); 
  bottom_left.set_attrs("internal",
                       "translation of bottom left corner of the grid cell",
                       "1", "");


  bottom_right.create(m_grid, "translate_bottom_right", WITH_GHOSTS,1); 
  bottom_right.set_attrs("internal",
                       "translation of bottom right corner of the grid cell",
                       "1", "");


  top_left.create(m_grid, "translate_top_left", WITH_GHOSTS,1); 
  top_left.set_attrs("internal",
                       "translation of top left corner of the grid cell",
                       "1", "");


  top_right.create(m_grid, "translate_top_right", WITH_GHOSTS,1); 
  top_right.set_attrs("internal",
                       "translation of top right corner of the grid cell",
                       "1", "");

  quad_area.create(m_grid, "quad_area", WITH_GHOSTS,1); 
  quad_area.set_attrs("internal",
                       "area of the quadralateral from the translated grid cell",
                       "1", "");


//  m_grid->variables().add(m_excess_water);

  m_theta.create(m_grid, "water_flow_direction", WITH_GHOSTS,1);
  m_theta.set_attrs("model_state",
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
    m_log->message(2, "* INIT_RESTART sediment thickness cover...\n");
    m_fraction_till.read(opts.filename, opts.record);
  } else if (opts.type == INIT_BOOTSTRAP) {
    m_log->message(2, "* INIT_BOOTSTRAP sediment thickness cover...\n");
    m_fraction_till.regrid(opts.filename, OPTIONAL, default_fraction_till);
  } else {
    m_log->message(2, "* setting sediment thickness cover to default_fraction_till...\n");
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
  double input_rate;
   //   m_log->message(2, "const_bmelt: %16.14f  \n",const_bmelt);
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
    input_rate=result(i,j);
    if(mask.icy(i,j)){
  //  m_log->message(2, " %5i %5i %16.14f %16.14f %16.14f %16.14f  \n", i, j, sqrt(pow(double(i-30),2)+pow(double(j-30),2)), m_bmelt_local(i,j), m_excess_water(i,j), input_rate );
   }
  }

  m_log->message(2, "* finished get_input_rate\n");
}


// Update the basal hydrology
// taken from NullTransportHydrology, but modified for my own purposes
void HydrologyMod::update_impl(double t, double dt) {


  m_log->message(2,
             "* calculating hydrology ...\n");

  double dx = m_grid->dx();
  double dy = m_grid->dy();

  if(fabs(dx - dy) > 0.00001) {

  m_log->message(2,
             "* warning, dx and dy are different, this means that the results of HydrologyMod are likely nonsense ...\n");

  }

  double m_t = t;
  double m_dt = dt;

  IceModelVec::AccessList list{&m_tillwat_flux, &m_pressure_gradient, &m_total_input,&m_theta,&m_gradient_magnitude};

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
  double flux_before;
  double flux_ratio;
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    flux_before= m_tillwat_flux(i,j);
 //     m_log->message(2, "* common vars...\n * dt: %10.8f \n * m_fraction_till: %7.2f \n * m_total_input %10.8f \n * m_tillwat_flux %10.8f \n", m_dt, m_fraction_till(i,j), m_total_input(i, j), m_tillwat_flux(i,j));

    if (mask.ocean(i,j) ) {

      // the other models say that till in the ocean will be have zero water layer thickness, of course it should be fully saturated
      m_Wtil(i, j) = tillwat_max;
      m_excess_water(i,j) = 0;
    } else if( mask.ice_free(i,j)) {

      // land without ice will be dry, will need to change it when there are lakes added
      m_Wtil(i,j) = 0.0;
      m_excess_water(i,j) = 0;

    } else {

      // till is filled if the amount of water entering the cell is greater than the till water flux

      if(m_tillwat_flux(i,j) > 0.0) {
        flux_ratio = m_total_input(i, j) / (m_tillwat_flux(i,j) * m_fraction_till(i,j));
      } else {
       if(m_total_input(i, j) > 0.0) {
        flux_ratio = 1.0;
       } else {
        flux_ratio = 0.0;
       }
      }

//      if(m_gradient_magnitude(i,j) < 1) { // any place where the gradient is zero, assume it is frozen to base with no till water. Probably later make this dependent on basal temperature
//       flux_ratio = 0.0;
//      }
      
      if(flux_ratio >= 1.0) {
        m_Wtil(i,j) = tillwat_max;
      } else {

        m_Wtil(i,j) = tillwat_max * flux_ratio;
      } // end if


      m_excess_water(i,j) = m_total_input(i, j);

 //     m_log->message(2, "%5i %5i %15.10f %15.10f %15.10f \n", i, j, m_total_input(i, j), m_excess_water(i,j), flux_ratio  );

    } // end if

/*
    } else if (m_total_input(i, j)/m_fraction_till(i,j) < m_tillwat_flux(i,j)) {  

      // the till water flux cannot be greater than the input

      m_Wtil(i,j) = m_dt * m_total_input(i, j)/m_fraction_till(i,j);
      

      m_tillwat_flux(i,j) = m_total_input(i, j)/m_fraction_till(i,j);

    } else {
      
      m_Wtil(i,j) = m_dt * m_total_input(i, j)/m_fraction_till(i,j);

    }

*/


//    m_log->message(2, "* m_Wtil: %5i %5i %10.8f %15.10f %15.10f %15.10f \n ", i, j, m_Wtil(i,j), m_total_input(i, j)*m_dt, -m_tillwat_flux(i,j)*m_dt,m_gradient_magnitude(i,j)  );

/*

      if (m_Wtil(i,j) > tillwat_max) { // add the extra water to m_excess_water
         if(m_gradient_magnitude(i,j) > 0.0) { // but only if the magnitude of the gradient is greater than 0, right not zero gradient means black hole

           m_excess_water(i,j) = (m_total_input(i, j) - tillwat_max * m_fraction_till(i,j) / m_dt);
         } else{
           m_excess_water(i,j) = 0;
         }
         m_Wtil(i,j) = tillwat_max;

      } else { // no excess water and the till is not full

       m_excess_water(i,j) = 0;


      }
     if(mask.icy(i,j)){
 //    m_log->message(2, "* m_wtil: %5i %5i %10.8f %15.10f %15.10f %15.10f \n ", i, j, m_Wtil(i,j), m_excess_water(i,j), m_total_input(i, j)*m_dt, -m_tillwat_flux(i,j)*m_dt );

    m_log->message(2, "* m_wtil: %5i %5i %15.10f %15.10f %15.10f \n ", i, j,m_total_input(i, j), flux_before, m_tillwat_flux(i,j) );
    }
*/

  } // end for



// finally, update the excess water to have it go into adjacent cells


  list.add(m_excess_water_playground);
  list.add(m_excess_water_removed);

  m_excess_water.update_ghosts();
  m_theta.update_ghosts();

//  m_excess_water_playground.copy_from(m_excess_water);
  m_excess_water_playground.set(0.0);
  m_excess_water_removed.set(0.0);
  double seconds_in_year = 365*24*3600;
   m_log->message(2,"* till_drainage before: \n");
  // next find the amount of water flow into each cell from adjacent cells with the updated ghosts
  double drainage, drainage_before;
  double left, right, up, down, topleft, topright, bottomleft, bottomright, center;
  double res_left, res_right, res_up, res_down;
  double right_fraction, top_fraction;
  double area, total_area, distance;
  double  pi = 3.14159265358979323846;

  double translate_center[3][3][2], transformation[2][2][2];

  int x_counter, y_counter;

  double quadrilateral[4][2];

  double a_x, a_y, b_x, b_y;

  const IceModelVec2S &thk = *m_grid->variables().get_2d_scalar("thk");
  list.add(thk);
   m_log->message(2,"* before cell translation: \n");

  list.add(bottom_left);
  list.add(top_left);
  list.add(top_right);
  list.add(bottom_right);
  list.add(quad_area);
  quad_area.set(1.0); // should be a unit cell outside of the icy area

  // calculate the translation of the grid cell's corners, to determine water content to surrounding cells
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
//   m_log->message(2,"* before mask %5i %5i\n", i, j);
   if(mask.icy(i,j) ) {

    // first find the translation of the center of all the cells
//    m_log->message(2,"* calculating the translation %5i %5i\n", i, j);

//   m_log->message(2,"* check_sin_and_cos %5i %5i\n", i, j);
    for(x_counter=0; x_counter<3; x_counter++ ){
      for(y_counter=0; y_counter<3; y_counter++ ){ 



       if(m_gradient_magnitude(i+x_counter-1,j+y_counter-1) > 0.0 && mask.icy(i+x_counter-1,j+y_counter-1)) {
        translate_center[x_counter][y_counter][0] = cos(m_theta(i+x_counter-1,j+y_counter-1)); // translate x
        translate_center[x_counter][y_counter][1] = sin(m_theta(i+x_counter-1,j+y_counter-1)); // translate y
       } else{
        translate_center[x_counter][y_counter][0] = 0; // translate x
        translate_center[x_counter][y_counter][1] = 0; // translate y  
       }

 //     m_log->message(2,"translate %5i %5i %15.10f %15.10f %15.10f\n",i+x_counter-1,j+y_counter-1, m_gradient_magnitude(i+x_counter-1,j+y_counter-1),  translate_center[x_counter][y_counter][0],translate_center[x_counter][y_counter][1]);

 //      m_log->message(2,"%5i %5i %15.10f %15.10f\n", i+x_counter-1,j+y_counter-1, translate_center[x_counter][y_counter][0],translate_center[x_counter][y_counter][1]);

      }  // end for
    } // end for

//    m_log->message(2,"* calculating the translation for each corner %5i %5i\n", i, j);
    // bottom left

//    m_log->message(2,"* bottom left %5i %5i\n", i, j);


    transformation[0][0][0] = translate_center[0][0][0];
    transformation[0][0][1] = translate_center[0][0][1];
    transformation[1][0][0] = translate_center[1][0][0];
    transformation[1][0][1] = translate_center[1][0][1];
    transformation[1][1][0] = translate_center[1][1][0];
    transformation[1][1][1] = translate_center[1][1][1];
    transformation[0][1][0] = translate_center[0][1][0];
    transformation[0][1][1] = translate_center[0][1][1];
/*
       m_log->message(2,"translate %15.10f %15.10f\n", transformation[0][0][0],transformation[0][0][1]);
       m_log->message(2,"translate %15.10f %15.10f\n", transformation[1][0][0],transformation[1][0][1]);
       m_log->message(2,"translate %15.10f %15.10f\n", transformation[1][1][0],transformation[1][1][1]);
       m_log->message(2,"translate %15.10f %15.10f\n", transformation[0][1][0],transformation[0][1][1]);
*/
    // calculate the affine transformation

    projection_transformation(transformation,quadrilateral[0][0], quadrilateral[0][1]);
    quadrilateral[0][0] = quadrilateral[0][0]-1.0;
    quadrilateral[0][1] = quadrilateral[0][1]-1.0;
    bottom_left(i,j).u = quadrilateral[0][0];
    bottom_left(i,j).v = quadrilateral[0][1];


    // top left
//    m_log->message(2,"* top left %5i %5i\n", i, j);

    transformation[0][0][0] = translate_center[0][1][0];
    transformation[0][0][1] = translate_center[0][1][1];
    transformation[1][0][0] = translate_center[1][1][0];
    transformation[1][0][1] = translate_center[1][1][1];
    transformation[1][1][0] = translate_center[1][2][0];
    transformation[1][1][1] = translate_center[1][2][1];
    transformation[0][1][0] = translate_center[0][2][0];
    transformation[0][1][1] = translate_center[0][2][1];

    projection_transformation(transformation,quadrilateral[1][0], quadrilateral[1][1]);
    quadrilateral[1][0] = quadrilateral[1][0]-1.0;
    top_left(i,j).u = quadrilateral[1][0];
    top_left(i,j).v = quadrilateral[1][1];


    // top right
//    m_log->message(2,"* top right %5i %5i\n", i, j);


    transformation[0][0][0] = translate_center[1][1][0];
    transformation[0][0][1] = translate_center[1][1][1];
    transformation[1][0][0] = translate_center[2][1][0];
    transformation[1][0][1] = translate_center[2][1][1];
    transformation[1][1][0] = translate_center[2][2][0];
    transformation[1][1][1] = translate_center[2][2][1];
    transformation[0][1][0] = translate_center[1][2][0];
    transformation[0][1][1] = translate_center[1][2][1];


    projection_transformation(transformation,quadrilateral[2][0], quadrilateral[2][1]);

    top_right(i,j).u = quadrilateral[2][0];
    top_right(i,j).v = quadrilateral[2][1];



    // bottom right
//    m_log->message(2,"* bottom right %5i %5i\n", i, j);

    transformation[0][0][0] = translate_center[1][0][0];
    transformation[0][0][1] = translate_center[1][0][1];
    transformation[1][0][0] = translate_center[2][0][0];
    transformation[1][0][1] = translate_center[2][0][1];
    transformation[1][1][0] = translate_center[2][1][0];
    transformation[1][1][1] = translate_center[2][1][1];
    transformation[0][1][0] = translate_center[1][1][0];
    transformation[0][1][1] = translate_center[1][1][1];
/*
      m_log->message(2,"translate %15.10f %15.10f %15.10f %15.10f\n",0.0, 0.0, transformation[0][0][0]+double(i),transformation[0][0][1]+double(j)-1.0);
      m_log->message(2,"translate %15.10f %15.10f %15.10f %15.10f\n",1.0, 0.0, transformation[1][0][0]+double(i)+1.0 ,transformation[1][0][1]+double(j)-1.0);
      m_log->message(2,"translate %15.10f %15.10f %15.10f %15.10f\n",1.0, 1.0, transformation[1][1][0]+double(i)+1.0 ,transformation[1][1][1]+double(j));
      m_log->message(2,"translate %15.10f %15.10f %15.10f %15.10f\n",0.0, 1.0, transformation[0][1][0]+double(i),transformation[0][1][1]+double(j));
*/
    projection_transformation(transformation,quadrilateral[3][0], quadrilateral[3][1]);

    quadrilateral[3][1] = quadrilateral[3][1]-1.0;
    bottom_right(i,j).u = quadrilateral[3][0];
    bottom_right(i,j).v = quadrilateral[3][1];



//    m_log->message(2,"%5i %5i %15.10f %15.10f\n", i, j, bottom_left(i,j).u, bottom_left(i,j).v);
//    m_log->message(2,"%5i %5i %15.10f %15.10f\n", i, j, top_left(i,j).u, top_left(i,j).v);
//    m_log->message(2,"%5i %5i %15.10f %15.10f\n", i, j, top_right(i,j).u, top_right(i,j).v);  
 //   m_log->message(2,"%5i %5i %15.10f %15.10f\n", i, j, bottom_right(i,j).u, bottom_right(i,j).v);  


    // TODO a check to make sure the translated quadralateral is regular? Probably shouldn't happen if the gradient is calculated correctly.


    // find the area of the quadralateral

    quad_area(i,j) = find_quad_area(quadrilateral);
/*
         m_log->message(2,"> %5i %5i %15.10f\n", i, j, quad_area(i,j));
         m_log->message(2,"%15.10f %15.10f\n", double(i) + bottom_left(i,j).u, double(j) + bottom_left(i,j).v);
         m_log->message(2,"%15.10f %15.10f\n", double(i) + top_left(i,j).u, double(j) + top_left(i,j).v);
         m_log->message(2,"%15.10f %15.10f\n", double(i) + top_right(i,j).u, double(j) + top_right(i,j).v);
         m_log->message(2,"%15.10f %15.10f\n", double(i) + bottom_right(i,j).u, double(j) + bottom_right(i,j).v);
*/
   } else {

//   m_log->message(2,"* not icy %5i %5i\n", i, j);
    // bottom right

    bottom_right(i,j).u = 0.5;
    bottom_right(i,j).v = -0.5;
//   m_log->message(2,"* finished bottom_right %5i %5i\n", i, j);

    // bottom left

    bottom_left(i,j).u = -0.5;
    bottom_left(i,j).v = -0.5;


    // top right

    top_right(i,j).u = 0.5;
    top_right(i,j).v = 0.5;


    // top left

    top_left(i,j).u = -0.5;
    top_left(i,j).v = 0.5;

    quad_area(i,j) = 1.0;
  // m_log->message(2,"* finished not icy %5i %5i\n", i, j);
   } // end if (hope this is in the right spot
 
  } // end for
   m_log->message(2,"* after cell translation: \n");


  // now that we know where the water spreads to, find how much goes into each cell

  // first update the ghosts

  top_left.update_ghosts();
  top_right.update_ghosts();
  bottom_left.update_ghosts();
  bottom_right.update_ghosts();
  quad_area.update_ghosts();
  
  double reference_cell[4][2], compare_cell[4][2];

  reference_cell[0][0] = -0.5; // bottom left
  reference_cell[0][1] = -0.5;
  reference_cell[1][0] = -0.5; // top left
  reference_cell[1][1] =  0.5;
  reference_cell[2][0] =  0.5; // top right
  reference_cell[2][1] =  0.5;
  reference_cell[3][0] =  0.5; // bottom right
  reference_cell[3][1] = -0.5;
  double sum_excess_water = 0.0, sum_excess_water_playground = 0.0;
  m_excess_water_playground.set(0.0);
   m_log->message(2,"* before calculate_water: \n");
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    sum_excess_water += m_excess_water(i,j);
//   if(mask.icy(i,j) ) {
//      m_log->message(2,"* icy: %5i %5i\n", i, j);
    // make a call to the formula 

    for(x_counter=0; x_counter<3; x_counter++ ){
      for(y_counter=0; y_counter<3; y_counter++ ){ 

       int x_index = i+x_counter-1;
       int y_index = j+y_counter-1;

       if(mask.icy(x_index,y_index)) {
        compare_cell[0][0] = double(x_counter-1) + bottom_left(x_index,y_index).u; // bottom left
        compare_cell[0][1] = double(y_counter-1) + bottom_left(x_index,y_index).v; 
        compare_cell[1][0] = double(x_counter-1) + top_left(x_index,y_index).u; // top left
        compare_cell[1][1] = double(y_counter-1) + top_left(x_index,y_index).v;
        compare_cell[2][0] = double(x_counter-1) + top_right(x_index,y_index).u; // top right
        compare_cell[2][1] = double(y_counter-1) + top_right(x_index,y_index).v; 
        compare_cell[3][0] = double(x_counter-1) + bottom_right(x_index,y_index).u; // bottom right
        compare_cell[3][1] = double(y_counter-1) + bottom_right(x_index,y_index).v;

/*
        if(x_counter == 1 && y_counter == 1) {

         m_log->message(2,">\n");
         m_log->message(2,"%15.10f %15.10f\n", double(i) + compare_cell[0][0], double(j) + compare_cell[0][1]);
         m_log->message(2,"%15.10f %15.10f\n", double(i) + compare_cell[1][0], double(j) + compare_cell[1][1]);
         m_log->message(2,"%15.10f %15.10f\n", double(i) + compare_cell[2][0], double(j) + compare_cell[2][1]);
         m_log->message(2,"%15.10f %15.10f\n", double(i) + compare_cell[3][0], double(j) + compare_cell[3][1]);
       }
*/
/*
        m_log->message(2,"* calling calculate_water: %5i %5i %5i %5i %15.10f %15.10f\n", i, j,i + (x_counter-1), j + (y_counter-1), m_excess_water(i + (x_counter-1), j + (y_counter-1)), quad_area(i + (x_counter-1), j + (y_counter-1)));
        m_log->message(2,"%15.10f %15.10f %15.10f %15.10f\n", compare_cell[0][0],  compare_cell[0][1], bottom_left(i,j).u, bottom_left(i,j).v);
        m_log->message(2,"%15.10f %15.10f %15.10f %15.10f\n", compare_cell[1][0],  compare_cell[1][1], top_left(i,j).u, top_left(i,j).v);
        m_log->message(2,"%15.10f %15.10f %15.10f %15.10f\n", compare_cell[2][0],  compare_cell[2][1], top_right(i,j).u, top_right(i,j).v);
        m_log->message(2,"%15.10f %15.10f %15.10f %15.10f\n", compare_cell[3][0],  compare_cell[3][1], bottom_right(i,j).u, bottom_right(i,j).v);
*/

    //    m_log->message(2, "+ %5i %5i\n", i, j); 
        m_excess_water_playground(i,j) += calculate_water(reference_cell, compare_cell,false) * m_excess_water(i + (x_counter-1), j + (y_counter-1)) / quad_area(i + (x_counter-1), j + (y_counter-1));
/*
      if(i  == 17 && j  == 30 && t > 2060.0 * 365.0 * 24.0 * 3600.0 ){

      m_log->message(2,"+ %5i %5i %5i %5i %15.10f %15.10f %15.10f %15.10f\n", i, j, x_counter-1, y_counter-1,calculate_water(reference_cell, compare_cell,false), m_excess_water(i + (x_counter-1), j + (y_counter-1)), quad_area(i + (x_counter-1), j + (y_counter-1)), calculate_water(reference_cell, compare_cell,false) * m_excess_water(i + (x_counter-1), j + (y_counter-1)) / quad_area(i + (x_counter-1), j + (y_counter-1)));
     double dummy = calculate_water(reference_cell, compare_cell,true);
     }
*/

       } // end if
      } // end for
    } // end for




/*
      if(m_excess_water_playground(i,j) > 0.0){
         m_log->message(2,"> -Z%15.10f\n",  m_excess_water_playground(i,j));
         m_log->message(2,"%15.10f %15.10f\n", double(i)-0.5, double(j)-0.5 );
         m_log->message(2,"%15.10f %15.10f\n", double(i)-0.5, double(j)+0.5 );
         m_log->message(2,"%15.10f %15.10f\n", double(i)+0.5, double(j)+0.5 );
         m_log->message(2,"%15.10f %15.10f\n", double(i)+0.5, double(j)-0.5 );
         m_log->message(2,"%15.10f %15.10f\n", double(i)-0.5, double(j)-0.5 );
      }
*/
/*
     if (m_excess_water_playground(i,j) > 1e-8) { // fudge
       m_excess_water_playground(i,j) = 1e-8;
     }
*/

    sum_excess_water_playground += m_excess_water_playground(i,j);
 if(mask.icy(i,j) && t > 1800.0 * 365.0 * 24.0 * 3600.0){
 //m_log->message(2,"%5i %5i %15.10f %15.10f\n", i, j, m_excess_water(i , j)*1e8,  m_excess_water_playground(i,j)*1e8);
}

//   } // end if
  } // end for

  m_excess_water.copy_from(m_excess_water_playground);
  m_log->message(2,
             "* finished calculating hydrology ...\n");
  
  m_log->message(2,
             "* sums: ... %15.10f %15.10f\n",sum_excess_water,sum_excess_water_playground );

}


void HydrologyMod::projection_transformation(double transformation[2][2][2],double& x, double& y){

 // from Heckbert 1989 (Fundamentals of texture mapping and image warping, page 20)

 double epsilon = 0.0000001;

 // make it a unit square


 transformation[0][1][1] = transformation[0][1][1] + 1.0;
 transformation[1][1][0] = transformation[1][1][0] + 1.0;
 transformation[1][1][1] = transformation[1][1][1] + 1.0;
 transformation[1][0][0] = transformation[1][0][0] + 1.0;


 // delta

 


 double delta[3][2];
 
 for(int counter=0; counter<=1; counter++){

  delta[0][counter] = transformation[1][0][counter] - transformation[1][1][counter];
  delta[1][counter] = transformation[0][1][counter] - transformation[1][1][counter];
  delta[2][counter] = transformation[0][0][counter] - transformation[1][0][counter] + transformation[1][1][counter] - transformation[0][1][counter];
 }

 // create transformation matrix

 double a_matrix[3][3];

 double denom_left = delta[0][0] * delta[1][1];
 double denom_right = delta[0][1] * delta[1][0];

 if(std::fabs(denom_left - denom_right) > epsilon) { // should not cause numerical errors

  a_matrix[2][0] = (delta[2][0] * delta[1][1] - delta[2][1] * delta[1][0]) / (denom_left - denom_right); // g
  a_matrix[2][1] = (delta[0][0] * delta[2][1] - delta[0][1] * delta[2][0]) / (denom_left - denom_right); // h
  a_matrix[2][2] = 1.0; // i

  a_matrix[0][0] = transformation[1][0][0] - transformation[0][0][0] + a_matrix[2][0] * transformation[1][0][0]; // a
  a_matrix[0][1] = transformation[0][1][0] - transformation[0][0][0] + a_matrix[2][1] * transformation[0][1][0]; // b
  a_matrix[0][2] = transformation[0][0][0]; // c
  a_matrix[1][0] = transformation[1][0][1] - transformation[0][0][1] + a_matrix[2][0] * transformation[1][0][1]; // d
  a_matrix[1][1] = transformation[0][1][1] - transformation[0][0][1] + a_matrix[2][1] * transformation[0][1][1]; // e
 a_matrix[1][2] = transformation[0][0][1]; // f
 



  x = (0.5 * a_matrix[0][0] + 0.5 * a_matrix[0][1] + a_matrix[0][2]) / (0.5 * a_matrix[2][0] + 0.5 * a_matrix[2][1] + a_matrix[2][2]);
  y = (0.5 * a_matrix[1][0] + 0.5 * a_matrix[1][1] + a_matrix[1][2]) / (0.5 * a_matrix[2][0] + 0.5 * a_matrix[2][1] + a_matrix[2][2]);

 } else {

  // do bilinear mapping instead since there would be a divide by zero error with a projection

  double a = transformation[0][0][0] - transformation[1][0][0] - transformation[0][1][0] + transformation[1][1][0];
  double b = -transformation[0][0][0] + transformation[1][0][0];
  double c = -transformation[0][0][0] + transformation[0][1][0];
  double d = transformation[0][0][0];

  double e = transformation[0][0][1] - transformation[1][0][1] - transformation[0][1][1] + transformation[1][1][1];
  double f = -transformation[0][0][1] + transformation[1][0][1];
  double g = -transformation[0][0][0] + transformation[0][1][1];
  double h = transformation[0][0][1];

  x = 0.5 * 0.5 * a + 0.5 * b + 0.5 * c + d;
  y = 0.5 * 0.5 * e + 0.5 * f + 0.5 * g + h;

 }

}


double HydrologyMod::find_quad_area(double quadrilateral[4][2]) {

//  m_log->message(2,
//             "* calculating quad area ...\n");
  double p_x = quadrilateral[3][0] - quadrilateral[1][0];
  double p_y = quadrilateral[3][1] - quadrilateral[1][1];
  double q_x = quadrilateral[2][0] - quadrilateral[0][0];
  double q_y = quadrilateral[2][1] - quadrilateral[0][1]; 


  return fabs(p_x*q_y - p_y*q_x) * 0.5;

}

double HydrologyMod::calculate_water(double reference_cell[4][2], double compare_cell[4][2], bool printing) {

  bool corner_inside1[4], corner_inside2[4];

  int counter;

  int polygon_size = 4;

  bool all_inside = true;
  bool all_outside = true;

  bool all_inside2 = true;
  bool all_outside2 = true;

  double water;

   double epsilon = 0.0000001;

//  m_log->message(2, "* within_calculate_water ...\n");

  polygon_linked_list reference;
//  polygon_linked_list compare(2);
 //   m_log->message(2, "* linked lists created ...\n");
  int reference_node_count = 0;
  int compare_node_count = 0;


  node * reference_node2;
  node * compare_node2;
  node * reference_node;

  node * created_node;

  bool on_edge = true; // used in point_in_polygon to check if the point is on the edge of the reference polygon

//  m_log->message(2, "* establish polygon lists ...\n");
  // first add all the points in the reference cell
  for(counter = 0; counter < polygon_size; counter++) {

//     m_log->message(2, "* assigning points %5i \n", counter);
    // check if reference corner is inside the reference
    corner_inside1[counter] = point_in_polygon(compare_cell, polygon_size, reference_cell[counter][0], reference_cell[counter][1], on_edge);

    if(!corner_inside1[counter]) {
//       m_log->message(2, "* turning off all inside %5i \n", corner_inside1[counter]);
      all_inside = false;

    }
//    created_node = new node;
    reference.create_node(created_node);
    reference_node =  created_node;
    reference_node -> x = reference_cell[counter][0];
    reference_node -> y = reference_cell[counter][1];
    reference_node -> inside = corner_inside1[counter];

    reference.insertNode(reference_node,counter+1,0);
    reference_node_count++;

//      m_log->message(2, "* adding reference_node %5i %2i %p %p \n", counter, corner_inside1[counter], reference_node, reference_node -> next[0]); 
   /*
    node * compare_node = new node;
    compare_node -> x = compare_cell[counter][0];
    compare_node -> y = compare_cell[counter][1];
    compare_node -> shared_node_number = 0;
    compare_node -> inside = corner_inside2[counter];
    compare_node -> next[0] = NULL;
    compare_node -> next[2] = NULL;
    m_log->message(2, "* adding compare_node %5i \n", counter); 
    compare.insertNode(compare_node,counter+1);
    compare_node_count++;
  */
//     m_log->message(2, "* assigned points %5i \n", counter);

    if(counter == 0) {  // create a full polygon this way 

    //  created_node = new node;
      reference.create_node(created_node);
      reference_node =  created_node;
      reference_node -> x = reference_cell[counter][0];
      reference_node -> y = reference_cell[counter][1];
      reference_node -> inside = corner_inside1[counter];
      reference.insertNode(reference_node,2,0);
  //    reference_node_count++;
 //     m_log->message(2, "* adding reference_node %5i %2i %p %p \n", counter, corner_inside1[counter], reference_node, reference_node -> next[0]); 
  
    }





  }


//  reference_node = new node;

//  m_log->message(2, "* object after insert reference: ... \n");
//  reference.print_polygon(0);

 //    m_log->message(2, "* adding the compare object **********************************:\n");


  // next add the compare cell

  for(counter = 0; counter < polygon_size; counter++) {
   node * compare_node;

//   m_log->message(2, "* compare assigning points %5i \n", counter);

   corner_inside2[counter] = point_in_polygon(reference_cell, polygon_size, compare_cell[counter][0], compare_cell[counter][1], on_edge);

 //    m_log->message(2, "* after point_in_polygon %5i \n", corner_inside2[counter]);

    if(corner_inside1[counter] || corner_inside2[counter]) {
 //     m_log->message(2, "* all_outside should be false %5i \n", corner_inside2[counter]);
      all_outside = false;
    }

    if(!corner_inside2[counter]) {
//      m_log->message(2, "* all_inside2 should be false %5i \n", corner_inside2[counter]);
      all_inside2 = false;

    }

//     m_log->message(2, "* testing the reference object before creating a new compare node:\n");
 //   reference.findNode(reference_node,1);
//     m_log->message(2, "* before add node: %p %15.10f \n", reference_node -> next[0], reference_node -> x);

/*
    m_log->message(2, "* before creating node: \n");
    bool dummy = reference.findNode(reference_node,2,0);
    dummy = reference.findNode(compare_node,2,1);

    created_node = NULL;
    m_log->message(2, "* address of reference and compare nodes: %p %p %p\n", reference_node, compare_node, created_node);
//    created_node = new node;
 //   created_node = new node;
 //   compare_node = created_node;
  //  compare_node = new node;
    reference.create_node(created_node);
    reference.create_node(compare_node);
    m_log->message(2, "* address of reference and compare nodes: %p %p %p\n", reference_node, compare_node, created_node);
    dummy = reference.findNode(compare_node,2,1);
    m_log->message(2, "* done checking, adding: %5i %15.10f\n", counter, compare_cell[counter][0] );
    m_log->message(2, "* check pointer to x created: \n");
    created_node -> x = compare_cell[counter][0];
    m_log->message(2, "* check pointer to x compare:\n" );
    compare_node -> x = compare_cell[counter][0];
    compare_node -> y = compare_cell[counter][1];
    m_log->message(2, "* adding pointer information:\n");
    compare_node -> shared_node_number = 0;
    compare_node -> inside = corner_inside2[counter];
    compare_node -> next[0] = NULL;
    compare_node -> next[1] = NULL;
    compare_node -> next[2] = NULL;
*/


    // call function that will replace the node with the equivalent reference node if they are at the same x,y point.

  //  reference.findNode(reference_node,1);


    // test before function

//    m_log->message(2, "* checking reference: %p \n", reference); 

//    bool found_shared = find_shared_node(compare_node, reference);

//     bool dummy = reference.findNode(reference_node,1,0);



    int check_counter = 0;


    bool continue_search = true;

    bool found_shared = false;
    while (continue_search) {

     check_counter++;

  //   m_log->message(2, "* attempting to find reference node: %5i %p \n", check_counter, reference); 
     continue_search = reference.findNode(reference_node,check_counter,0);

     if(continue_search && std::fabs(compare_cell[counter][0] - reference_node -> x) < epsilon && std::fabs(compare_cell[counter][1] - reference_node -> y) < epsilon) {

      // compare_node = reference_node;
       found_shared = true;
       continue_search = false;
      }

    }

    

    if(found_shared) {
     corner_inside2[counter] =  true;
     reference_node -> inside = true;
 //    m_log->message(2, "* found shared \n");
     reference.insertNode(reference_node,counter+1,1);
     compare_node_count++;

    } else {
     reference.create_node(compare_node);
     compare_node -> x = compare_cell[counter][0];
     compare_node -> y = compare_cell[counter][1];
     compare_node -> inside = corner_inside2[counter];

   //  m_log->message(2, "* found non-overlapping %p %p %15.10f %15.10f \n", compare_node -> next[0], compare_node -> next[1], compare_node -> x, compare_node -> y);
     reference.insertNode(compare_node,counter+1,1);
     compare_node_count++;

     // check if node also is on the reference polygon
     if( corner_inside2[counter]) {

      
      int reference_count = 1;
      bool before_done = true;

      while(before_done) {
       before_done = reference.findNode(reference_node,reference_count,0);
       before_done = reference.findNode(reference_node2,reference_count+1,0);
       if(before_done) {
  //      m_log->message(2, "* checking: %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f \n", reference_node -> x, reference_node -> y, reference_node2 -> x, reference_node2 -> y, compare_node -> x, compare_node -> y  );
        if(std::fabs(reference_node -> x - reference_node2 -> x) < epsilon) { // vertical line
         if(std::fabs(reference_node -> x - compare_node -> x) < epsilon) { // falls on the line
          if (std::fabs(reference_node -> y - compare_node -> y) > epsilon && std::fabs(reference_node2 -> y - compare_node -> y) > epsilon) { // should fall between the lines
           reference.insertNode(compare_node,reference_count+1,0);
          // m_log->message(2, "* found overlap:\n");
           before_done = false;
          }
         }

         } else if(std::min(reference_node -> x, reference_node2 -> x) < compare_node -> x && std::max(reference_node -> x, reference_node2 -> x) > compare_node -> x) { // could overlap

          double slope = (reference_node -> y - reference_node2 -> y) /  (reference_node -> x - reference_node2 -> x);
          double intercept = reference_node -> y - slope * reference_node -> x;
          double temp_y = slope * compare_node -> x + intercept;
          if(std::fabs(compare_node -> y - temp_y) < epsilon ) { // sould overlap
           reference.insertNode(compare_node,reference_count+1,0);
         //  m_log->message(2, "* found overlap:\n");
           before_done = false;
          }
         }
        }
       
       reference_count++;
      }
     }
    }



    if(counter == 0) {  // create a full polygon this way 

  //    created_node = new node;
      reference.create_node(created_node);
      compare_node = created_node;
      compare_node -> x = compare_cell[counter][0];
      compare_node -> y = compare_cell[counter][1];
      compare_node -> inside = corner_inside2[counter];
      reference.insertNode(compare_node,2,1);
      compare_node_count++;
 //    m_log->message(2, "* added last node %p %p %15.10f %15.10f \n", compare_node -> next[0], compare_node -> next[1], compare_node -> x, compare_node -> y);
 //   dummy = reference.findNode(reference_node,1,0);
 //   dummy = reference.findNode(reference_node,1,1);
    }

  }


//  m_log->message(2, "* object after insert compare: ... \n");
//  reference.print_polygon(0);
//  reference.print_polygon(1);

  node * compare_node;

//  compare_node = new node;
  // find the crossover points
 
  // go back to the head
/*
  m_log->message(2, "* find the crossover points ...\n");
     m_log->message(2, "* found the first node ...\n");
   // this should work assuming there isn't multiple overlapping polygons
     m_log->message(2, "* reference ...\n");
  reference.print_polygon(0);
     m_log->message(2, "* compare ...\n");
  reference.print_polygon(1);
*/
  bool not_finished = true;
  bool not_finished2;
  bool is_crossover;

  int reference_point = 1;
  int compare_point;

  int shared_node_number = 0;

  node * crossover_node;

  while (not_finished) {

   not_finished = reference.findNode(reference_node,reference_point,0);
 //  m_log->message(2, "* proceed through the while loop %d ...\n", not_finished);
 //  if (not_finished) {

    not_finished = reference.findNode(reference_node2,reference_point+1,0);

    if (not_finished) {

     // This should have found the reference line segement, now find crossover points with the compare polygon

     not_finished2 = true;
     compare_point = 1;
     while(not_finished2) {
 
      not_finished2 = reference.findNode(compare_node,compare_point,1);

      //if(not_finished2) {

       not_finished2 = reference.findNode(compare_node2,compare_point+1,1);

       if(not_finished2) {

        // now have two line segments, find out if they intersect
    //     m_log->message(2, "* do the lines overlap? ...\n");
   //      m_log->message(2, "* %15.10f %15.10f %15.10f %15.10f \n", reference_node -> x, reference_node -> y, compare_node -> x, compare_node -> y);
   //      m_log->message(2, "* %15.10f %15.10f %15.10f %15.10f \n", reference_node2 -> x, reference_node2 -> y, compare_node2 -> x, compare_node2 -> y);
     //    created_node = new node;
      //   reference.create_node(created_node);
        // crossover_node = created_node;
         double x, y;
         is_crossover = find_crossover(reference_node, reference_node2, compare_node, compare_node2, x, y);

      if(printing){
         m_log->message(2, "* %15.10f %15.10f %15.10f %15.10f  %15.10f %15.10f %5i \n", reference_node -> x, reference_node -> y, compare_node -> x, compare_node -> y, x, y, is_crossover);
         m_log->message(2, "* %15.10f %15.10f %15.10f %15.10f \n", reference_node2 -> x, reference_node2 -> y, compare_node2 -> x, compare_node2 -> y);
      }

         if(is_crossover) {

           shared_node_number++;
           // first check that they aren't shared nodes
           if(std::fabs(x - reference_node -> x) < epsilon && std::fabs(y - reference_node -> y) < epsilon) { // the crossover matches the reference node
             reference_node -> inside = true;
             reference_node -> shared_node_number = shared_node_number;
             reference.insertNode(reference_node,compare_point+1,1);      
             compare_node_count++;
           } else if(std::fabs(x - reference_node2 -> x) < epsilon && std::fabs(y - reference_node2 -> y) < epsilon) { // the crossover matches the reference node2
             reference_node2 -> inside = true;
             reference_node2 -> shared_node_number = shared_node_number;
             reference.insertNode(reference_node2,compare_point+1,1);      
             compare_node_count++;
           } else { // shouldn't need to check the compare nodes, because that was done during the creation of the compare polygon
           // add the node to the polygons and restart the search

           
            reference.create_node(crossover_node);
            crossover_node -> x = x;
            crossover_node -> y = y;
            crossover_node -> inside = true;
            crossover_node -> shared_node_number = shared_node_number;

            reference.insertNode(crossover_node,reference_point+1,0);
            reference_node_count++;

            reference.insertNode(crossover_node,compare_point+1,1);
            compare_node_count++;


           }
           compare_point = 1;
           reference_point = 0;
           not_finished2 = false;

         } else {
   //        m_log->message(2, "* no, they don't ...\n");
          compare_point++;

         } // end is_crossover





       } // end if not_finished2

      //}

     }
      reference_point++;

    }

  // }

  }

//  reference.findNode(reference_node,1,0);
//  reference.findNode(compare_node, 1,1);


//  m_log->message(2, "* object after insert nodes: ... \n");
//  reference.print_polygon(0);
//  reference.print_polygon(1);


  /// Doing something like the Weiler Atherton clipping algorithm https://en.wikipedia.org/wiki/Weiler%E2%80%93Atherton_clipping_algorithm 
  // will require a linked list

 // m_log->message(2, "* find the overlap ...\n");
//  polygon_linked_list overlapping_polygon(3);

  node * final_node;
  node * append_node;
  node * overlap_node;
//  m_log->message(2, "* find water, %2i %2i ...\n", all_inside, all_inside2);

  all_inside = false; // fudge
  all_inside2 = false;  // fudge
  if (all_inside && ! all_inside2) {

    // the reference cell is entirely contained within the compare cell
//  m_log->message(2, "* the reference cell is entirely contained within the compare cell ...\n");
    water = find_quad_area(reference_cell) / find_quad_area(compare_cell);

  } else if (all_inside2 && ! all_inside){

    // the compare cell is entirely contained within the reference cell
//  m_log->message(2, "* the compare cell is entirely contained within the reference cell ...\n");
    water = find_quad_area(compare_cell);

  } else {

    // part of the compare cell is in the reference cell

    // start by grabbing the first node of the reference cell
 //    m_log->message(2, "* attempting to find the first node ...\n");

    int counter = 1, overlap_points = 1;
    bool found_first = false;
    bool found_node_reference, found_node_compare;
    bool use_reference;
    bool end_found;
    while(! found_first) {

 //    m_log->message(2, "* trying to find first node ...\n");
     found_node_reference = reference.findNode(reference_node, counter,0);
     found_node_compare = reference.findNode(compare_node, counter,1);

     if(found_node_reference && found_node_compare) { // have to decide which one comes first
      if(reference_node -> inside && ! compare_node -> inside) { // reference node comes first
       reference.insertNode(reference_node,overlap_points,2);
       final_node = reference_node;
       overlap_node = reference_node;
       found_first = true;
       use_reference = true;
       end_found = false;

      } else if (! reference_node -> inside &&  compare_node -> inside) { // compare node comes first

       reference.insertNode(compare_node,overlap_points,2);
       final_node = compare_node;
       overlap_node = compare_node;
       found_first = true;
       use_reference = false;
       end_found = false;
      } else if(reference_node -> inside &&  compare_node -> inside) { // both nodes are possible candidates, I'm assuming it doesn't really matter which one you choose, so I just pick the reference one. Keeping this if it this turns out to be a bad assumption.

       reference.insertNode(reference_node,overlap_points,2);
       final_node = reference_node;
       overlap_node = reference_node;
       found_first = true;
       use_reference = true;
       end_found = false;
      } else {
         counter++;
      }

     } else { // assume that there is no overlap
       found_first = true;
       water = 0.0; // no overlap
       end_found = true;
     }

    }
/*
     m_log->message(2, "* found the first node ...\n");
   // this should work assuming there isn't multiple overlapping polygons
     m_log->message(2, "* reference ...\n");
  reference.print_polygon(0);
     m_log->message(2, "* compare ...\n");
  reference.print_polygon(1);
*/
   bool found_counter;

  // overlap_node = new node; 
   while(! end_found) {

 //    m_log->message(2, "* searching ...\n");
     reference_node = overlap_node -> next[0];
     compare_node = overlap_node -> next[1];

 //    m_log->message(2, "* nodes: %p %p %p %p\n", reference_node, compare_node, overlap_node, final_node);
     

     if(reference_node != NULL && reference_node ->next[0] !=NULL && reference_node -> inside) {
//       m_log->message(2, "* reference: %15.10f %15.10f %p\n", reference_node ->x, reference_node ->y, reference_node);
 //      m_log->message(2, "* final: %15.10f %15.10f %p\n", final_node ->x, final_node ->y, final_node);
 //       m_log->message(2, "* reference node inside: %5i\n",reference_node -> inside);
       if (reference_node != final_node) {
         overlap_points++;
         overlap_node = reference_node;
         reference.insertNode(overlap_node,overlap_points,2);
       } else {

       end_found = true;
       }

     } else if(compare_node != NULL && compare_node ->next[1] !=NULL && compare_node -> inside) {
//       m_log->message(2, "* compare: %15.10f %15.10f %p\n", compare_node ->x, compare_node ->y, reference_node);
 //      m_log->message(2, "* final: %15.10f %15.10f %p\n", final_node ->x, final_node ->y, final_node);
//        m_log->message(2, "* compare node inside: %5i\n",compare_node -> inside);
      if (compare_node != final_node ) {

 //       m_log->message(2, "* compare node inside: %5i\n",compare_node -> inside);
       overlap_points++;
       overlap_node = compare_node;
       reference.insertNode(overlap_node,overlap_points,2);
      } else {

       end_found = true;
      }

     } else {

      end_found = true;
     }


   } // end while

   // TODO add function to calculate area of resulting polygon


   // for now, just set water to equal 1
// m_log->message(2, "* overlapping points: ... %5i\n", overlap_points);
   if( overlap_points < 3) {
    water = 0.0;
   }else {
    water = reference.polygon_area(2);
   }

  }

 if (printing){
 // reference.print_polygon(2);
 // m_log->message(2, "* amount of water: ... %15.10f\n", water);
  reference.print_polygon(0);
  reference.print_polygon(1); 
  reference.print_polygon(2);

 }
  // deallocate memory
 // m_log->message(2, "* deallocate overlapping_polygon ...\n");



  reference.~polygon_linked_list();

  return water;


}

bool HydrologyMod::find_crossover(node *reference1, node *reference2, node *compare1, node *compare2, double& x, double& y) {


  // if the nodes are equal, then by definition, there is no crossover
  double epsilon = 0.0000001; // If the difference between the x values are sufficiently small, I consider the line to be essentially vertical, I hope 10^-6 is good enough

  if((std::fabs(reference1 -> x - compare1 -> x) < epsilon && std::fabs(reference1 ->y - compare1 -> y) < epsilon) || 
     (std::fabs(reference1 -> x - compare2 -> x) < epsilon && std::fabs(reference1 ->y - compare2 -> y) < epsilon) || 
     (std::fabs(reference2 -> x - compare1 -> x) < epsilon && std::fabs(reference2 ->y - compare1 -> y) < epsilon) || 
     (std::fabs(reference2 -> x - compare2 -> x) < epsilon && std::fabs(reference2 ->y - compare2 -> y) < epsilon)) {

//   m_log->message(2, "*lines have a shared node: ...\n");
   return false;
  }



  // good reference for the problems with floating point math: https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/

  double slope_reference, slope_compare, intercept_reference, intercept_compare;

  // first check if the slope of the reference points are infinite
  if (fabs(reference1 -> x - reference2 -> x) < epsilon) {
   // also check the reference
   if (fabs(compare1 -> x - compare2 -> x) < epsilon) { // lines are likely parallel
//     m_log->message(2, "*lines are considered to be parallel 1: ...\n");
     return false; // right now no check of the lines completely overlap, should do that in another part of the program
   } else {

    slope_compare = (compare1 -> y - compare2 -> y) / (compare1 -> x - compare2 -> x);
    intercept_compare = compare1 -> y - compare1 -> x * slope_compare;

    x = reference1 -> x;
    y = reference1 -> x * slope_compare + intercept_compare;

   }


  } else {

   slope_reference = (reference1 -> y - reference2 -> y) / (reference1 -> x - reference2 -> x);
   intercept_reference = reference1 -> y - reference1 -> x * slope_reference;

   if (fabs(compare1 -> x - compare2 -> x) < epsilon) { // assume that the line is vertical
 //    m_log->message(2, "*reference line is vertical: ...\n");
    x = compare1 -> x;
    y = compare1 -> x * slope_reference + intercept_reference;
     
   } else { 

    slope_compare = (compare1 -> y - compare2 -> y) / (compare1 -> x - compare2 -> x);
    intercept_compare = compare1 -> y - compare1 -> x * slope_compare;

    if(fabs(slope_reference - slope_compare) < epsilon) { // lines are likely parallel
 //    m_log->message(2, "*lines are considered to be parallel 2: ...\n");
      return false; // right now no check of the lines completely overlap, should do that in another part of the program
    }

    x = (intercept_reference - intercept_compare) / (slope_compare - slope_reference);
    y = x * slope_compare + intercept_compare;

   }

  }

/*
 &&
     std::fabs( x - compare1 -> x) > epsilon && std::fabs( x - compare2 -> x) > epsilon &&
     std::fabs( x - reference1 -> x) > epsilon && std::fabs( x - reference2 -> x) > epsilon &&
     std::fabs( y - compare1 -> y) > epsilon && std::fabs( y - compare2 -> y) > epsilon &&
     std::fabs( y - reference1 -> y) > epsilon && std::fabs( y - reference2 -> y) > epsilon
*/

 //    m_log->message(2, "*x and y: %15.10f %15.10f\n", x, y);
  // check if the crossover is between the two lines, and not lying directly on one of the other nodes
  // had to add an epsilon factor because it was failing if x1 == x2 or y1 == y2
  if(x <= std::max(compare1 -> x, compare2 -> x)+epsilon &&  x >= std::min(compare1 -> x, compare2 -> x)-epsilon && 
      y <= std::max(compare1 -> y, compare2 -> y)+epsilon &&  y >= std::min(compare1 -> y, compare2 -> y)-epsilon &&
      x <= std::max(reference1 -> x, reference2 -> x)+epsilon &&  x >= std::min(reference1 -> x, reference2 -> x)-epsilon && 
      y <= std::max(reference1 -> y, reference2 -> y)+epsilon &&  y >= std::min(reference1 -> y, reference2 -> y)-epsilon) {
 //    m_log->message(2, "* lines overlap: ...\n");
     return true;
  } else {
//     m_log->message(2, "* lines do not overlap: ...\n");
    return false;

  }

  // if it got this far without returning, there is obviously something wrong with the function

  m_log->message(2,
             "* warning: error in find_crossover()\n");

  return false;

}
/*
 bool HydrologyMod::find_shared_node(struct node *& node_check, polygon_linked_list& list_to_check) {

  m_log->message(2,
             "* searching for shared node\n");
  double epsilon = 0.000001; // If the difference between the x and y values is sufficiently small, I consider the points to be equivalent

  node * reference_node;

  bool continue_search = true;

  int counter = 0;

   m_log->message(2,"just checking things out: %p %p %p\n", list_to_check, node_check, reference_node);
  continue_search = list_to_check.findNode(reference_node,1);
 // list_to_check.findNode(reference_node,2);

  while (continue_search) {

   counter++;

   m_log->message(2,"Checking point: %5i\n", counter);
   continue_search = list_to_check.findNode(reference_node,counter);
   m_log->message(2,"x: %15.10f %15.10f \n", node_check -> x, reference_node -> x);
   m_log->message(2,"y: %15.10f %15.10f \n", node_check -> y, reference_node -> y);
   m_log->message(2,"inside: %2i %2i \n", node_check -> inside, reference_node -> inside);
   m_log->message(2,"next: %p %p \n", node_check -> next[0], reference_node -> next[0]);

   if(std::fabs(node_check -> x - reference_node -> x) < epsilon && std::fabs(node_check -> y - reference_node -> y) < epsilon) {

    node_check = reference_node;
    return true;

   }
  }

  return false;
 }

*/

bool HydrologyMod::point_in_polygon(double polygon[][2], int polygon_size, double x, double y, bool on_edge ) {

  // second index of polygon should be size 2 (for x and y), the other index is polygon_size

  bool inside = false;

  int current_point, next_point;



  for(current_point=0; current_point<polygon_size; current_point++) {

    if(current_point == polygon_size - 1) {
      next_point = 0;
    } else {
      next_point = current_point + 1;
    }

     // even-odd rule algorithm to determine if the point is inside or outside

    if (std::min(polygon[current_point][1], polygon[next_point][1]) < y && std::max(polygon[current_point][1], polygon[next_point][1]) >= y) {

     if (polygon[current_point][0] + (y - polygon[current_point][1]) / (polygon[next_point][1] - polygon[current_point][1]) * (polygon[next_point][0] - polygon[current_point][0]) < x) {

      inside = ! inside;

     }

    }

  }

  if(inside) {
//   m_log->message(2,"node is inside: %15.10f %15.10f\n", x, y);
  } else {
  //    m_log->message(2,"node is outside, checking to see if it is on the edge\n");
  }

  if(!inside && on_edge) { // find out if the point is on the edge of the polygon

   double epsilon = 0.0000001; // If the difference between the x values are sufficiently small, I consider the line to be essentially vertical, I hope 10^-6 is good enough
   double slope, intercept;

   for(current_point=0; current_point<polygon_size; current_point++) {

    if(current_point == polygon_size - 1) {
      next_point = 0;
    } else {
      next_point = current_point + 1;
    }

    if(std::fabs(polygon[current_point][0] - polygon[next_point][0]) < epsilon) { // probably a vertical line
      if (std::fabs(x - std::fabs(polygon[current_point][0])) < epsilon) { // point has the same x, but is it in the range of y?
        if (y >= std::min(polygon[current_point][1],polygon[next_point][1]) && y <= std::min(polygon[current_point][1],polygon[next_point][1])) { // falls on the line

    //    m_log->message(2,"node is presumably on the edge (vertical)\n");
   //       m_log->message(2,"%15.10f %15.10f\n", x, y );
        for(int current_point2=0; current_point2<polygon_size; current_point2++) {
     //     m_log->message(2,"%15.10f %15.10f\n", polygon[current_point2][0], polygon[current_point2][1] );
        }
    //   m_log->message(2,"node is on edge 1: %15.10f %15.10f\n", x, y);
         return true;
        }
      }
    } else {

      slope = (polygon[current_point][1] - polygon[next_point][1]) / (polygon[current_point][0] - polygon[next_point][0]);
      intercept = polygon[current_point][1] - polygon[current_point][0] * slope;
      double temp_y = slope*x+intercept;

      if(std::fabs(y - temp_y) < epsilon && x >= std::min(polygon[current_point][0],polygon[next_point][0]) && x <= std::max(polygon[current_point][0],polygon[next_point][0])) {


  //           m_log->message(2,"node is on edge 2: %15.10f %15.10f\n", x, y);

  //        m_log->message(2,"points: %15.10f %15.10f %15.10f %15.10f %15.10f\n", x, y, temp_y, slope, intercept );
  //        m_log->message(2,"points: %15.10f %15.10f %15.10f %15.10f \n", polygon[current_point][0], polygon[current_point][1], polygon[next_point][0], polygon[next_point][1]);

       return true;
      } // end if

    } // end if


   } // end for

  } // end if

  if(inside){
     //   m_log->message(2,"node is presumably inside\n");
     //     m_log->message(2,"%15.10f %15.10f\n", x, y);
        for(int current_point2=0; current_point2<polygon_size; current_point2++) {
    //      m_log->message(2,"%15.10f %15.10f\n", polygon[current_point2][0], polygon[current_point2][1] );
        }
  }
  return inside;
}

// put variables you want to output to file in these functions
void HydrologyMod::define_model_state_impl(const PIO &output) const {
  Hydrology::define_model_state_impl(output);
  m_tillwat_flux.define(output);
  m_excess_water.define(output);
  m_fraction_till.define(output);
  m_total_input.define(output);
  m_gradient_magnitude.define(output);
  m_theta.define(output);
  m_Pover_ghosts.define(output);
  m_pressure_gradient.define(output);
}

void HydrologyMod::write_model_state_impl(const PIO &output) const {
  Hydrology::write_model_state_impl(output);
  m_tillwat_flux.write(output);
  m_excess_water.write(output);
  m_fraction_till.write(output);
  m_total_input.write(output);
  m_gradient_magnitude.write(output);
  m_theta.write(output);
  m_Pover_ghosts.write(output);
  m_pressure_gradient.write(output);
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

//  m_log->message(2,"**************************************************** ...\n");
//  m_log->message(2,"dx: %15.10f \n", dx);
//  m_log->message(2,"dx: %15.10f \n", dy);

//  m_log->message(2,"**************************************************** ...\n");


  IceModelVec::AccessList list{&m_Pover_ghosts};

  Hydrology::overburden_pressure(m_Pover_ghosts); // grab the overburden pressure
  m_Pover_ghosts.update_ghosts();
  const IceModelVec2CellType &mask = *m_grid->variables().get_2d_cell_type("mask");
  list.add(result);
  list.add(result_mag);
  list.add(result_angle);
  list.add(mask);


  // third order finite difference method for calculating gradient, should give good results (Skidmore, 1989)
  for (Points p(*m_grid); p; p.next()) {										
    const int i = p.i(), j = p.j();


    // I think I need to change the gradient calculation
       

    if (mask.icy(i, j) && (m_Pover_ghosts(i+1,j+1) > m_Pover_ghosts(i,j) ||
       m_Pover_ghosts(i,j+1) > m_Pover_ghosts(i,j) ||
       m_Pover_ghosts(i-1,j+1) > m_Pover_ghosts(i,j) ||
       m_Pover_ghosts(i+1,j) > m_Pover_ghosts(i,j) ||
       m_Pover_ghosts(i-1,j) > m_Pover_ghosts(i,j) ||
       m_Pover_ghosts(i+1,j-1) > m_Pover_ghosts(i,j) ||
       m_Pover_ghosts(i,j-1) > m_Pover_ghosts(i,j) ||
       m_Pover_ghosts(i-1,j-1) > m_Pover_ghosts(i,j))){
      result(i,j).u = ((m_Pover_ghosts(i+1,j+1) + 2.0 * m_Pover_ghosts(i+1,j) + m_Pover_ghosts(i+1,j-1)) -
			    (m_Pover_ghosts(i-1,j+1) + 2.0 * m_Pover_ghosts(i-1,j) + m_Pover_ghosts(i-1,j-1))) / (8.0 * dx);
  //
/*
      if((i== 30 || i==29) && j==28){
      m_log->message(2,"> %5i %5i %5i\n", i-1, i, i+1 );
      m_log->message(2,"%5i %12.2f %12.2f %12.2f \n",j+1, m_Pover_ghosts(i-1,j+1), m_Pover_ghosts(i,j+1), m_Pover_ghosts(i+1,j+1));
      m_log->message(2,"%5i %12.2f %12.2f %12.2f \n",j, m_Pover_ghosts(i-1,j), m_Pover_ghosts(i,j), m_Pover_ghosts(i-1,j));
      m_log->message(2,"%5i %12.2f %12.2f %12.2f \n",j-1, m_Pover_ghosts(i-1,j-1), m_Pover_ghosts(i,j-1), m_Pover_ghosts(i+1,j-1));
      }
*/
      result(i,j).v = ((m_Pover_ghosts(i+1,j+1) + 2.0 * m_Pover_ghosts(i,j+1) + m_Pover_ghosts(i-1,j+1)) -
			    (m_Pover_ghosts(i+1,j-1) + 2.0 * m_Pover_ghosts(i,j-1) + m_Pover_ghosts(i-1,j-1))) / (8.0 * dy);

//      m_log->message(2,"gradient u and v: %5i %5i %12.2f %12.2f \n", i, j, result(i,j).u, result(i,j).v );

     

    } else {
      result(i,j).u = 0.0;

      result(i,j).v = 0.0;
    }

  }

  m_log->message(2,"* calculated gradient, now removing holes ...\n");

  // check to see if there are adjacent grid cells with gradient vector components that oppose each other. If so, the one
  // that has a larger magnitude wins, and the smaller one gets set to zero. This is why the ghosts are reset.
  // Obviously, not exactly mass conserving, basically any water that goes into one of these grid cells will become a black
  // hole for water using my "drain completely" scheme.
  // After, the magnitude of the gradient is calculated. Also, the direction of water flow is calculated.

  result.update_ghosts();
  double mag_before;
  for (Points p(*m_grid); p; p.next()) {										
    const int i = p.i(), j = p.j();

  //  m_log->message(2,"gradient u and v before black hole: %12.2f %12.2f \n", result(i,j).u, result(i,j).v );
    mag_before = sqrt(pow(result(i,j).u,2) + pow(result(i,j).v,2));
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

//    m_log->message(2,"gradient u and v after black hole: %12.2f %12.2f \n", result(i,j).u, result(i,j).v );

    result_mag(i,j) = sqrt(pow(result(i,j).u,2) + pow(result(i,j).v,2));
    if(mask.icy(i,j)){
  //  m_log->message(2,"%5i %5i %15.10f %15.10f %15.10f \n", i, j,  mag_before, result_mag(i,j),  mag_before-result_mag(i,j));
    }

    // test fudge
    // 3 is stable
    // 2 is stable
    // 1 is stable, takes a bit longer to stabilize
    // 0.5 takes longer to become stable, but eventually gets there
    // 0.25 is not stable, it reaches a non-symettric state in EISMINT
/*
    if(result_mag(i,j) < 2) {
      result_mag(i,j) = 0.0;
    }
*/

  //  m_log->message(2,"mag_grad: %5i % 5i %12.2f %12.2f %12.2f %12.2f %12.2f \n", i, j, result_mag(i,j), m_Pover_ghosts(i-1,j), m_Pover_ghosts(i,j+1),  m_Pover_ghosts(i+1,j),  m_Pover_ghosts(i,j-1) );
//    m_log->message(2,"gradient u and v after black hole: %12.2f %12.2f %12.2f \n", result(i,j).u, result(i,j).v, result_mag(i,j) );
    if(result_mag(i,j) > 0) {
      result_angle(i,j) = atan2 (-result(i,j).v,-result(i,j).u);
   //   m_log->message(2,"gradient u and v and direction: %5i %5i %15.10f %15.10f %15.10f \n", i, j, result(i,j).u, result(i,j).v, result_angle(i,j) );
   


    }
    else {
      result_angle(i,j) = 4.0;
    }

   if(mask.icy(i,j)) {
  //  m_log->message(2,"angle: %5i %5i %15.10f  \n", i, j, result_angle(i,j) );
 //  m_log->message(2,"gradient: %5i %5i %15.10f %15.10f %15.10f %15.10f %15.10f \n", i, j, result(i,j).u, result(i,j).v, result_mag(i,j), result_angle(i,j),  m_Pover_ghosts(i,j) );
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
  const IceModelVec2CellType &mask = *m_grid->variables().get_2d_cell_type("mask");
  list.add(mask);
  
  const double water_viscosity = m_config->get_double("hydrology.water_viscosity");
  double drainage;
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // Darcy's Law, combined with fraction sediments
 //   result(i,j) =  m_till_permeability(i,j) / water_viscosity * m_gradient_magnitude(i,j)* m_fraction_till(i,j);  // flux out of the cell from sediments
     result(i,j) =  m_till_permeability(i,j) / water_viscosity * m_gradient_magnitude(i,j);  // flux out of the cell from sediments
   // of course, it is possible that the rate of discharge cannot be so much that it completely drains the till layer before the end of the next time step
   // I'm setting it so that the rate can only be at maximum complete drainage during the time step

   double drainage_thickness = result(i,j) * m_dt;
//   m_log->message(2,"%5i %5i %14.10f %14.10f %14.10f \n", i, j, result(i,j), drainage_thickness,  m_gradient_magnitude(i,j));
   if(mask.icy(i,j)) {
 //    m_log->message(2,"%5i %5i %14.10f %14.10f \n", i, j, result(i,j),  m_gradient_magnitude(i,j));

   }


   if(m_Wtil(i,j)*m_fraction_till(i,j) < drainage_thickness) {
//     result(i,j) = m_Wtil(i,j)*m_fraction_till(i,j) / m_dt;  // commenting out for now
   }
//   drainage = result(i,j)*m_dt;
//   m_log->message(2,"* till_drainage: %14.10f \n", drainage);
  }  // end for

  m_log->message(2,"* finished calculating basic input, now find how much drains out ...\n");

  result.update_ghosts();


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


// adding in the linked list class
//class polygon_linked_list;

 polygon_linked_list::polygon_linked_list(){ // constructor

//  std::cout << "creating linked list " << "\n";
//   log.message(2,"* creating linked list: %1i ...\n", number);

   listLength[0] = 0;
   listLength[1] = 0;
   listLength[2] = 0;

//  std::cout << "before head " << listLength << "\n";
   create_node(head);
/*
   head = new node;
   head -> x = 0.0;
   head -> y = 0.0;

   head -> next[0] = NULL;
   head -> next[1] = NULL;
   head -> next[2] = NULL;
   head -> shared_node_number = 0;
*/

//  std::cout << "finished initializing polygon_linked_list "  << "\n";
 }

 polygon_linked_list::~polygon_linked_list(){ // destructor

/*
  for(int polygon_number = 0; polygon_number<=2; polygon_number++){
   node * p = head;
   node * q = head;

   int check1, check2, counter;
   std::cout << "attempting to destruct " << polygon_number << "\n";
   counter = 0;
   bool continue_loop = true;
   while (continue_loop){
    counter++;
    if(counter <= listLength[polygon_number]) {

    p = q;
    q = p -> next[polygon_number];

    if(p -> next[polygon_number]) {
      std::cout << "not pointing to NULL " << counter << " " << polygon_number << "\n";



     if(polygon_number == 0) {
      check1 = 1;
      check2 = 2;
     } else if (polygon_number == 1) {
      check1 = 0;
      check2 = 2;
     } else {
      check1 = 0;
      check2 = 1;
     }

     if ( p -> next[check1] == NULL &&  p -> next[check2] == NULL ) { // safe to delete
      std::cout << "destroying node " << counter << " " << listLength[polygon_number] << "\n";
       delete p;
     } else { // get rid of the pointer
      std::cout << "leaving node where it is " << counter << " " << listLength[polygon_number] << "\n";
      p -> next[polygon_number] = NULL; 

     }

    }else{
      std::cout << "pointing to NULL " << counter << " " << polygon_number << "\n";
      continue_loop = false;
    }

    } else{
     continue_loop = false;
    }

   }

  std::cout << "finished destructing " << polygon_number << "\n";

  }
*/
 }


 void polygon_linked_list::insertNode(struct node * newNode, int position, int polygon_number ){
  if(polygon_number == 0) {
//    std::cout << "inside insertNode reference " << polygon_number << " " << position << " "  << listLength[polygon_number] << "\n";
  } else if (polygon_number == 1) {
 //   std::cout << "inside insertNode compare " << polygon_number << " " << position << " "  << listLength[polygon_number] << "\n";
  } else {
 //   std::cout << "inside insertNode overlap " << polygon_number << " " << position << " "  << listLength[polygon_number] << "\n";
  }
//  std::cout << "check head node: " << head << " " << head -> next[0] << " " << head -> next[1] << "\n";
  if ((position <= 0) || (position > listLength[polygon_number] + 1)){
   //     m_log_local->message(2,"* There is something wrong with code to insert a node (1)!!!!\n");     
  } // end if

  if (head -> next[polygon_number] == NULL){
//    std::cout << "attempting to add new pointer at head\n";
    head -> next[polygon_number] = newNode;
    listLength[polygon_number]++;
//    std::cout << "added new pointer at head\n";
    return;
  } // end if
//    std::cout << "attempting to add new pointer after head, checking first:\n";
  int count = 0;
  node * p = head;
  node * q = head -> next[polygon_number];
  for (count=1; count<=listLength[polygon_number]; count++){
 //   std::cout << "pointer_check insert: " << q << " " << count << " " << q -> x << " " << q ->y << " " << q -> next[polygon_number] << "\n";
    q = q -> next[polygon_number];
  }

  q = head -> next[polygon_number];

 // std::cout << "adding node after head " << polygon_number << " " << position << " "  << listLength[polygon_number] << "\n";
//  std::cout << "x, y, inside: " << newNode -> x << " "  << newNode -> y << " "  << newNode -> inside << "\n";

  for (count=1; count<=listLength[polygon_number]; count++){
 
 //   std::cout << "count " << count << " " << position << " " <<  q <<  "\n";
    if (count == position){
 //    std::cout << "attempting to add the pointer \n";
      p -> next[polygon_number] = newNode;
      newNode -> next[polygon_number] = q;
 //    std::cout << "should be added \n";
      listLength[polygon_number]++;
      return;
    } // end if
    p = q;
    q = p -> next[polygon_number];

  } // end for

//  std::cout << "adding node at the end of the list " << polygon_number << " " << position << " "  << listLength[polygon_number] << "\n";
  if (position == listLength[polygon_number]+1){
    p -> next[polygon_number] = newNode;
 //   newNode -> next[polygon_number] = q;
    listLength[polygon_number]++;
    return;
  } // end if
  // couldn't add point for some reason
//  m_log_local->message(2,"* There is something wrong with code to insert a node (2)!!!!\n");     
 
 } // end polygon_linked_list::insertNode

 bool polygon_linked_list::findNode(struct node *& node_out, int position, int polygon_number) {

  node_out = NULL;
//  std::cout << "inside findNode " << polygon_number << " "  << position << " " << listLength[polygon_number] << "\n";
 // std::cout << "check head node: " << head << " " << head -> next[0] << " " << head -> next[1] << "\n";
//   std::cout << "also check passed node: " << node_out <<  "\n"; 
  if ((position <= 0) || (position > listLength[polygon_number] )){
//     std:: cout << "Returning false\n";
      return false;  
  } // end if

  int count = 0;
// std::cout << "checking the head node " << head << " " << head -> next[polygon_number] << "\n";

  node * q = head -> next[polygon_number];



//  std::cout << "assigned q " << polygon_number << " " << q << "\n";

  for(count = 1; count <= listLength[polygon_number]; count++){

//  std::cout << "pointer_check: " << q << " " << count << " " << q -> x << " " << q ->y << " " << q -> next[polygon_number] << "\n";
  q = q -> next[polygon_number];

  }

  q = head -> next[polygon_number];

  for(count = 1; count <= listLength[polygon_number]; count++){
//   std::cout << "searching " << count << " " << position << " " << q -> x << "\n";
   if(count == position) {
 //    std::cout << "found, or at least it should be: \n";
 //    std::cout << "found: " << count << " " << q -> x << " " << q ->y << " " << q -> next[polygon_number] << "\n";
      node_out = q;

 //    std::cout << "found: " << count << " " << node_out -> x << " " << node_out ->y << " " << node_out -> next[polygon_number] << "\n";
      return true;
    } // end if
    q = q -> next[polygon_number];

  } // end while

  //   std::cout << "didn't find: " << count << position << "\n";
  return false;

 } // end polygon_linked_list::findNode

  void polygon_linked_list::create_node(struct node *& node_out){

   node_out =  new node;

   node_out -> x = 0.0;
   node_out -> y = 0.0;
   node_out -> inside = false;
   node_out -> next[0] = NULL;
   node_out -> next[1] = NULL;
   node_out -> next[2] = NULL;
   node_out -> shared_node_number = 0;

 }

 void polygon_linked_list::print_polygon(int polygon_number) {

  if(polygon_number == 0) {
   // std::cout << "reference: " << polygon_number << " "  << listLength[polygon_number] << "\n";
   std::cout << "> -Z-1" << "\n";
  } else if (polygon_number == 1) {
   // std::cout << "compare: " << polygon_number << " "  << listLength[polygon_number] << "\n";
   std::cout << "> -Z0" << "\n";
  } else {
   std::cout << "> -Z1" << "\n";
    if(listLength[polygon_number] == 0) {
     return;
    }

//    std::cout << "overlap: " << polygon_number << " "  << listLength[polygon_number]+1 << "\n";
  }

  if(listLength[polygon_number] == 0) {
 //  std:: cout << "no points \n";
  } else {

   node * q = head -> next[polygon_number];

   for(int count = 1; count <= listLength[polygon_number]; count++){

    if(q -> inside) {
     std::cout << q -> x << " " << q ->y << " true " << "\n";
//     std::cout << q -> x << " " << q ->y << "\n";
    } else {
     std::cout << q -> x << " " << q ->y << " false "  << "\n";
//     std::cout << q -> x << " " << q ->y << "\n";  
    }
    q = q -> next[polygon_number];

   }
  }
 
  if(polygon_number==2) {
   node * q = head -> next[polygon_number];
   std::cout << q -> x << " " << q ->y << "\n";
  }
 }

 double polygon_linked_list::polygon_area(int polygon_number) {

  // http://mathworld.wolfram.com/PolygonArea.html

  double area = 0.0;

  double final_x, final_y;
  if(listLength[polygon_number] > 1) { // by definition, the last node is not included

   node * q = head -> next[polygon_number];
   final_x = q->x;
   final_y = q->y;
   node * p = q -> next[polygon_number];
   
   for(int count = 1; count < listLength[polygon_number]; count++){

    area = area + q->x * p->y - q->y * p->x;
  //  std::cout << "area sum: " << count << " " << area << " " <<  "\n";
    q = p;
    p = q -> next[polygon_number];   

   }

   area = area + q->x * final_y - q->y * final_x;
  //  std::cout << "area sum: " << listLength[polygon_number] << " " << area << " " <<  "\n";
   area = 0.5 * std::fabs(area);

   // std::cout << "final sum: " << listLength[polygon_number] << " " << area << " " <<  "\n";

  }
  
  return area;

 }

} // end namespace hydrology
} // end namespace pism+;
