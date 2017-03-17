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

//#ifndef _PISMHYDROLOGY_H_
//#define _PISMHYDROLOGY_H_

#ifndef _PISMHYDROLOGYMOD_H_
#define _PISMHYDROLOGYMOD_H_


#include "PISMHydrology.hh"
#include "base/util/iceModelVec.hh"
#include "base/util/PISMComponent.hh"

namespace pism {

class IceModelVec2T;

namespace stressbalance {
class StressBalance;
}

//! @brief Sub-glacial hydrology models and related diagnostics.
namespace hydrology {


class HydrologyMod : public Hydrology {
public:

  HydrologyMod(IceGrid::ConstPtr g);
  virtual ~HydrologyMod();

  virtual void init();

  virtual void subglacial_water_thickness(IceModelVec2S &result) const;

  virtual void subglacial_water_pressure(IceModelVec2S &result) const;



protected:
  IceModelVec2S m_fraction_till, m_fraction_channel, m_Pover_ghosts, m_gradient_magnitude, m_tillwat_flux, m_excess_water, m_till_permeability;

  IceModelVec2V m_pressure_gradient;
  virtual void get_input_rate(double hydro_t, double hydro_dt, IceModelVec2S &result);
  virtual void update_impl(double t, double dt);
//  virtual std::map<std::string, Diagnostic::Ptr> diagnostics_impl() const;

  virtual void define_model_state_impl(const PIO &output) const;
  virtual void write_model_state_impl(const PIO &output) const;

  virtual void pressure_gradient(IceModelVec2V &result, IceModelVec2S &result_mag, IceModelVec2S &result_angle);
  virtual void till_drainage(IceModelVec2S &result, double dt);
  virtual void tunnels(IceModelVec2S &result);

  virtual MaxTimestep max_timestep_impl(double t) const;

private:
  IceModelVec2S m_theta, m_excess_water_playground,  m_excess_water_removed, quad_area; 
  IceModelVec2V bottom_left, top_left, bottom_right, top_right;

};

// linked list for creating a polygon
// I used this nice tutorial as a basis for this code: http://pumpkinprogrammer.com/2014/06/13/c-tutorial-intro-to-linked-lists/

struct node
{
    double x;
    double y;
    bool inside;
    int shared_node_number;
    node * next;
};


class polygon_linked_list{

public:

   polygon_linked_list(); // constructor
   ~polygon_linked_list(); // destructor
   void insertNode( node * newNode, int position );
   removeNode( int position );

private:
    // Head of the list contains no song data, 
    // but points to the song playlist.
    node * head;
    int listLength;

}

} // end namespace hydrology


} // end namespace pism


#endif /* _PISMHYDROLOGYMOD_H_ */