// Copyright (C) 2011, 2012, 2013, 2014, 2015 PISM Authors
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

#ifndef _PISMMOHRCOULOMBYIELDSTRESSMOD_H_
#define _PISMMOHRCOULOMBYIELDSTRESSMOD_H_

#include "PISMYieldStress.hh"
#include "PISMMohrCoulombYieldStress.hh"
#include "base/util/iceModelVec.hh"

namespace pism {

namespace hydrology {
class Hydrology;
class HydrologyMod;
}

class MohrCoulombYieldStressMod : public MohrCoulombYieldStress {
public:

  MohrCoulombYieldStressMod(IceGrid::ConstPtr g, hydrology::Hydrology *hydro);
  virtual ~MohrCoulombYieldStressMod();



protected:

  virtual void update_impl();
  virtual void init_impl();

  virtual void define_model_state_impl(const PIO &output) const;
  virtual void write_model_state_impl(const PIO &output) const;

protected:

  IceModelVec2S m_fraction_till_internal, m_fraction_channel_internal;

  hydrology::HydrologyMod *m_hydrology;

};

} // end of namespace pism

#endif /* _PISMMOHRCOULOMBYIELDSTRESSMOD_H_ */
