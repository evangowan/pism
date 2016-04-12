// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2016 Ed Bueler and Constantine Khroulev
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

#include "PISMBedSmoother.hh"
#include "base/util/Mask.hh"
#include "base/util/IceGrid.hh"
#include "base/util/petscwrappers/Vec.hh"
#include "base/util/IceModelVec2CellType.hh"

#include "base/util/error_handling.hh"
#include "base/util/pism_const.hh"
#include "base/util/Logger.hh"

namespace pism {
namespace stressbalance {

BedSmoother::BedSmoother(IceGrid::ConstPtr g, int MAX_GHOSTS)
    : grid(g), config(g->ctx()->config()) {

  const Logger &log = *grid->ctx()->log();

  {
    // allocate Vecs that live on all procs; all have to be as "wide" as any of
    //   their prospective uses
    topgsmooth.create(grid, "topgsmooth", WITH_GHOSTS, MAX_GHOSTS);
    topgsmooth.set_attrs("bed_smoother_tool",
                         "smoothed bed elevation, in bed roughness parameterization",
                         "m", "");
    maxtl.create(grid, "maxtl", WITH_GHOSTS, MAX_GHOSTS);
    maxtl.set_attrs("bed_smoother_tool",
                    "maximum elevation in local topography patch, in bed roughness parameterization",
                    "m", "");
    C2.create(grid, "C2bedsmooth", WITH_GHOSTS, MAX_GHOSTS);
    C2.set_attrs("bed_smoother_tool",
                 "polynomial coeff of H^-2, in bed roughness parameterization",
                 "m2", "");
    C3.create(grid, "C3bedsmooth", WITH_GHOSTS, MAX_GHOSTS);
    C3.set_attrs("bed_smoother_tool",
                 "polynomial coeff of H^-3, in bed roughness parameterization",
                 "m3", "");
    C4.create(grid, "C4bedsmooth", WITH_GHOSTS, MAX_GHOSTS);
    C4.set_attrs("bed_smoother_tool",
                 "polynomial coeff of H^-4, in bed roughness parameterization",
                 "m4", "");

    // allocate Vecs that live on processor 0:
    topgp0 = topgsmooth.allocate_proc0_copy();
    topgsmoothp0 = topgsmooth.allocate_proc0_copy();
    maxtlp0 = maxtl.allocate_proc0_copy();
    C2p0 = C2.allocate_proc0_copy();
    C3p0 = C3.allocate_proc0_copy();
    C4p0 = C4.allocate_proc0_copy();
  }

  m_Glen_exponent = config->get_double("sia_Glen_exponent"); // choice is SIA; see #285
  m_smoothing_range = config->get_double("bed_smoother_range");

  if (m_smoothing_range > 0.0) {
    log.message(2,
                "* Initializing bed smoother object with %.3f km half-width ...\n",
                units::convert(grid->ctx()->unit_system(), m_smoothing_range, "m", "km"));
  }

  // Make sure that Nx and Ny are initialized. In most cases SIAFD::update() will call
  // preprocess_bed() and set appropriate values, but in a zero-length (-y 0) run IceModel does not
  // call SIAFD::update()... We may need to re-structure this class so that everything is
  // initialized right after construction and users don't have to call preprocess_bed() manually.
  Nx = -1;
  Ny = -1;
}


BedSmoother::~BedSmoother() {
  // empty
}

/*!
Input lambda gives physical half-width (in m) of square over which to do the
average.  Only square smoothing domains are allowed with this call, which is the
default case.
 */
void BedSmoother::preprocess_bed(const IceModelVec2S &topg) {

  if (m_smoothing_range <= 0.0) {
    // smoothing completely inactive.  we transfer the original bed topg,
    //   including ghosts, to public member topgsmooth ...
    topg.update_ghosts(topgsmooth);
    // and we tell get_theta() to return theta=1
    Nx = -1;
    Ny = -1;
    return;
  }

  // determine Nx, Ny, which are always at least one if m_smoothing_range > 0
  Nx = static_cast<int>(ceil(m_smoothing_range / grid->dx()));
  Ny = static_cast<int>(ceil(m_smoothing_range / grid->dy()));
  if (Nx < 1) {
    Nx = 1;
  }
  if (Ny < 1) {
    Ny = 1;
  }

  preprocess_bed(topg, Nx, Ny);
}

const IceModelVec2S& BedSmoother::get_smoothed_bed() {
  return topgsmooth;
}

/*!
Inputs Nx,Ny gives half-width in number of grid points, over which to do the
average.
 */
void BedSmoother::preprocess_bed(const IceModelVec2S &topg,
                                 unsigned int Nx_in, unsigned int Ny_in) {

  if ((Nx_in >= grid->Mx()) || (Ny_in >= grid->My())) {
    throw RuntimeError("input Nx, Ny in bed smoother is too large because\n"
                       "domain of smoothing exceeds IceGrid domain");
  }
  Nx = Nx_in; Ny = Ny_in;

  topg.put_on_proc0(*topgp0);
  smooth_the_bed_on_proc0();
  // next call *does indeed* fill ghosts in topgsmooth
  topgsmooth.get_from_proc0(*topgsmoothp0);

  compute_coefficients_on_proc0();
  // following calls *do* fill the ghosts
  maxtl.get_from_proc0(*maxtlp0);
  C2.get_from_proc0(*C2p0);
  C3.get_from_proc0(*C3p0);
  C4.get_from_proc0(*C4p0);
}


/*!
Call preprocess_bed() first.
 */
void BedSmoother::get_smoothing_domain(int &Nx_out, int &Ny_out) {
  Nx_out = Nx;
  Ny_out = Ny;
}


//! Computes the smoothed bed by a simple average over a rectangle of grid points.
void BedSmoother::smooth_the_bed_on_proc0() {

  ParallelSection rank0(grid->com);
  try {
    if (grid->rank() == 0) {
      petsc::VecArray2D
        b0(*topgp0,       grid->Mx(), grid->My()),
        bs(*topgsmoothp0, grid->Mx(), grid->My());

      for (int j=0; j < (int)grid->My(); j++) {
        for (int i=0; i < (int)grid->Mx(); i++) {
          // average only over those points which are in the grid; do
          // not wrap periodically
          double sum = 0.0, count = 0.0;
          for (int r = -Nx; r <= Nx; r++) {
            for (int s = -Ny; s <= Ny; s++) {
              if ((i+r >= 0) and (i+r < (int)grid->Mx()) and
                  (j+s >= 0) and (j+s < (int)grid->My())) {
                sum += b0(i+r, j+s);
                count += 1.0;
              }
            }
          }
          // unprotected division by count but r=0,s=0 case guarantees count>=1
          bs(i, j) = sum / count;
        }
      }
    }
  } catch (...) {
    rank0.failed();
  }
  rank0.check();
}


void BedSmoother::compute_coefficients_on_proc0() {

  const unsigned int Mx = grid->Mx(), My = grid->My();

  ParallelSection rank0(grid->com);
  try {
    if (grid->rank() == 0) {
      petsc::VecArray2D
        b0(*topgp0,       Mx, My),
        bs(*topgsmoothp0, Mx, My),
        mt(*maxtlp0,      Mx, My),
        c2(*C2p0,         Mx, My),
        c3(*C3p0,         Mx, My),
        c4(*C4p0,         Mx, My);

      for (int j=0; j < (int)My; j++) {
        for (int i=0; i < (int)Mx; i++) {
          // average only over those points which are in the grid
          // do not wrap periodically
          double
            topgs     = bs(i, j),
            maxtltemp = 0.0,
            sum2      = 0.0,
            sum3      = 0.0,
            sum4      = 0.0,
            count     = 0.0;

          for (int r = -Nx; r <= Nx; r++) {
            for (int s = -Ny; s <= Ny; s++) {
              if ((i+r >= 0) && (i+r < (int)Mx) && (j+s >= 0) && (j+s < (int)My)) {
                // tl is elevation of local topography at a pt in patch
                const double tl  = b0(i+r, j+s) - topgs;
                maxtltemp = std::max(maxtltemp, tl);
                // accumulate 2nd, 3rd, and 4th powers with only 3 multiplications
                const double tl2 = tl * tl;
                sum2 += tl2;
                sum3 += tl2 * tl;
                sum4 += tl2 * tl2;
                count += 1.0;
              }
            }
          }
          mt(i, j) = maxtltemp;

          // unprotected division by count but r=0,s=0 case guarantees count>=1
          c2(i, j) = sum2 / count;
          c3(i, j) = sum3 / count;
          c4(i, j) = sum4 / count;
        }
      }

      // scale the coeffs in Taylor series
      const double
        n = m_Glen_exponent,
        k  = (n + 2) / n,
        s2 = k * (2 * n + 2) / (2 * n),
        s3 = s2 * (3 * n + 2) / (3 * n),
        s4 = s3 * (4 * n + 2) / (4 * n);

      PetscErrorCode ierr;
      ierr = VecScale(*C2p0,s2);
      PISM_CHK(ierr, "VecScale");

      ierr = VecScale(*C3p0,s3);
      PISM_CHK(ierr, "VecScale");

      ierr = VecScale(*C4p0,s4);
      PISM_CHK(ierr, "VecScale");
    }
  } catch (...) {
    rank0.failed();
  }
  rank0.check();
}


//! Computes a smoothed thickness map.
/*!
The result `thksmooth` is the difference between the given upper surface
elevation (`usurf`) and the stored smoothed bed topography (`topgsmooth`),
except where the given original thickness (`thk`) is zero.  In places where
the original thickness is zero, the result `thksmooth` is also set to zero.

Ghosted values are updated directly and no communication occurs.  In fact,
we \e assume `usurf`, `thk`, and `thksmooth` all have stencil width at least
equal to GHOSTS.  We \e check whether `topgsmooth`, which has stencil width
maxGHOSTS, has at least GHOSTS stencil width, and throw an error if not.

Call preprocess_bed() first.
 */
void BedSmoother::get_smoothed_thk(const IceModelVec2S &usurf,
                                   const IceModelVec2S &thk,
                                   const IceModelVec2CellType &mask,
                                   IceModelVec2S &result) {

  IceModelVec::AccessList list;
  list.add(mask);
  list.add(maxtl);
  list.add(result);
  list.add(thk);
  list.add(topgsmooth);
  list.add(usurf);

  unsigned int GHOSTS = result.get_stencil_width();
  assert(mask.get_stencil_width()       >= GHOSTS);
  assert(maxtl.get_stencil_width()      >= GHOSTS);
  assert(thk.get_stencil_width()        >= GHOSTS);
  assert(topgsmooth.get_stencil_width() >= GHOSTS);
  assert(usurf.get_stencil_width()      >= GHOSTS);

  ParallelSection loop(grid->com);
  try {
    for (PointsWithGhosts p(*grid, GHOSTS); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (thk(i, j) < 0.0) {
        throw RuntimeError::formatted("BedSmoother detects negative original thickness\n"
                                      "at location (i, j) = (%d, %d) ... ending", i, j);
      } else if (thk(i, j) == 0.0) {
        result(i, j) = 0.0;
      } else if (maxtl(i, j) >= thk(i, j)) {
        result(i, j) = thk(i, j);
      } else {
        if (mask.grounded(i, j)) {
          // if grounded, compute smoothed thickness as the difference of ice
          // surface elevation and smoothed bed elevation
          const double thks_try = usurf(i, j) - topgsmooth(i, j);
          result(i, j) = (thks_try > 0.0) ? thks_try : 0.0;
        } else {
          // if floating, use original thickness (note: surface elevation was
          // computed using this thickness and the sea level elevation)
          result(i, j) = thk(i, j);
        }
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}


/*!
Implements the strategy for computing \f$\theta(h,x,y)\f$ from previously-
stored coefficients, described on [Bed roughness parameterization](@ref bedrough) page and in [\ref
Schoofbasaltopg2003].

Specifically, \f$\theta = \omega^{-n}\f$ where \f$\omega\f$ is a local average
of a rational function of surface elevation, approximated here by a Taylor polynomial:
  \f[ \omega = \fint \left(1 - \frac{\tilde b(x_1,x_2,\xi_1,\xi_2)}{H}
                           \right)^{-(n+2)/n}\,d\xi_1\,d\xi_2
             \approx 1 + C_2 H^{-2} + C_3 H^{-3} + C_4 H^{-4} \f]
where \f$h =\f$ usurf, \f$H = h -\f$ topgsmooth and \f$\tilde b\f$ is the local
bed topography, a function with mean zero.  The coefficients \f$C_2,C_3,C_4\f$,
which depend on \f$x,y\f$, are precomputed by `preprocess_bed()`.

Ghosted values are updated directly and no communication occurs.  In fact,
we \e assume `usurf` and `theta` have stencil width at least
equal to GHOSTS.  We \e check whether `topgsmooth`, which has stencil width
maxGHOSTS, has at least GHOSTS stencil width, and throw an error if not.

Call preprocess_bed() first.
 */
void BedSmoother::get_theta(const IceModelVec2S &usurf, IceModelVec2S &result) {

  if ((Nx < 0) || (Ny < 0)) {
    result.set(1.0);
    return;
  }

  IceModelVec::AccessList list;
  list.add(C2);
  list.add(C3);
  list.add(C4);
  list.add(maxtl);
  list.add(result);
  list.add(topgsmooth);
  list.add(usurf);

  unsigned int GHOSTS = result.get_stencil_width();
  assert(C2.get_stencil_width()         >= GHOSTS);
  assert(C3.get_stencil_width()         >= GHOSTS);
  assert(C4.get_stencil_width()         >= GHOSTS);
  assert(maxtl.get_stencil_width()      >= GHOSTS);
  assert(topgsmooth.get_stencil_width() >= GHOSTS);
  assert(usurf.get_stencil_width()      >= GHOSTS);

  ParallelSection loop(grid->com);
  try {
    for (PointsWithGhosts p(*grid, GHOSTS); p; p.next()) {
      const int i = p.i(), j = p.j();

      const double H = usurf(i, j) - topgsmooth(i, j);
      if (H > maxtl(i, j)) {
        // thickness exceeds maximum variation in patch of local topography,
        // so ice buries local topography; note maxtl >= 0 always
        const double Hinv = 1.0 / std::max(H, 1.0);
        double omega = 1.0 + Hinv*Hinv * (C2(i, j) + Hinv * (C3(i, j) + Hinv*C4(i, j)));
        if (omega <= 0) {  // this check *should not* be necessary: p4(s) > 0
          throw RuntimeError::formatted("omega is negative for i=%d, j=%d\n"
                                        "in BedSmoother.get_theta() ... ending", i, j);
        }

        if (omega < 0.001) {      // this check *should not* be necessary
          omega = 0.001;
        }

        result(i, j) = pow(omega,-m_Glen_exponent);
        // now guarantee in [0,1]; this check *should not* be necessary, by convexity of p4
        if (result(i, j) > 1.0) {
          result(i, j) = 1.0;
        }
        if (result(i, j) < 0.0) {
          result(i, j) = 0.0;
        }
      } else {
        result(i, j) = 0.00;  // FIXME = min_theta; make configurable
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}


} // end of namespace stressbalance
} // end of namespace pism
