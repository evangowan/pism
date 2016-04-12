// Copyright (C) 2010--2016 PISM Authors
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

static char help[] = "\nBEDROUGH_TEST\n"
  "  Simple testing program for Schoof (2003)-type bed smoothing and roughness-\n"
  "  parameterization schemes.  Allows comparison of computed theta to result\n"
  "  from Matlab/Octave code exampletheta.m in src/base/bedroughplay.  Also\n"
  "  used in PISM software (regression) test.\n\n";

#include "base/util/Context.hh"
#include <cmath>
#include <gsl/gsl_math.h>       // M_PI
#include <cstdio>
#include "base/util/pism_options.hh"
#include "base/util/IceGrid.hh"
#include "base/util/iceModelVec.hh"
#include "base/stressbalance/sia/PISMBedSmoother.hh"

#include "base/util/petscwrappers/PetscInitializer.hh"
#include "base/util/error_handling.hh"

using namespace pism;

int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;
  MPI_Comm com = MPI_COMM_WORLD;

  petsc::Initializer petsc(argc, argv, help);

  com = PETSC_COMM_WORLD;

  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  try {
    Context::Ptr ctx = context_from_options(com, "bedrough_test");
    Config::Ptr config = ctx->config();

    GridParameters P(config);

    P.Lx = 1200e3;
    P.Ly = P.Lx;
    P.Mx = 81;
    P.My = P.Mx;
    P.vertical_grid_from_options(config);
    P.ownership_ranges_from_options(ctx->size());
    P.periodicity = NOT_PERIODIC;

    // create grid
    IceGrid::Ptr grid(new IceGrid(ctx, P));

    ierr = PetscPrintf(grid->com,"BedSmoother TEST\n");
    PISM_CHK(ierr, "PetscPrintf");

    bool show = options::Bool("-show", "turn on diagnostic viewers");

    IceModelVec2S topg, usurf, theta;
    topg.create(grid, "topg", WITH_GHOSTS, 1);
    topg.set_attrs("trybedrough_tool", "original topography",
                   "m", "bedrock_altitude");
    usurf.create(grid, "usurf", WITH_GHOSTS, 1);
    usurf.set_attrs("trybedrough_tool", "ice surface elevation",
                    "m", "surface_altitude");
    theta.create(grid, "theta", WITH_GHOSTS, 1);
    theta.set_attrs("trybedrough_tool",
                    "coefficient theta in Schoof (2003) bed roughness parameterization",
                    "", "");

    // put in bed elevations, a la this Matlab:
    //    topg0 = 400 * sin(2 * pi * xx / 600e3) + ...
    //            100 * sin(2 * pi * (xx + 1.5 * yy) / 40e3);
    IceModelVec::AccessList list(topg);
    for (Points p(*grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      topg(i,j) = 400.0 * sin(2.0 * M_PI * grid->x(i) / 600.0e3) +
        100.0 * sin(2.0 * M_PI * (grid->x(i) + 1.5 * grid->y(j)) / 40.0e3);
    }

    usurf.set(1000.0);  // compute theta for this constant thk

    // actually use the smoother/bed-roughness-parameterizer
    config->set_double("Glen_exponent", 3.0);
    config->set_double("bed_smoother_range", 50.0e3);
    stressbalance::BedSmoother smoother(grid, 1);
    smoother.preprocess_bed(topg);
    int Nx,Ny;
    smoother.get_smoothing_domain(Nx,Ny);

    ierr = PetscPrintf(grid->com,"  smoothing domain:  Nx = %d, Ny = %d\n",Nx,Ny);
    PISM_CHK(ierr, "PetscPrintf");

    smoother.get_theta(usurf, theta);

    const IceModelVec2S &topg_smoothed = smoother.get_smoothed_bed();
    if (show) {
      const int  window = 400;
      topg.view(window);
      topg_smoothed.view(window);
      theta.view(window);
      printf("[showing topg, topg_smoothed, theta in X windows for 10 seconds ...]\n");

      ierr = PetscSleep(10);
      PISM_CHK(ierr, "PetscSleep");
    }

    double topg_min, topg_max, topgs_min, topgs_max, theta_min, theta_max;
    topg_min = topg.min();
    topg_max = topg.max();
    topgs_min = topg_smoothed.min();
    topgs_max = topg_smoothed.max();
    theta_min = theta.min();
    theta_max = theta.max();
    ierr = PetscPrintf(grid->com,
                       "  original bed    :  min elev = %12.6f m,  max elev = %12.6f m\n",
                       topg_min, topg_max);
    PISM_CHK(ierr, "PetscPrintf");

    ierr = PetscPrintf(grid->com,
                       "  smoothed bed    :  min elev = %12.6f m,  max elev = %12.6f m\n",
                       topgs_min, topgs_max);
    PISM_CHK(ierr, "PetscPrintf");

    ierr = PetscPrintf(grid->com,
                       "  Schoof's theta  :  min      = %12.9f,    max      = %12.9f\n",
                       theta_min, theta_max);
    PISM_CHK(ierr, "PetscPrintf");

    bool dump = options::Bool("-dump", "dump bed roughness data");
    if (dump) {
      topg.dump("bedrough_test_topg.nc");
      topg_smoothed.dump("bedrough_test_topg_smoothed.nc");
      theta.dump("bedrough_test_theta.nc");
    }

  }
  catch (...) {
    handle_fatal_errors(com);
    return 1;
  }
  return 0;
}
