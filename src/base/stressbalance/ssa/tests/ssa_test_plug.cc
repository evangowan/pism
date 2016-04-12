// Copyright (C) 2010--2016 Ed Bueler, Constantine Khroulev, and David Maxwell
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

/* This file implements a test case for the ssa: plug flow. The geometry
   consists of a constant surface slope in the positive x-direction, and the
   ice is pinned on the y-boundaries. There is no basal shear stress, and hence
   the the only nonzero terms in the SSA are the "p-laplacian" and the driving
   stress.
*/

static char help[] =
  "\nSSA_TEST_PLUG\n"
  "  Testing program for the finite element implementation of the SSA.\n"
  "  Does a time-independent calculation.  Does not run IceModel or a derived\n"
  "  class thereof.\n\n";

#include <cmath>

#include "base/basalstrength/basal_resistance.hh" // IceBasalResistancePlasticLaw
#include "base/stressbalance/ssa/SSAFD.hh"
#include "base/stressbalance/ssa/SSAFEM.hh"
#include "base/stressbalance/ssa/SSATestCase.hh"
#include "base/util/Context.hh"
#include "base/util/VariableMetadata.hh"
#include "base/util/error_handling.hh"
#include "base/util/iceModelVec.hh"
#include "base/util/io/PIO.hh"
#include "base/util/petscwrappers/PetscInitializer.hh"
#include "base/util/pism_const.hh"
#include "base/util/pism_options.hh"
#include "verif/tests/exactTestsIJ.h"

namespace pism {
namespace stressbalance {

class SSATestCasePlug: public SSATestCase {
public:
  SSATestCasePlug(Context::Ptr ctx, double n)
    : SSATestCase(ctx) {
    H0    = 2000.;              //m
    L     = 50.e3;              // 50km half-width
    dhdx  = 0.001;              // pure number, slope of surface & bed
    tauc0 = 0.;                 // No basal shear stress

    // Pa s^{1/3}; hardness given on p. 239 of Schoof; why so big?
    B0           = 3.7e8;

    this->glen_n = n;
  }

protected:
  virtual void initializeGrid(int Mx,int My);

  virtual void initializeSSAModel();

  virtual void initializeSSACoefficients();

  virtual void exactSolution(int i, int j,
    double x, double y, double *u, double *v);


  double H0; // Thickness
  double L;  // Half-width
  double dhdx; // surface slope
  double tauc0; // zero basal shear stress
  double B0;  // hardness
  double glen_n;

  bool dimensionless;

};


void SSATestCasePlug::initializeGrid(int Mx,int My) {
  double Lx=L, Ly = L;
  m_grid = IceGrid::Shallow(m_ctx, Lx, Ly,
                            0.0, 0.0, // center: (x0,y0)
                            Mx, My, NOT_PERIODIC);
}


void SSATestCasePlug::initializeSSAModel() {
  // Basal sliding law parameters are irrelevant because tauc=0

  // Enthalpy converter is irrelevant (but still required) for this test.
  m_enthalpyconverter = EnthalpyConverter::Ptr(new EnthalpyConverter(*m_config));

  // Use constant hardness
  m_config->set_string("ssa_flow_law", "isothermal_glen");
  m_config->set_double("ice_softness", pow(B0, -glen_n));
}

void SSATestCasePlug::initializeSSACoefficients() {

  // The finite difference code uses the following flag to treat the non-periodic grid correctly.
  m_config->set_boolean("compute_surf_grad_inward_ssa", true);
  m_config->set_double("epsilon_ssa", 0.0);

  // Ensure we never use the strength extension.
  m_ssa->strength_extension->set_min_thickness(H0/2);

  // Set constant coefficients.
  m_thickness.set(H0);
  m_tauc.set(tauc0);


  // Set boundary conditions (Dirichlet all the way around).
  m_bc_mask.set(0.0);

  IceModelVec::AccessList list;
  list.add(m_bc_values);
  list.add(m_bc_mask);
  list.add(m_bed);
  list.add(m_surface);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double myu, myv;
    const double myx = m_grid->x(i), myy=m_grid->y(j);

    m_bed(i,j) = -myx*(dhdx);
    m_surface(i,j) = m_bed(i,j) + H0;

    bool edge = ((j == 0) || (j == (int)m_grid->My() - 1) ||
                 (i == 0) || (i == (int)m_grid->Mx() - 1));
    if (edge) {
      m_bc_mask(i,j) = 1;
      exactSolution(i,j,myx,myy,&myu,&myv);
      m_bc_values(i,j).u = myu;
      m_bc_values(i,j).v = myv;
    }
  }

  m_bc_values.update_ghosts();
  m_bc_mask.update_ghosts();
  m_bed.update_ghosts();
  m_surface.update_ghosts();

  m_ssa->set_boundary_conditions(m_bc_mask, m_bc_values);
}

void SSATestCasePlug::exactSolution(int /*i*/, int /*j*/,
                                              double /*x*/, double y,
                                              double *u, double *v) {
  double earth_grav = m_config->get_double("standard_gravity"),
    ice_rho = m_config->get_double("ice_density");
  double f = ice_rho * earth_grav * H0* dhdx;
  double ynd = y/L;

  *u = 0.5*pow(f,3)*pow(L,4)/pow(B0*H0,3)*(1-pow(ynd,4));
  *v = 0;
}

} // end of namespace stressbalance
} // end of namespace pism

int main(int argc, char *argv[]) {

  using namespace pism;
  using namespace pism::stressbalance;

  MPI_Comm com = MPI_COMM_WORLD;  // won't be used except for rank,size
  petsc::Initializer petsc(argc, argv, help);
  PetscErrorCode ierr;

  com = PETSC_COMM_WORLD;

  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  try {
    verbosityLevelFromOptions();
    Context::Ptr ctx = context_from_options(com, "ssa_test_plug");
    Config::Ptr config = ctx->config();

    setVerbosityLevel(5);

    bool
      usage_set = options::Bool("-usage", "show the usage message"),
      help_set  = options::Bool("-help", "show the help message");
    if (usage_set or help_set) {
      ierr = PetscPrintf(com,
                         "\n"
                         "usage of SSA_TEST_PLUG:\n"
                         "  run ssa_test_plug -Mx <number> -My <number> -ssa_method <fd|fem>\n"
                         "\n");
      PISM_CHK(ierr, "PetscPrintf");
    }

    // Parameters that can be overridden by command line options

    options::Integer Mx("-Mx", "Number of grid points in the X direction", 11);
    options::Integer My("-My", "Number of grid points in the Y direction", 61);

    options::Keyword method("-ssa_method", "Algorithm for computing the SSA solution",
                            "fem,fd", "fem");

    options::String output_file("-o", "Set the output file name", "ssa_test_plug.nc");
    options::Real glen_n("-ssa_glen_n", "Glen exponent for the SSA", 3.0);

    options::Integer my_verbosity_level("-verbose", "Verbosity level", 2);
    if (my_verbosity_level.is_set()) {
      setVerbosityLevel(my_verbosity_level);
    }

    // Determine the kind of solver to use.
    SSAFactory ssafactory = NULL;
    if (method == "fem") {
      ssafactory = SSAFEMFactory;
    } else if (method == "fd") {
      ssafactory = SSAFDFactory;
    } else {
      /* can't happen */
    }

    SSATestCasePlug testcase(ctx, glen_n);
    testcase.init(Mx,My,ssafactory);
    testcase.run();
    testcase.report("plug");
    testcase.write(output_file);
  }
  catch (...) {
    handle_fatal_errors(com);
    return 1;
  }

  return 0;
}
