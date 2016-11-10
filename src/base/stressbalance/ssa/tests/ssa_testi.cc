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

static char help[] =
  "\nSSA_TESTI\n"
  "  Testing program for the finite element implementation of the SSA.\n"
  "  Does a time-independent calculation.  Does not run IceModel or a derived\n"
  "  class thereof. Uses verification test I. Also may be used in a PISM\n"
  "  software (regression) test.\n\n";

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

class SSATestCaseI: public SSATestCase {
public:
  SSATestCaseI(Context::Ptr ctx)
    : SSATestCase(ctx) {
    // empty
  }

protected:
  virtual void initializeGrid(int Mx,int My);

  virtual void initializeSSAModel();

  virtual void initializeSSACoefficients();

  virtual void exactSolution(int i, int j,
    double x, double y, double *u, double *v);

};

const double m_schoof = 10; // (pure number)
const double L_schoof = 40e3; // meters
const double aspect_schoof = 0.05; // (pure)
const double H0_schoof = aspect_schoof * L_schoof;
                                       // = 2000 m THICKNESS
const double B_schoof = 3.7e8; // Pa s^{1/3}; hardness
                                     // given on p. 239 of Schoof; why so big?

void SSATestCaseI::initializeGrid(int Mx,int My) {
  double Ly = 3*L_schoof;  // 300.0 km half-width (L=40.0km in Schoof's choice of variables)
  double Lx = std::max(60.0e3, ((Mx - 1) / 2) * (2.0 * Ly / (My - 1)));
  m_grid = IceGrid::Shallow(m_ctx, Lx, Ly,
                            0.0, 0.0, // center: (x0,y0)
                            Mx, My, NOT_PERIODIC);
}


void SSATestCaseI::initializeSSAModel() {
  m_enthalpyconverter = EnthalpyConverter::Ptr(new EnthalpyConverter(*m_config));

  m_config->set_boolean("basal_resistance.pseudo_plastic.enabled", false);

  m_config->set_string("stress_balance.ssa.flow_law", "isothermal_glen");
  m_config->set_double("flow_law.isothermal_Glen.ice_softness", pow(B_schoof, -m_config->get_double("stress_balance.ssa.Glen_exponent")));
}

void SSATestCaseI::initializeSSACoefficients() {

  m_bc_mask.set(0);
  m_thickness.set(H0_schoof);

  // ssa->strength_extension->set_min_thickness(2*H0_schoof);

  // The finite difference code uses the following flag to treat the non-periodic grid correctly.
  m_config->set_boolean("stress_balance.ssa.compute_surface_gradient_inward", true);
  m_config->set_double("stress_balance.ssa.epsilon", 0.0);  // don't use this lower bound

  IceModelVec::AccessList list;
  list.add(m_tauc);

  double standard_gravity = m_config->get_double("constants.standard_gravity"),
    ice_rho = m_config->get_double("constants.ice.density");

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double y = m_grid->y(j);
    const double theta = atan(0.001);   /* a slope of 1/1000, a la Siple streams */
    const double f = ice_rho * standard_gravity * H0_schoof * tan(theta);
    m_tauc(i,j) = f * pow(fabs(y / L_schoof), m_schoof);
  }
  m_tauc.update_ghosts();

  list.add(m_bc_values);
  list.add(m_bc_mask);
  list.add(m_surface);
  list.add(m_bed);
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double myx = m_grid->x(i), myy=m_grid->y(j);
    // eval exact solution; will only use exact vels if at edge
    struct TestIParameters I_parameters = exactI(m_schoof, myx, myy);
    m_bed(i, j) = I_parameters.bed;
    m_surface(i,j) = m_bed(i,j) + H0_schoof;

    bool edge = ((j == 0) || (j == (int)m_grid->My() - 1) ||
                 (i == 0) || (i == (int)m_grid->Mx() - 1));
    if (edge) {
      m_bc_mask(i,j) = 1;
      m_bc_values(i,j).u = I_parameters.u;
      m_bc_values(i,j).v = I_parameters.v;
    }
  }

  // communicate what we have set
  m_surface.update_ghosts();
  m_bed.update_ghosts();
  m_bc_mask.update_ghosts();
  m_bc_values.update_ghosts();
}


void SSATestCaseI::exactSolution(int /*i*/, int /*j*/,
                                 double x, double y,
                                 double *u, double *v) {
  struct TestIParameters I_parameters = exactI(m_schoof, x, y);
  *u = I_parameters.u;
  *v = I_parameters.v;
}

} // end of namespace stressbalance
} // end of namespace pism

int main(int argc, char *argv[]) {

  using namespace pism;
  using namespace pism::stressbalance;

  MPI_Comm com = MPI_COMM_WORLD;
  petsc::Initializer petsc(argc, argv, help);
  PetscErrorCode ierr;

  com = PETSC_COMM_WORLD;

  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  try {
    Context::Ptr ctx = context_from_options(com, "ssa_testi");
    Config::Ptr config = ctx->config();

    bool
      usage_set = options::Bool("-usage", "print usage info"),
      help_set  = options::Bool("-help", "print help info");
    if (usage_set or help_set) {
      ierr = PetscPrintf(com,
                         "\n"
                         "usage of SSA_TESTi:\n"
                         "  run ssa_testi -Mx <number> -My <number> -ssa_method <fd|fem>\n"
                         "\n");
      PISM_CHK(ierr, "PetscPrintf");
    }

    // Parameters that can be overridden by command line options
    unsigned int Mx = config->get_double("grid.Mx");
    unsigned int My = config->get_double("grid.My");
    options::Keyword method("-ssa_method", "Algorithm for computing the SSA solution",
                            "fem,fd", "fem");

    options::String output_file("-o", "Set the output file name", "ssa_test_i.nc");

    // Determine the kind of solver to use.
    SSAFactory ssafactory = NULL;
    if (method == "fem") {
      ssafactory = SSAFEMFactory;
    } else if (method == "fd") {
      ssafactory = SSAFDFactory;
    } else {
      /* can't happen */
    }

    SSATestCaseI testcase(ctx);
    testcase.init(Mx,My,ssafactory);
    testcase.run();
    testcase.report("I");
    testcase.write(output_file);
  }
  catch (...) {
    handle_fatal_errors(com);
    return 1;
  }

  return 0;
}
