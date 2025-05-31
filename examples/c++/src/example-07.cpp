
#include "examples.h"

#include "cea.h"


void run_example_07() {

  CEA::Problem prob;

  // void CEA::Problem::set_problem
  //   - mode:              std::string (required)
  //   - name:              std::string (required)
  //   - mole_ratios:       bool (optional, default: false)
  //   - equilibrium:       bool (optional, default: false)
  //   - ions:              bool (optional, default: false)
  //   - frozen:            bool (optional, default: false)
  //   - frozen_at_throat:  bool (optional, default: false)
  //   - incident:          bool (optional, default: false)
  //   - reflected:         bool (optional, default: false)
  //   - thermo_lib:        std::string (optional, default: auto detect)
  //   - trans_lib:         std::string (optional, default: auto detect)
  prob.set_problem("shock", "Example 7", true, true, false, true, false, true);

  // void CEA::Problem::add_reactant
  //   - type:            std::string (required, must start with "fu", "ox" or "na")
  //   - name:            std::string (required)
  //   - formula:         std::string (optional, required when h or u is given)
  //   - ratio:           double (optional)
  //   - T:               double (optional)
  //   - rho:             double (optional)
  //   - h:               double (optional)
  //   - u:               double (optional)
  //   - T_unit:          std::string (optional, default: K)
  //   - rho_unit:        std::string (optional, default: kg/m^3)
  //   - h_unit:          std::string (optional, default: J/mol)
  //   - u_unit:          std::string (optional, default: J/mol)
  prob.add_reactant("name", "H2", std::nullopt, 0.05, 300.);
  prob.add_reactant("name", "O2", std::nullopt, 0.05, 300.);
  prob.add_reactant("name", "Ar", std::nullopt, 0.90, 300.);

  prob.set_chamber_pressures({10.0, 20.0}, "mmHg");

  prob.set_initial_velocities({1000., 1100., 1200., 1250., 1300., 1350., 1400.});

  prob.set_legacy_mode(true);

  prob.run("example-07.out");

}
