
#include "examples.h"

#include "cea.h"


void run_example_14() {

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
  prob.set_problem("tp", "Example 14", true);

  // void CEA::Problem::set_output_options
  //   - SI:              bool (optional, default: true)
  //   - debug_points:    std::vector<size_t> (optional, default: {})
  //   - mass_fractions:  bool (optional, default: false)
  //   - _short:          bool (optional, default: false)
  //   - trace_tol:       double (optional, default: 0.0)
  //   - transport:       bool (optional, default: false)
  prob.set_output_options(true, {4});

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
  prob.add_reactant("name", "H2(L)", "H 2", 100., std::nullopt,  874.);
  prob.add_reactant("name", "O2(L)", "O 2",  60., std::nullopt, 1431.);

  prob.set_chamber_pressures({0.05}, "atm");
  prob.set_chamber_temperatures({1000., 500., 350., 305., 304.3, 304.2, 304., 300.});

  prob.set_legacy_mode(true);

  prob.run("example-14.out");

}
