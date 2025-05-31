
#include "examples.h"

#include "cea.h"


void run_example_02() {

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
  prob.set_problem("tv", "Example 2");

  // void CEA::Problem::set_output_options
  //   - SI:              bool (optional, default: true)
  //   - debug_points:    std::vector<size_t> (optional, default: {})
  //   - mass_fractions:  bool (optional, default: false)
  //   - _short:          bool (optional, default: false)
  //   - trace_tol:       double (optional, default: 0.0)
  //   - transport:       bool (optional, default: false)
  prob.set_output_options(false, {}, false, false, 0., true);

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
  prob.add_reactant("fuel", "H2",  std::nullopt, 100.);
  prob.add_reactant("oxyd", "Air", std::nullopt, 100.);

  prob.set_chamber_temperatures({3000});
  prob.set_chamber_densities({9.1864e-5, 8.0877e-6, 6.6054e-7}, "g/cc");
  prob.set_mixture_ratios({1.0}, "eq.ratio");

  prob.set_only_species({"Ar",   "C",    "CO",   "CO2",  "H",    "H2",   "H2O",  "HNO",  "HO2",  "HNO2", "HNO3",
                         "N",    "NH",   "NO",   "N2",   "N2O3", "O",    "O2",   "OH",   "O3"});

  prob.set_legacy_mode(true);

  prob.run("example-02.out");

}
