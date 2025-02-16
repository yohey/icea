
#include "examples.h"

#include "cea.h"


void run_example_12() {

  CEA::Problem prob;

  // void CEA::Problem::set_problem
  //   - mode:              std::string (required)
  //   - name:              std::string (required)
  //   - mole_ratios:       bool (optional, default: false)
  //   - equilibrium:       bool (optional, default: false)
  //   - ions:              bool (optional, default: false)
  //   - frozen:            bool (optional, default: false)
  //   - frozen_at_throat:  bool (optional, default: false)
  //   - thermo_lib:        std::string (optional, default: auto detect)
  //   - trans_lib:         std::string (optional, default: auto detect)
  prob.set_problem("rocket", "Example 12", false, true, false, true, true);

  // void CEA::Problem::set_output_options
  //   - SI:              bool (optional, default: true)
  //   - debug_points:    std::vector<size_t> (optional, default: {})
  //   - mass_fractions:  bool (optional, default: false)
  //   - _short:          bool (optional, default: false)
  //   - trace_tol:       double (optional, default: 0.0)
  //   - transport:       bool (optional, default: false)
  prob.set_output_options(true, {}, true);

  // void CEA::Problem::add_reactant
  //   - type:            std::string (required, must start with "fu", "ox" or "na")
  //   - name:            std::string (required)
  //   - ratio:           double (optional)
  //   - T:               double (optional)
  //   - rho:             double (optional)
  //   - T_unit:          std::string (optional, default: K)
  //   - rho_unit:        std::string (optional, default: g/cc)
  prob.add_reactant("fuel", "CH6N2(L)", std::nullopt, std::nullopt, 0.874);
  prob.add_reactant("oxyd", "N2O4(L)",  std::nullopt, std::nullopt, 1.431);

  prob.set_chamber_pressures({1000}, "legacy-psi");
  prob.set_mixture_ratios({2.5});
  prob.set_pressure_ratios({68.0457});
  prob.set_supersonic_area_ratios({5, 10, 25, 50, 75, 100, 150, 200});

  // prob.set_only_species({"CO",     "CO2",    "H",      "HNO",    "HNO2",   "HO2",    "H2",     "H2O",    "H2O2",
  //                        "N",      "NO",     "NO2",    "N2",     "N2O",    "O",      "OH",     "O2",     "HCO",
  //                        "NH",     "CH4",    "NH2",    "NH3",    "H2O(L)", "C(gr)"});

  prob.set_legacy_mode(true);

  prob.run("example-12.out");

}
