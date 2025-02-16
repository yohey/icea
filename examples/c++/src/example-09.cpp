
#include "examples.h"

#include "cea.h"


void run_example_09() {

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
  prob.set_problem("rocket", "Example 9", false, false);

  // void CEA::Problem::set_output_options
  //   - SI:              bool (optional, default: true)
  //   - debug_points:    std::vector<size_t> (optional, default: {})
  //   - mass_fractions:  bool (optional, default: false)
  //   - _short:          bool (optional, default: false)
  //   - trace_tol:       double (optional, default: 0.0)
  //   - transport:       bool (optional, default: false)
  prob.set_output_options(true);

  // void CEA::Problem::add_reactant
  //   - type:            std::string (required, must start with "fu", "ox" or "na")
  //   - name:            std::string (required)
  //   - ratio:           double (optional)
  //   - T:               double (optional)
  //   - rho:             double (optional)
  //   - T_unit:          std::string (optional, default: K)
  //   - rho_unit:        std::string (optional, default: g/cc)
  prob.add_reactant("fuel", "H2(L)", 100, 20.27);
  prob.add_reactant("oxyd", "O2(L)", 100, 90.17);

  prob.set_finite_area_combustor(1.58);
  prob.set_chamber_pressures({53.3172});
  prob.set_mixture_ratios({5.55157});
  prob.set_pressure_ratios({10, 100, 1000});
  prob.set_supersonic_area_ratios({25, 50, 75});

  prob.set_legacy_mode(true);

  prob.run("example-09.out");

}
