
#include "cea.h"


void run_legacy_examples() {
  std::vector<CEA::Problem> cea;

  cea = CEA::read_legacy_input("cea2.inp");
  CEA::run_all_cases(cea, "cea2.out", "cea2.plt");
  cea.clear();

  cea = CEA::read_legacy_input("test-thermp.inp");
  CEA::run_all_cases(cea, "test-thermp.out", "test-thermp.plt");
  cea.clear();

  cea = CEA::read_legacy_input("test-deton.inp");
  CEA::run_all_cases(cea, "test-deton.out", "test-deton.plt");
  cea.clear();

  cea = CEA::read_legacy_input("test-shock.inp");
  CEA::run_all_cases(cea, "test-shock.out");
  cea.clear();

  cea = CEA::read_legacy_input("test-rocket.inp");
  CEA::run_all_cases(cea, "test-rocket.out", "test-rocket.plt");
  cea.clear();
}
