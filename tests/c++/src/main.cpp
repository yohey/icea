
#include "cea.h"


int main() {

  // CEA::Problem cea;
  // cea.print_case();

  std::vector<CEA::Problem> cea = CEA::read_legacy_input("../../../tests/orig/cea2.inp");

  return 0;
}
