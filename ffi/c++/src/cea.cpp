
#include "cea.h"


namespace CEA {

  Problem::Problem() {
    this->_ffi = ffi_cea_new_problem();
  }

  Problem::Problem(FFI_CEA_Problem_Ptr ptr) {
    this->_ffi = ptr;
  }

  Problem::~Problem() {
    ffi_cea_del_problem(this->_ffi);
  }


  void Problem::print_case() {
    ffi_cea_print_case(this->_ffi);
  }


  std::vector<Problem> read_legacy_input(const std::string& filename) {
    FFI_CEA_Problem_Array array = ffi_cea_read_legacy_input(filename.c_str());

    std::vector<Problem> problems;
    problems.reserve(array.size);

    for (int i = 0; i < array.size; ++i) {
      problems.emplace_back(array.addr[i]);
    }

    return problems;
  }

}
