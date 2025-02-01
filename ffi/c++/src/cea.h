#pragma once

#include <string>
#include <vector>

#include "mod_ffi.h"


namespace CEA {

  class Problem {
  public:
    Problem();
    Problem(FFI_CEA_Problem_Ptr);
    ~Problem();

  private:
    FFI_CEA_Problem_Ptr _ffi;

  public:
    void print_case();
  };


  std::vector<Problem> read_legacy_input(const std::string&);
}
