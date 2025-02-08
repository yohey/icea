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

    friend void run_all_cases(std::vector<Problem>&, const std::string&, const std::string&);
    friend bool ffi_need_aligned(const std::vector<Problem>&);

  public:
    void print_case();
  };


  std::vector<Problem> read_legacy_input(const std::string&);
  void run_all_cases(std::vector<Problem>&, const std::string& out_filename = "", const std::string& plt_filename = "");

  bool ffi_need_aligned(const std::vector<Problem>&);
}
