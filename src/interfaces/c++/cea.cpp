
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


  void run_all_cases(std::vector<Problem>& problems, const std::string& out_filename, const std::string& plt_filename) {
    FFI_CEA_Problem_Array array;
    array.size = problems.size();

    bool _need_aligned_buffer = ffi_need_aligned(problems);
    size_t obj_size = ffi_cea_sizeof(problems[0]._ffi);

    if (_need_aligned_buffer) {
      array.addr = static_cast<FFI_CEA_Problem_Ptr*>(malloc(problems.size() * obj_size));

      for (size_t i = 0; i < problems.size(); ++i) {
        FFI_CEA_Problem_Ptr ptr = reinterpret_cast<FFI_CEA_Problem_Ptr>((reinterpret_cast<uintptr_t>(array.addr) + i * obj_size));
        memcpy(ptr, problems[i]._ffi, obj_size);
      }

    } else {
      array.addr = static_cast<FFI_CEA_Problem_Ptr*>(problems[0]._ffi);
    }

    if (out_filename.length() > 0 && plt_filename.length() > 0) {
      ffi_cea_run_all_cases(array, out_filename.c_str(), plt_filename.c_str());
    } else if (out_filename.length() > 0) {
      ffi_cea_run_all_cases(array, out_filename.c_str());
    } else if (plt_filename.length() > 0) {
      ffi_cea_run_all_cases(array, nullptr, plt_filename.c_str());
    } else {
      ffi_cea_run_all_cases(array);
    }

    if (_need_aligned_buffer) {
      for (size_t i = 0; i < problems.size(); ++i) {
        FFI_CEA_Problem_Ptr ptr = reinterpret_cast<FFI_CEA_Problem_Ptr>((reinterpret_cast<uintptr_t>(array.addr) + i * obj_size));
        memcpy(problems[i]._ffi, ptr, obj_size);
      }

      free(array.addr);
    }
  }


  bool ffi_need_aligned(const std::vector<Problem>& problems) {
    if (problems.size() < 2) return false;

    size_t obj_size = ffi_cea_sizeof(problems[0]._ffi);
    uintptr_t head_addr = reinterpret_cast<uintptr_t>(problems[0]._ffi);
    uintptr_t prev_tail_addr = head_addr + obj_size;

    for (size_t i = 1; i < problems.size(); ++i) {
      head_addr = reinterpret_cast<uintptr_t>(problems[i]._ffi);
      if (head_addr != prev_tail_addr) {
        return true;
      } else {
        prev_tail_addr = head_addr + obj_size;
      }
    }

    return false;
  }
}
