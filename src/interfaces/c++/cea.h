#pragma once

#include <optional>
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

    friend void run_all_cases(std::vector<Problem>&, const std::optional<const char*>&, const std::optional<const char*>&);
    friend bool ffi_need_aligned(const std::vector<Problem>&);

  public:
    void set_problem(const char* mode, const std::optional<const char*>& name = std::nullopt, const bool& mole_ratios = false, const bool& equilibrium = false, const bool& ions = false,
                     const bool& frozen = false, const bool& frozen_at_throat = false, const std::optional<const char*>& thermo_lib = std::nullopt, const std::optional<const char*>& trans_lib = std::nullopt);
    void set_output_options(const bool& SI = true, const std::vector<size_t>& debug_points = {}, const bool& mass_fractions = false, const bool& _short = false,
                            const double& trace_tol = 0, const bool& transport = false);
    void set_chamber_pressures(const std::vector<double>& pressure_list, const std::optional<const char*>& unit = std::nullopt);
    void set_mixture_ratios(const std::vector<double>& ratio_list, const std::optional<const char*>& type = std::nullopt);
    void set_pressure_ratios(const std::vector<double>& ratio_list);
    void set_subsonic_area_ratios(const std::vector<double>& ratio_list);
    void set_supersonic_area_ratios(const std::vector<double>& ratio_list);
    void set_finite_area_combustor(const std::optional<double>& contraction_ratio = std::nullopt, const std::optional<double>& mass_flow_ratio = std::nullopt);
    void add_reactant(const char* type, const char* name, const std::optional<const char*>& formula = std::nullopt, const std::optional<double>& ratio = std::nullopt,
                      const std::optional<double>& T = std::nullopt, const std::optional<double>& rho = std::nullopt,
                      const std::optional<double>& h = std::nullopt, const std::optional<double>& u = std::nullopt,
                      const std::optional<const char*>& T_unit = std::nullopt, const std::optional<const char*>& rho_unit = std::nullopt,
                      const std::optional<const char*>& h_unit = std::nullopt, const std::optional<const char*> & u_unit = std::nullopt);
    void set_insert_species(const std::vector<std::string>&);
    void set_omit_species(const std::vector<std::string>&);
    void set_only_species(const std::vector<std::string>&);
    void set_legacy_mode(const bool& legacy_mode = true);
    void run(const std::optional<const char*>& out_filename = std::nullopt);
    void write_debug_output(const char* out_filename);
  };


  std::vector<Problem> read_legacy_input(const char*);
  void run_all_cases(std::vector<Problem>&, const std::optional<const char*>& out_filename = std::nullopt, const std::optional<const char*>& plt_filename = std::nullopt);

  bool ffi_need_aligned(const std::vector<Problem>&);
}
