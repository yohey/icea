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

    // friend void run_all_cases(std::vector<Problem>&, const std::optional<std::string>&, const std::optional<std::string>&);
    // friend bool ffi_need_aligned(const std::vector<Problem>&);

  public:
    void set_problem(const std::string& mode, const std::optional<std::string>& name = std::nullopt, const bool& mole_ratios = false, const bool& equilibrium = false, const bool& ions = false,
                     const bool& frozen = false, const bool& frozen_at_throat = false, const std::optional<std::string>& thermo_lib = std::nullopt, const std::optional<std::string>& trans_lib = std::nullopt);
    void set_output_options(const bool& SI = true, const std::vector<size_t>& debug_points = {}, const bool& mass_fractions = false, const bool& _short = false,
                            const double& trace_tol = 0, const bool& transport = false);
    void set_chamber_pressures(const std::vector<double>& pressure_list, const std::optional<std::string>& unit = std::nullopt);
    void set_mixture_ratios(const std::vector<double>& ratio_list, const std::optional<std::string>& type = std::nullopt);
    void set_pressure_ratios(const std::vector<double>& ratio_list);
    void set_subsonic_area_ratios(const std::vector<double>& ratio_list);
    void set_supersonic_area_ratios(const std::vector<double>& ratio_list);
    void set_finite_area_combustor(const std::optional<double>& contraction_ratio = std::nullopt, const std::optional<double>& mass_flow_ratio = std::nullopt);
    void add_reactant(const std::string& type, const std::string& name, const std::optional<double>& ratio = std::nullopt, const std::optional<double>& T = std::nullopt,
                      const std::optional<double>& rho = std::nullopt, const std::optional<std::string>& T_unit = std::nullopt, const std::optional<std::string>& rho_unit = std::nullopt);
    void insert_species(const std::vector<std::string>&);
    void set_omit_species(const std::vector<std::string>&);
    void set_only_species(const std::vector<std::string>&);
    void set_legacy_mode(const bool& legacy_mode = true);
    void run(const std::optional<std::string>& out_filename = std::nullopt);
    void write_debug_output(const std::string& out_filename);
  };


  // std::vector<Problem> read_legacy_input(const std::string&);
  // void run_all_cases(std::vector<Problem>&, const std::optional<std::string>& out_filename = std::nullopt, const std::optional<std::string>& plt_filename = std::nullopt);

  // bool ffi_need_aligned(const std::vector<Problem>&);
}
