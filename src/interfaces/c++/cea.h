#pragma once

#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "mod_ffi.h"


namespace CEA {

  class Problem {
  public:
    Problem();
    Problem(Problem&&);
    Problem(const Problem&);
    Problem(_FFI_CEA_Problem_Ptr);
    ~Problem() = default;

  private:
    std::unique_ptr<void, std::function<void(const _FFI_CEA_Problem_Ptr&)>> _ffi;

    friend void run_all_cases(std::vector<Problem>&, const std::optional<const char*>&, const std::optional<const char*>&);
    friend bool ffi_need_aligned(const std::vector<Problem>&);

  public:
    void set_problem(const char* mode, const std::optional<const char*>& name = std::nullopt, const bool& mole_ratios = false, const bool& equilibrium = false, const bool& ions = false,
                     const bool& frozen = false, const bool& frozen_at_throat = false, const bool& incident = false, const bool& reflected = false,
                     const std::optional<const char*>& thermo_lib = std::nullopt, const std::optional<const char*>& trans_lib = std::nullopt);
    void set_output_options(const bool& SI = true, const std::vector<size_t>& debug_points = {}, const bool& mass_fractions = false, const bool& _short = false,
                            const double& trace_tol = 0, const bool& transport = false, const std::vector<std::string>& plot = {});
    void set_chamber_pressures(const std::vector<double>& pressure_list, const std::optional<const char*>& unit = std::nullopt);
    void set_chamber_temperatures(const std::vector<double>& temperature_list, const std::optional<const char*>& unit = std::nullopt);
    void set_chamber_densities(const std::vector<double>& density_list, const std::optional<const char*>& unit = std::nullopt);
    void set_chamber_internal_energy(const double& uByR);
    void set_mixture_ratios(const std::vector<double>& ratio_list, const std::optional<const char*>& type = std::nullopt);
    void set_pressure_ratios(const std::vector<double>& ratio_list);
    void set_subsonic_area_ratios(const std::vector<double>& ratio_list);
    void set_supersonic_area_ratios(const std::vector<double>& ratio_list);
    void set_finite_area_combustor(const std::optional<double>& contraction_ratio = std::nullopt, const std::optional<double>& mass_flow_ratio = std::nullopt);
    void set_initial_velocities(const std::vector<double>& velocity_list, const bool& is_mach = false);
    void add_reactant(const char* type, const char* name, const std::optional<const char*>& formula = std::nullopt, const std::optional<double>& ratio = std::nullopt,
                      const std::optional<double>& T = std::nullopt, const std::optional<double>& rho = std::nullopt,
                      const std::optional<double>& h = std::nullopt, const std::optional<double>& u = std::nullopt,
                      const std::optional<const char*>& T_unit = std::nullopt, const std::optional<const char*>& rho_unit = std::nullopt,
                      const std::optional<const char*>& h_unit = std::nullopt, const std::optional<const char*> & u_unit = std::nullopt);
    void set_reactant(const size_t& index, const std::optional<double>& ratio = std::nullopt,
                      const std::optional<double>& T = std::nullopt, const std::optional<double>& rho = std::nullopt,
                      const std::optional<double>& h = std::nullopt, const std::optional<double>& u = std::nullopt,
                      const std::optional<const char*>& T_unit = std::nullopt, const std::optional<const char*>& rho_unit = std::nullopt,
                      const std::optional<const char*>& h_unit = std::nullopt, const std::optional<const char*> & u_unit = std::nullopt);
    void set_insert_species(const std::vector<std::string>&);
    void set_omit_species(const std::vector<std::string>&);
    void set_only_species(const std::vector<std::string>&);
    void set_legacy_mode(const bool& legacy_mode = true);
    void run(const std::optional<const char*>& out_filename = std::nullopt, const std::optional<const char*>& plt_filename = std::nullopt);

    double get_pressure(const size_t& iOF, const size_t& ipt);
    double get_temperature(const size_t& iOF, const size_t& ipt);
    double get_chamber_temperature();
    double get_molecular_weight(const size_t& iOF, const size_t& ipt);
    double get_specific_heat(const size_t& iOF, const size_t& ipt);
    double get_specific_heat_ratio(const size_t& iOF, const size_t& ipt);
    double get_characteristic_velocity();
    double get_specific_impulse(const size_t& iOF, const size_t& ipt, const bool& vacuum = false);

    void calc_frozen_exhaust(const double& T, const std::optional<const double*>& cp_ptr = std::nullopt, const std::optional<const double*>& gamma_ptr = std::nullopt,
                             const std::optional<const double*>& mu_ptr = std::nullopt, const std::optional<const double*>& k_ptr = std::nullopt, const std::optional<const double*>& Pr_ptr = std::nullopt);
    void get_thermo_reference_properties(const char* name, double* const M_ptr, double* const T_ref_ptr, double* const h0_ref_ptr);

    void write_debug_output(const char* out_filename);
  };


  std::vector<Problem> read_legacy_input(const char*);
  void run_all_cases(std::vector<Problem>&, const std::optional<const char*>& out_filename = std::nullopt, const std::optional<const char*>& plt_filename = std::nullopt);

  bool ffi_need_aligned(const std::vector<Problem>&);
}
