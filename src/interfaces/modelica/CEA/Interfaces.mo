within CEA;
package Interfaces
  extends Modelica.Icons.Package;
  import Modelica.Units.SI;

  class CEA_Problem
    extends ExternalObject;

    function constructor
      output CEA_Problem prob;
    external "C" prob = _ffi_cea_new_problem_C_impl()
        annotation(Include = "#include \"ffi_cea_modelica.h\"",
                   Library = {"cea_modelica", "-lc++", "-lc++abi", "-lgfortran"},
                   IncludeDirectory = "modelica://CEA/Resources/Include",
                   LibraryDirectory = "modelica://CEA/Resources/Library");
    end constructor;

    function destructor
      input CEA_Problem prob;
    external "C" _ffi_cea_del_problem_C_impl(prob);
    end destructor;

  end CEA_Problem;

  function ffi_cea_set_problem
    input CEA_Problem prob;
    input String mode;
    input String name = "Modelica Interface";
    input Boolean mole_ratios = false;
    input Boolean equilibrium = false;
    input Boolean ions = false;
    input Boolean frozen = false;
    input Boolean frozen_at_throat = false;
    input Boolean incident = false;
    input Boolean reflected = false;
  external "C" _ffi_cea_set_problem_C_impl(prob, mode, name, mole_ratios, equilibrium, ions, frozen, frozen_at_throat, incident, reflected);
  end ffi_cea_set_problem;

  function ffi_cea_set_output_options
    input CEA_Problem prob;
    input Boolean SI = false;
    input Integer debug_points[:] = fill(0, 0);
    input Boolean mass_fractions = false;
    input Boolean short = false;
    input Real trace_tol = 0;
    input Boolean transport = false;
    input String plot[:] = fill("", 0);
  external "C" _ffi_cea_set_output_options_C_impl(prob, SI, debug_points, size(debug_points, 1), mass_fractions, short, trace_tol, transport, plot, size(plot, 1));
  end ffi_cea_set_output_options;

  function ffi_cea_set_chamber_pressures
    input CEA_Problem prob;
    input SI.Pressure pressures[:];
  external "C" _ffi_cea_set_chamber_pressures_C_impl(prob, pressures, size(pressures, 1), "Pa");
  end ffi_cea_set_chamber_pressures;

  function ffi_cea_set_mixture_ratios
    input CEA_Problem prob;
    input SI.DimensionlessRatio ratios[:];
    input String _type = "";
  external "C" _ffi_cea_set_mixture_ratios_C_impl(prob, ratios, size(ratios, 1), _type);
  end ffi_cea_set_mixture_ratios;

  function ffi_cea_set_pressure_ratios
    input CEA_Problem prob;
    input SI.DimensionlessRatio ratios[:];
  external "C" _ffi_cea_set_pressure_ratios_C_impl(prob, ratios, size(ratios, 1));
  end ffi_cea_set_pressure_ratios;

  function ffi_cea_set_subsonic_area_ratios
    input CEA_Problem prob;
    input SI.DimensionlessRatio ratios[:];
  external "C" _ffi_cea_set_subsonic_area_ratios_C_impl(prob, ratios, size(ratios, 1));
  end ffi_cea_set_subsonic_area_ratios;

  function ffi_cea_set_supersonic_area_ratios
    input CEA_Problem prob;
    input SI.DimensionlessRatio ratios[:];
  external "C" _ffi_cea_set_supersonic_area_ratios_C_impl(prob, ratios, size(ratios, 1));
  end ffi_cea_set_supersonic_area_ratios;

  function ffi_cea_add_reactant
    input CEA_Problem prob;
    input String _type;
    input String name;
    input String formula = "";
    input Real ratio = ffi_cea_quiet_NaN();
    input SI.Temperature T = ffi_cea_quiet_NaN();
    input SI.Density rho = ffi_cea_quiet_NaN();
    input SI.SpecificEnthalpy h = ffi_cea_quiet_NaN();
    input SI.SpecificInternalEnergy u = ffi_cea_quiet_NaN();
    input String T_unit = "";
    input String rho_unit = "";
    input String h_unit = "";
    input String u_unit = "";
  external "C" _ffi_cea_add_reactant_C_impl(prob, _type, name, formula, ratio, T, rho, h, u, T_unit, rho_unit, h_unit, u_unit);
  end ffi_cea_add_reactant;

  function ffi_cea_set_reactant
    input CEA_Problem prob;
    input Integer index;
    input Real ratio = ffi_cea_quiet_NaN();
    input SI.Temperature T = ffi_cea_quiet_NaN();
    input SI.Density rho = ffi_cea_quiet_NaN();
    input SI.SpecificEnthalpy h = ffi_cea_quiet_NaN();
    input SI.SpecificInternalEnergy u = ffi_cea_quiet_NaN();
    input String T_unit = "";
    input String rho_unit = "";
    input String h_unit = "";
    input String u_unit = "";
  external "C" _ffi_cea_set_reactant_C_impl(prob, index, ratio, T, rho, h, u, T_unit, rho_unit, h_unit, u_unit);
  end ffi_cea_set_reactant;

  function ffi_cea_set_omit_species
    input CEA_Problem prob;
    input String species[:] = fill("", 0);
  external "C" _ffi_cea_set_omit_species_C_impl(prob, species, size(species, 1));
  end ffi_cea_set_omit_species;

  function ffi_cea_set_legacy_mode
    input CEA_Problem prob;
    input Boolean legacy_mode = true;
  external "C" _ffi_cea_set_legacy_mode_C_impl(prob, legacy_mode);
  end ffi_cea_set_legacy_mode;

  function ffi_cea_run
    input CEA_Problem prob;
    input String out_filename = "";
    input String plt_filename = "";
  external "C" _ffi_cea_run_C_impl(prob, out_filename, plt_filename);
  end ffi_cea_run;

  function ffi_cea_get_chamber_temperature
    input CEA_Problem prob;
    output SI.Temperature T;
  external "C" T = _ffi_cea_get_chamber_temperature_C_impl(prob);
  end ffi_cea_get_chamber_temperature;

  function ffi_cea_write_debug_output
    input CEA_Problem prob;
    input String filename;
  external "C" _ffi_cea_write_debug_output_C_impl(prob, filename);
  end ffi_cea_write_debug_output;

  function ffi_cea_get_thermo_reference_properties
    input CEA_Problem prob;
    input String name;
    output SI.MolarMass M;
    output SI.Temperature T_ref;
    output SI.SpecificEnthalpy h0_ref;
  external "C" _ffi_cea_get_thermo_reference_properties_C_impl(prob, name, M, T_ref, h0_ref);
  end ffi_cea_get_thermo_reference_properties;

  function ffi_cea_sizeof
    input CEA_Problem prob;
    output Integer n;
  external "C" n = _ffi_cea_sizeof_C_impl(prob);
  end ffi_cea_sizeof;

  function ffi_cea_quiet_NaN
    output Real NaN;
  external "C" NaN = _ffi_cea_quiet_NaN_C_impl();
  end ffi_cea_quiet_NaN;

end Interfaces;
