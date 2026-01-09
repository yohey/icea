within CEA.Examples;
model Example08
  extends Modelica.Icons.Example;
  import CEA.Interfaces.*;

  CEA_Problem prob = CEA_Problem();

initial algorithm
  ffi_cea_set_problem(prob, mode = "rocket", name = "Example 8", equilibrium = true);

  ffi_cea_set_output_options(prob, SI = true);

  ffi_cea_add_reactant(prob, "fuel", "H2(L)", ratio = 100, T = 20.27);
  ffi_cea_add_reactant(prob, "oxyd", "O2(L)", ratio = 100, T = 90.17);

  ffi_cea_set_chamber_pressures(prob, {5.33172e6});
  ffi_cea_set_mixture_ratios(prob, {5.55157});
  ffi_cea_set_pressure_ratios(prob, {10, 100, 1000});
  ffi_cea_set_subsonic_area_ratios(prob, {1.58});
  ffi_cea_set_supersonic_area_ratios(prob, {25, 50, 75});

  ffi_cea_set_legacy_mode(prob, true);

  ffi_cea_run(prob, "example-08.out");
end Example08;
