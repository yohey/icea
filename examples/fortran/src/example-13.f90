subroutine run_example_13()
  use cea
  implicit none

  type(CEA_Problem):: prob

  call prob%set_problem(mode = 'rocket', name = 'Example 13', equilibrium = .true.)
  call prob%set_output_options(SI = .false., trace_tol = 1d-10)

  call prob%add_reactant('fuel', 'N2H4(L)', ratio = 0.8d0, T = 298.15d0, rho =  874d0)
  call prob%add_reactant('fuel', 'Be(a)',   ratio = 0.2d0, T = 298.15d0, rho = 1431d0)
  call prob%add_reactant('oxyd', 'H2O2(L)', ratio = 1.0d0, T = 298.15d0)

  call prob%set_chamber_pressures([3000d0], unit = 'legacy-psi')
  call prob%set_mixture_ratios([67d0], type = '%fuel')
  call prob%set_pressure_ratios([3d0, 10d0, 30d0, 300d0])

  call prob%set_insert_species(['BeO(L)'])

  call prob%set_legacy_mode(.true.)

  call prob%run(out_filename = 'example-13.out')

  return
end subroutine run_example_13
