subroutine run_example_06()
  use cea
  implicit none

  type(CEA_Problem):: prob

  call prob%set_problem(mode = 'deton', name = 'Example 6')
  call prob%set_output_options(SI = .false., transport = .true.)

  call prob%add_reactant('oxyd', 'O2', ratio = 100d0, T = 298.15d0)
  call prob%add_reactant('fuel', 'H2', ratio = 100d0, T = 298.15d0)

  call prob%set_chamber_temperatures([298.15d0, 500d0])
  call prob%set_mixture_ratios([1.0d0], type = 'eq.ratio')
  call prob%set_chamber_pressures([1d0, 20d0], unit = 'bar')

  call prob%set_legacy_mode(.true.)

  call prob%run(out_filename = 'example-06.out')

  return
end subroutine run_example_06
