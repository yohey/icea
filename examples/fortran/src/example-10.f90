subroutine run_example_10()
  use cea
  implicit none

  type(CEA_Problem):: prob

  call prob%set_problem(mode = 'rocket', name = 'Example 10', equilibrium = .false.)
  call prob%set_output_options(SI = .true., short = .true.)

  call prob%add_reactant('fuel', 'H2(L)', 100d0, T = 20.27d0)
  call prob%add_reactant('oxyd', 'O2(L)', 100d0, T = 90.17d0)

  call prob%set_finite_area_combustor(mass_flow_ratio = 1333.9d0)
  call prob%set_chamber_pressures([53.3172d0])
  call prob%set_mixture_ratios([5.55157d0])
  call prob%set_pressure_ratios([1d1, 1d2, 1d3])
  call prob%set_supersonic_area_ratios([25d0, 50d0, 75d0])

  call prob%set_legacy_mode(.true.)

  call prob%run(out_filename = 'example-10.out')

  return
end subroutine run_example_10
