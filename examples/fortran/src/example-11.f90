subroutine run_example_11()
  use cea
  implicit none

  type(CEA_Problem):: prob

  call prob%set_problem(mode = 'rocket', name = 'Example 11', equilibrium = .true., mole_ratios = .true., ions = .true.)
  call prob%set_output_options(SI = .true., transport = .true.)

  call prob%add_reactant('fuel', 'Li(cr)', ratio = 1.0000d0, T = 298.15d0)
  call prob%add_reactant('oxyd', 'F2(L)',  ratio = 0.5556d0, T =  85.02d0)

  call prob%set_chamber_pressures([1000d0], unit = 'legacy-psi')
  call prob%set_pressure_ratios([68.0457d0])
  call prob%set_subsonic_area_ratios([10d0])
  call prob%set_supersonic_area_ratios([10d0, 20d0, 100d0])

  call prob%set_legacy_mode(.true.)

  call prob%run(out_filename = 'example-11.out')

  return
end subroutine run_example_11
