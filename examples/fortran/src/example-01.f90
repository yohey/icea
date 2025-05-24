subroutine run_example_01()
  use cea
  implicit none

  type(CEA_Problem):: prob

  call prob%set_problem(mode = 'tp', name = 'Example 1', mole_ratios = .true.)
  call prob%set_output_options(SI = .false.)

  call prob%add_reactant('fuel', 'H2',  ratio = 1.0d0)
  call prob%add_reactant('oxyd', 'Air', ratio = 1.0d0)

  call prob%set_chamber_pressures([1.0d0, 0.1d0, 0.01d0], unit = 'atm')
  call prob%set_chamber_temperatures([3000d0, 2000d0])
  call prob%set_mixture_ratios([1.0d0, 1.5d0], type = 'eq.ratio')

  call prob%set_only_species(['Ar  ', 'C   ', 'CO  ', 'CO2 ', 'H   ', 'H2  ', 'H2O ', 'HNO ', 'HO2 ', 'HNO2', 'HNO3', &
                              'N   ', 'NH  ', 'NO  ', 'N2  ', 'N2O3', 'O   ', 'O2  ', 'OH  ', 'O3  '])

  call prob%set_legacy_mode(.true.)

  call prob%run(out_filename = 'example-01.out')

  return
end subroutine run_example_01
