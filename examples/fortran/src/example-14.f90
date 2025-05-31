subroutine run_example_14()
  use cea
  implicit none

  type(CEA_Problem):: prob

  call prob%set_problem(mode = 'tp', name = 'Example 14', mole_ratios = .true.)
  call prob%set_output_options(SI = .true., debug_points = [5])

  call prob%add_reactant('name', 'H2(L)', 'H 2', ratio = 100d0, rho =  874d0)
  call prob%add_reactant('name', 'O2(L)', 'O 2', ratio =  60d0, rho = 1431d0)

  call prob%set_chamber_pressures([0.05d0], unit = 'atm')
  call prob%set_chamber_temperatures([1000d0, 500d0, 350d0, 305d0, 304.3d0, 304.2d0, 304d0, 300d0])

  call prob%set_legacy_mode(.true.)

  !! copy example-13 results to emulate known bug
  prob%Mu(10) = -115.20613749687224d0
  prob%Mu(11) = -172.80920624530833d0
  prob%Mu(12) = -230.41227499374443d0

  call prob%run(out_filename = 'example-14.out')

  return
end subroutine run_example_14
