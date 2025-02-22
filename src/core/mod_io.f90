module mod_io
  implicit none

  interface set_library_paths
     procedure:: set_library_paths_scalar
     procedure:: set_library_paths_array
  end interface set_library_paths

contains

  subroutine set_library_paths_scalar(cea, thermo_lib_candidate, trans_lib_candidate)
    use mod_constants, only: MAX_FILENAME
    use mod_types, only: CEA_Core_Problem

    class(CEA_Core_Problem), intent(inout):: cea
    character(*), intent(in), optional:: thermo_lib_candidate
    character(*), intent(in), optional:: trans_lib_candidate

    character(MAX_FILENAME):: thermo_lib_found
    character(MAX_FILENAME):: trans_lib_found

    call find_library_paths(thermo_lib_found, trans_lib_found, thermo_lib_candidate, trans_lib_candidate)

    cea%filename_thermo_lib = trim(thermo_lib_found)
    cea%filename_trans_lib = trim(trans_lib_found)

    return
  end subroutine set_library_paths_scalar

  subroutine set_library_paths_array(cea, thermo_lib_candidate, trans_lib_candidate)
    use mod_constants, only: MAX_FILENAME
    use mod_types, only: CEA_Core_Problem

    class(CEA_Core_Problem), dimension(:), intent(inout):: cea
    character(*), intent(in), optional:: thermo_lib_candidate
    character(*), intent(in), optional:: trans_lib_candidate

    character(MAX_FILENAME):: thermo_lib_found
    character(MAX_FILENAME):: trans_lib_found
    integer:: icase

    call find_library_paths(thermo_lib_found, trans_lib_found, thermo_lib_candidate, trans_lib_candidate)

    do icase = 1, size(cea)
       cea(icase)%filename_thermo_lib = trim(thermo_lib_found)
       cea(icase)%filename_trans_lib = trim(trans_lib_found)
    end do

    return
  end subroutine set_library_paths_array


  subroutine find_library_paths(thermo_lib_found, trans_lib_found, thermo_lib_candidate, trans_lib_candidate)
    use mod_constants, only: MAX_FILENAME

    character(MAX_FILENAME), intent(out):: thermo_lib_found
    character(MAX_FILENAME), intent(out):: trans_lib_found
    character(*), intent(in), optional:: thermo_lib_candidate
    character(*), intent(in), optional:: trans_lib_candidate
    character(MAX_FILENAME):: ENV_CEA_ROOT, tmp_path
    logical:: thermo_lib_exists, trans_lib_exists

    thermo_lib_exists = .false.
    trans_lib_exists = .false.

    if (present(thermo_lib_candidate)) then
       inquire(file = thermo_lib_candidate, exist = thermo_lib_exists)
       if (thermo_lib_exists) thermo_lib_found = trim(thermo_lib_candidate)
    end if

    if (present(trans_lib_candidate)) then
       inquire(file = trans_lib_candidate, exist = trans_lib_exists)
       if (trans_lib_exists) trans_lib_found = trim(trans_lib_candidate)
    end if

    if (thermo_lib_exists .and. trans_lib_exists) return

    if (.not. thermo_lib_exists) then
       tmp_path = './thermo.lib'
       inquire(file = tmp_path, exist = thermo_lib_exists)
       if (thermo_lib_exists) thermo_lib_found = tmp_path
    end if

    if (.not. trans_lib_exists) then
       tmp_path = './trans.lib'
       inquire(file = tmp_path, exist = trans_lib_exists)
       if (trans_lib_exists) trans_lib_found = tmp_path
    end if

    if (thermo_lib_exists .and. trans_lib_exists) return

    call get_environment_variable('CEA_ROOT', ENV_CEA_ROOT)

    if (len_trim(ENV_CEA_ROOT) > 0) then
       if (.not. thermo_lib_exists) then
          tmp_path = trim(ENV_CEA_ROOT) // '/lib/cea/thermo.lib'
          inquire(file = tmp_path, exist = thermo_lib_exists)
          if (thermo_lib_exists) thermo_lib_found = tmp_path
       end if

       if (.not. trans_lib_exists) then
          tmp_path = trim(ENV_CEA_ROOT) // '/lib/cea/trans.lib'
          inquire(file = tmp_path, exist = trans_lib_exists)
          if (trans_lib_exists) trans_lib_found = tmp_path
       end if
    end if

    if (thermo_lib_exists .and. trans_lib_exists) return

#ifdef CEA_ROOT
    if (.not. thermo_lib_exists) then
       tmp_path = trim(CEA_ROOT) // '/lib/cea/thermo.lib'
       inquire(file = tmp_path, exist = thermo_lib_exists)
       if (thermo_lib_exists) thermo_lib_found = tmp_path
    end if

    if (.not. trans_lib_exists) then
       tmp_path = trim(CEA_ROOT) // '/lib/cea/trans.lib'
       inquire(file = tmp_path, exist = trans_lib_exists)
       if (trans_lib_exists) trans_lib_found = tmp_path
    end if
#endif

    if (.not. thermo_lib_exists) then
       thermo_lib_found = ''
       write(0, *) '[WARNING] thermo.lib is not found.'
    end if

    if (.not. trans_lib_exists) then
       trans_lib_found = ''
       write(0, *) '[WARNING] trans.lib is not found.'
    end if

    return
  end subroutine find_library_paths


  subroutine write_debug_output(cea, filename)
    use mod_constants
    use mod_types, only: CEA_Core_Problem, CEA_Point

    class(CEA_Core_Problem), intent(in):: cea
    character(*), intent(in):: filename

    integer:: io_out, iof, ipt, j, max_points
    real(8), allocatable:: X(:, :)
    character(12), allocatable:: labels(:)
    type(CEA_Point), pointer:: p


    if (cea%Shock) then
       max_points = cea%Nsk
    else
       max_points = cea%Nt * cea%Np
       if (cea%Rkt) then
          max_points = max_points * (cea%Npp + cea%Nsub + cea%Nsup)
       end if
    end if

    write(0, *) '[DEBUG] max_points = ', max_points

    open(newunit = io_out, file = trim(filename), status = 'unknown', form = 'formatted')

    write(io_out, '("Case = ", a)') trim(cea%Case)

    write(io_out, '(/"OPTIONS:")')
    write(io_out, '("Nof  = ", i5)') cea%Nof
    write(io_out, '("Nt   = ", i5)') cea%Nt
    write(io_out, '("Np   = ", i5)') cea%Np
    write(io_out, '("Npp  = ", i5)') cea%Npp
    write(io_out, '("Nsub = ", i5)') cea%Nsub
    write(io_out, '("Nsup = ", i5)') cea%Nsup

    write(io_out, '("TP = ", l1)') cea%Tp
    write(io_out, '("HP = ", l1)') (cea%Hp .and. .not. cea%Vol)
    write(io_out, '("SP = ", l1)') cea%Sp
    write(io_out, '("TV = ", l1)') (cea%Tp .and. cea%Vol)
    write(io_out, '("UV = ", l1)') (cea%Hp .and. cea%Vol)
    write(io_out, '("SV = ", l1)') (cea%Sp .and. cea%Vol)

    write(io_out, '("DETN   = ", l1)') cea%Detn
    write(io_out, '("SHOCK  = ", l1)') cea%Shock
    write(io_out, '("REFL   = ", l1)') (cea%Shock .and. (cea%Reflfz .or. cea%Refleq))
    write(io_out, '("INCD   = ", l1)') (cea%Shock .and. (cea%Incdfz .or. cea%Incdeq))
    write(io_out, '("Rkt    = ", l1)') cea%Rkt
    write(io_out, '("FROZ   = ", l1)') cea%Froz
    write(io_out, '("EQL    = ", l1)') cea%Eql
    write(io_out, '("IONS   = ", l1)') cea%Ions
    write(io_out, '("SIUNIT = ", l1)') cea%SIunit
    write(io_out, '("DEBUGF = ", l1)') cea%Debugf
    write(io_out, '("SHKDBG = ", l1)') cea%Shkdbg
    write(io_out, '("DETDBG = ", l1)') cea%Detdbg
    write(io_out, '("TRNSPT = ", l1)') cea%Trnspt

    if (max_points == 0) return

    allocate(X(20, cea%Nof * max_points))

    if (cea%Rkt) then
       allocate(labels(cea%Nof * max_points))
       labels(:) = 'EXIT'
       if (cea%Fac) then
          labels(1) = 'INJECTOR'
          labels(2) = 'COMB END'
          labels(3) = 'THROAT'
       else
          labels(1) = 'CHAMBER'
          labels(2) = 'THROAT'
       end if
    end if

    do concurrent (iof = 1:cea%Nof, ipt = 1:max_points)
       p => cea%points(iof, ipt)
       j = (iof - 1) * max_points + ipt
       X(1, j) = cea%Oxf(iof)
       X(2, j) = p%Ppp * 0.1d0
       X(3, j) = p%Ttt
       if (p%Vlm == 0) then
          X(4, j) = 0
       else
          X(4, j) = 1d5 / p%Vlm
       end if
       X(5, j) = p%Hsum * cea%R
       X(6, j) = (p%Hsum - p%Ppp * p%Vlm / R0) * cea%R
       X(7, j) = (p%Hsum - p%Ttt * p%Ssum) * cea%R
       X(8, j) = p%Ssum * cea%R
       X(9, j) = p%Wm
       X(10, j) = p%Dlvpt
       X(11, j) = p%Dlvtp
       X(12, j) = p%Cpr * cea%R
       X(13, j) = p%Gammas
       X(14, j) = p%Sonvel
       X(15, j) = p%Vmoc
    end do

    do concurrent (iof = 1:cea%Nof, ipt = 2:max_points)
       p => cea%points(iof, ipt)
       j = (iof - 1) * max_points + ipt
       X(16, j) = p%AeAt
       X(17, j) = cea%Cstr
       if (cea%Cstr == 0) then
          X(18, j) = 0
       else
          X(18, j) = p%Spim / cea%Cstr
       end if
       if (p%Wm * p%Spim == 0) then
          X(19, j) = 0
       else
          X(19, j) = (p%Spim + R0 * p%Ttt / (p%Wm * p%Spim)) / g0
       end if
       X(20, j) = p%Spim / g0
    end do

    if (.not. cea%SIunit) then
       X(5:8, :) = X(5:8, :) * cal_to_J
       X(12, :) = X(12, :) * cal_to_J
    end if

    write(io_out, *)
    if (cea%Rkt) write(io_out, '(17x, *(a12))') adjustr(labels(:))
    write(io_out, '("O/F            = ", *(f12.5))') X(1, :)
    write(io_out, '("P [MPa]        = ", *(f12.7))') X(2, :)
    write(io_out, '("T [K]          = ", *(f12.2))') X(3, :)
    write(io_out, '("ρ [kg/m^3]     = ", *(f12.6))') X(4, :)
    write(io_out, '("H [kJ/kg]      = ", *(f12.2))') X(5, :)
    write(io_out, '("U [kJ/kg]      = ", *(f12.2))') X(6, :)
    write(io_out, '("G [kJ/kg]      = ", *(f12.1))') X(7, :)
    write(io_out, '("S [kJ/(kg·K)]  = ", *(f12.4))') X(8, :)

    write(io_out, '("M [g/mol]      = ", *(f12.3))') X(9, :)
    write(io_out, '("(dlnV/dlnP)T   = ", *(f12.4))') X(10, :)
    write(io_out, '("(dlnV/dlnT)P   = ", *(f12.4))') X(11, :)
    write(io_out, '("cp [kJ/(kg·K)] = ", *(f12.4))') X(12, :)
    write(io_out, '("γ [-]          = ", *(f12.4))') X(13, :)
    write(io_out, '("a [m/s]        = ", *(f12.2))') X(14, :)
    write(io_out, '("Ma [-]         = ", *(f12.4))') X(15, :)

    write(io_out, '("Ae/At [-]      = ", 12x, *(f12.4))') X(16, 2:)
    write(io_out, '("C* [m/s]       = ", 12x, *(f12.2))') X(17, 2:)
    write(io_out, '("CF [-]         = ", 12x, *(f12.4))') X(18, 2:)
    write(io_out, '("Ivac [s]       = ", 12x, *(f12.3))') X(19, 2:)
    write(io_out, '("Isp [s]        = ", 12x, *(f12.3))') X(20, 2:)

    deallocate(X)
    if (cea%Rkt) deallocate(labels)

    close(io_out)

    return
  end subroutine write_debug_output

end module mod_io
