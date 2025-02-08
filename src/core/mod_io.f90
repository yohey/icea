module mod_io
  implicit none

contains

  subroutine write_debug_output(cea, filename)
    use mod_constants
    use mod_cea_types, only: CEA_Problem, CEA_Point

    type(CEA_Problem), intent(in):: cea
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
