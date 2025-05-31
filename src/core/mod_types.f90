module mod_types
  use mod_constants
  implicit none

  integer, parameter:: maxNgc = 600
  integer, parameter:: maxNc  = 300
  integer, parameter:: Ncol   = 8
  integer, parameter:: maxMat = 50
  integer, parameter:: maxTr  = 40
  integer, parameter:: maxR   = 24
  integer, parameter:: maxEl  = 20
  integer, parameter:: maxNg  = 500
  integer, parameter:: maxMix = 52
  integer, parameter:: maxT   = 51
  integer, parameter:: maxPv  = 26

  ! The following parameters set the input/output unit numbers.
  ! These numbers are also defined in the manual, part 2 p. 39,
  ! and may be adjusted as desired.
  integer:: IOINP = 0
  integer:: IOOUT = 0
  !***********************************************************************

  type:: ThermoProperty
     integer:: type
     character(15):: name
     integer:: ntl
     character(6):: date
     character(2):: sym(5)
     real(8):: fno(5)
     integer:: ifaz
     real(8):: tl(2)
     real(8):: mwt
     real(8):: thermo(9, 3)
  end type ThermoProperty


  type:: TransportProperty
     integer:: j(2)
     real(8):: data(6, 3, 2)
  end type TransportProperty


  type:: CEA_Point
     logical:: Debug
     real(8):: En(maxNgc)
     real(8):: Cpr, Dlvpt, Dlvtp, Gammas, Hsum
     real(8):: Ppp, Ssum, Totn, Ttt, Vlm, Wm
     real(8):: Sonvel, Spim, Vmoc
     ! ROCKET
     real(8):: AeAt, App, Ivac
     ! SHCK
     real(8):: U1, Mach1
     ! TRANP
     real(8):: Coneql, Confro, Cpeql, Cpfro, Preql, Prfro, Vis
  end type CEA_Point


  type:: CEA_Core_Problem
     integer:: io_log = 0

     logical:: invalid_case = .false.
     logical:: legacy_mode = .false.

     character(MAX_FILENAME):: filename_thermo_lib
     character(MAX_FILENAME):: filename_trans_lib

     type(ThermoProperty), pointer:: thermo_properties(:) => null()
     type(TransportProperty), allocatable:: transport_properties(:)
     type(CEA_Point), pointer:: points(:, :) => null()

     integer:: iOF !< index of O/F
     integer:: ipt !< index of point (new)
     integer:: max_points

     integer:: Nonly_in, Nreac_in
     logical:: Eql_in
     character(15):: Prod_in(0:maxNgc)
     integer:: Npp_in, Nsub_in, Nsup_in
     real(8):: AcAt_in, mdotByAc_in
     integer, allocatable:: Debug_in(:)
     real(8), allocatable:: U1_in(:), Ma1_in(:)
     logical:: eqrats_in
     logical:: incd_in = .false.
     logical:: refl_in = .false.
     real(8), allocatable:: Oxf_in(:)

     character(15):: ensert(20)

     ! for subroutine SETEN
     integer:: Isv, lsave, lsav
     real(8):: Enlsav, Ensave, Tsave

     real(8):: Enn, Ennl, Sumn
     real(8):: Deln(maxNgc), Enln(maxNgc), Sln(maxNgc)

     integer:: Iplt, Nc, Ng, Ngp1, Nlm, Nplt, Nof, Nomit, Nonly, Np, Npr, Ngc, Nsert, Nspr, Nspx, Nt
     integer:: Jcond(45), Jx(maxEl), Nfla(maxR), Ifz(maxNc)

     real(8):: Cpmix, Wmix, Bcheck
     real(8):: Am(2), Hpp(2), Vmin(2), Vpls(2), Wp(2), Oxf(maxMix), P(maxPv), Rh(2), T(maxT), V(maxPv)
     real(8):: B0p(maxEl, 2)

     integer:: Imat, Iq1, Jliq, Jsol, Msing

     logical:: Convg, Detdbg, Detn, Eql, Gonly, Hp, Ions, Massf, Moles
     logical:: Pderiv, Shock, Short, SIunit, Sp, Tp, Trnspt, Vol

     real(8):: Eqrat, Hsub0, Oxfl, Pp, R, Size, S0, Tln, Tm, Trace, Tt, Viscns, Vv
     real(8):: Atwt(maxEl), B0(maxEl), X(maxMat)
     real(8):: A(maxEl, maxNgc), G(maxMat, maxMat+1)

     character(2):: Elmt(maxEl), Ratom(maxR, 12)
     character(8):: Fox(maxR)
     character(10):: Thdate
     character(15):: Case, Energy(maxR), Omit(0:maxNgc), Pltvar(20), Prod(0:maxNgc), Rname(maxR)

     real(8):: Pltout(500, 20)

     integer:: Nreac
     integer:: Jray(maxR)
     real(8):: Dens(maxR), Enth(maxR), Pecwt(maxR), Rmw(maxR), Rtemp(maxR), Rnum(maxR, 12)

     real(8):: Cpsum
     real(8):: Cft(9, maxNc), Coef(9, maxNg, 3), Temp(2, maxNc)
     real(8):: Cp(maxNgc), H0(maxNgc), Mu(maxNgc), Mw(maxNgc), S(maxNgc), Tg(4)

     integer:: Iopt, Isup, Nfz, Npp, Nsub, Nsup
     logical:: Area, Debugf, Fac, Froz, Page1, Rkt
     real(8):: AcAt, Awt, Cstr, Tcest, mdotByAc
     real(8):: Pcp(2*Ncol), Subar(13), Supar(13)

     integer:: Nsk
     logical:: Incdeq, Incdfz, Refleq, Reflfz, Shkdbg
     real(8):: A1, Gamma1

     integer:: Nm, Nr, Ntape
     integer:: Ind(maxTr), Jcm(maxEl)
     real(8):: Cprr(maxTr), Con(maxTr), Wmol(maxTr), Xs(maxTr)
     real(8):: Eta(maxTr, maxTr), Stc(maxTr, maxTr)

     ! for subroutine EQLBRM
     logical:: newcom
     real(8):: ensol, pisave(maxMat-2)

     ! for subroutine TRANIN
     integer:: Nms, inds(maxTr)
     real(8):: xss(maxTr), wmols(maxTr)

     ! for subroutine OUT3
     integer:: num_omitted

     ! Information used in variable output format
     character(4):: fmt(30) = [character(4):: '(1X', ',A15', ',', 'F9.', '0,', 'F9.', '0,', 'F9.', '0,', &
                               'F9.', '0,', 'F9.', '0,', 'F9.', '0,', 'F9.', '0,', 'F9.', '0,', &
                               'F9.', '0,', 'F9.', '0,', 'F9.', '0,', 'F9.', '0,', 'F9.', '0', ')']

   contains
     final:: del_problem

  end type CEA_Core_Problem

contains

  subroutine init_case(cea)
    class(CEA_Core_Problem), intent(inout):: cea

    cea%Case = 'New Case'

    cea%Detn = .false.
    cea%Rkt = .false.
    cea%Shock = .false.
    cea%Hp = .false.
    cea%Sp = .false.
    cea%Tp = .false.
    cea%Vol = .false.

    cea%Moles = .false.
    cea%Eql_in = .false.
    cea%Eql = .false.
    cea%Ions = .false.
    cea%Froz = .false.

    cea%Fac = .false.

    cea%SIunit = .true.
    cea%Trnspt = .false.
    cea%Debugf = .false.
    cea%Shkdbg = .false.
    cea%Incdeq = .false.
    cea%Incdfz = .false.
    cea%Refleq = .false.
    cea%Reflfz = .false.

    cea%Nt = 1
    cea%Nreac_in = 0
    cea%Nonly_in = 0
    cea%Nsub_in = 0
    cea%Nsup_in = 0
    cea%AcAt_in = 0
    cea%mdotByAc_in = 0

    cea%Hsub0 = 1d30

    cea%Jray(:) = 0

    cea%ensert(:) = ""

    cea%newcom = .false.

    return
  end subroutine init_case


  subroutine del_problem(cea)
    type(CEA_Core_Problem), intent(inout):: cea

#ifndef NDEBUG
    write(0, *) '[DEBUG] CEA_Core_Problem (mod_types.f90) destructor is called: ', trim(cea%Case)
#endif

    if (associated(cea%points)) deallocate(cea%points)
    if (associated(cea%thermo_properties)) deallocate(cea%thermo_properties)
    if (allocated(cea%transport_properties)) deallocate(cea%transport_properties)
    if (allocated(cea%Debug_in)) deallocate(cea%Debug_in)
    if (allocated(cea%U1_in)) deallocate(cea%U1_in)
    if (allocated(cea%Ma1_in)) deallocate(cea%Ma1_in)

    return
  end subroutine del_problem


  subroutine allocate_points(cea)
    class(CEA_Core_Problem), intent(inout):: cea
    integer:: i, iof, ipt, Nof_tmp

    if (cea%Shock) then
       cea%max_points = cea%Nsk
    else
       cea%max_points = cea%Nt * cea%Np
       if (cea%Rkt) then
          cea%max_points = cea%max_points * (cea%Npp_in + cea%Nsub_in + cea%Nsup_in)
          if (cea%Fac) then
             cea%max_points = cea%max_points + 3
          else
             cea%max_points = cea%max_points + 2
          end if
       end if
    end if

    if (cea%max_points > 0) then
       Nof_tmp = max(1, cea%Nof)

       allocate(cea%points(Nof_tmp, cea%max_points))

       do concurrent (iof = 1:Nof_tmp, ipt = 1:cea%max_points)
          cea%points(iof, ipt)%Debug = .false.
          cea%points(iof, ipt)%En(:) = 0
       end do

       if (allocated(cea%U1_in)) then
          do concurrent (iof = 1:Nof_tmp, ipt = 1:cea%Nsk)
             cea%points(iof, ipt)%U1 = cea%U1_in(ipt)
             cea%points(iof, ipt)%Mach1 = 0
          end do
       else if (allocated(cea%Ma1_in)) then
          do concurrent (iof = 1:Nof_tmp, ipt = 1:cea%Nsk)
             cea%points(iof, ipt)%U1 = 0
             cea%points(iof, ipt)%Mach1 = cea%Ma1_in(ipt)
          end do
       end if

       if (allocated(cea%Debug_in)) then
          do concurrent (iof = 1:Nof_tmp, i = 1:size(cea%Debug_in))
             cea%points(iof, cea%Debug_in(i))%Debug = .true.
          end do
       end if
    end if

    return
  end subroutine allocate_points


  subroutine reset_case(cea)
    class(CEA_Core_Problem), intent(inout):: cea
    integer:: i
    real(8):: xi, xyz, denmtr

    cea%Eql = cea%Eql_in
    cea%Npp = cea%Npp_in
    cea%Nsub = cea%Nsub_in
    cea%Nsup = cea%Nsup_in
    cea%AcAt = cea%AcAt_in
    cea%mdotByAc = cea%mdotByAc_in


    !! subroutine INPUT
    !! if (code(1:3) == 'end')
    if (cea%Nof == 0) then
       cea%Nof = 1
       if (cea%Wp(2) > 0) then
          cea%Oxf(1) = cea%Wp(1) / cea%Wp(2)
       else
          cea%Oxf(1) = 0
       end if
    else if (cea%eqrats_in) then
       do concurrent (i = 1:cea%Nof, cea%Oxf_in(i) > 0)
          xyz = -cea%Oxf_in(i) * cea%Vmin(2) - cea%Vpls(2)
          denmtr = cea%Oxf_in(i) * cea%Vmin(1) + cea%Vpls(1)
          if (abs(denmtr) < 1d-30) then
             cea%Oxf(i) = 0
          else
             cea%Oxf(i) = xyz / denmtr
          end if
       end do
    else
       cea%Oxf(1:cea%Nof) = cea%Oxf_in(1:cea%Nof)
    end if

    if (cea%Trnspt) cea%Viscns = 0.3125d0 * sqrt(1d5 * Boltz / (pi * Avgdr))

    if (cea%SIunit) then
       cea%R = R0 / 1d3
    else
       cea%R = R0 / (cal_to_J * 1d3)
    end if

    cea%iOF = 0
    cea%num_omitted = 0

    !! migrate from INIT
    if (cea%Shock) then
       cea%Incdfz = cea%incd_in .and. cea%Froz
       cea%Incdeq = cea%incd_in .and. cea%Eql_in
       cea%Reflfz = cea%refl_in .and. cea%Froz
       cea%Refleq = cea%refl_in .and. cea%Eql_in
    end if

    !! migrate from read_libraries
    if (associated(cea%points)) then
       if (cea%Ng == 0) then
          xi = cea%Enn
       else
          xi = cea%Enn / cea%Ng
       end if

       do concurrent (i = 1:cea%Nof)
          cea%points(i, 1)%En(cea%Ng+1:cea%Ng+cea%Nc) = 0
          cea%points(i, 1)%En(1:cea%Ng) = xi
       end do

       cea%Enln(cea%Ng+1:cea%Ng+cea%Nc) = 0
       cea%Enln(1:cea%Ng) = log(xi)
    end if

    return
  end subroutine reset_case

end module mod_types
