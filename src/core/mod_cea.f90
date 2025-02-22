!-------------------------------------------------------------------------------
!> @mainpage CEA: Chemical Equlibrium with Applications
!!
!-------------------------------------------------------------------------------
module mod_cea
  implicit none

contains

  subroutine run_all_cases(cea, out_filename, plt_filename)
    use mod_constants, only: stderr
    use mod_types, only: CEA_Core_Problem, MAX_FILENAME, IOOUT, allocate_points, reset_case
    use mod_legacy_io
    implicit none

    class(CEA_Core_Problem), intent(inout):: cea(:)
    character(*), intent(in), optional:: out_filename
    character(*), intent(in), optional:: plt_filename
    integer:: icase, num_cases
    logical:: any_legacy_mode, is_opened

    num_cases = size(cea)

    any_legacy_mode = any([(cea(icase)%legacy_mode, icase = 1, num_cases)])

    do icase = 1, num_cases
       if (cea(icase)%legacy_mode) then
          inquire(cea(icase)%io_log, opened = is_opened)
          if (cea(icase)%io_log == 0 .or. .not. is_opened) open(newunit = cea(icase)%io_log, status = 'scratch', form = 'formatted')
       end if

       if (cea(icase)%Nlm == 0) call REACT(cea(icase))

       if (.not. associated(cea(icase)%points)) call allocate_points(cea(icase))

       call read_libraries(cea(icase))
    end do

    if (any_legacy_mode) then
       call open_legacy_output(IOOUT, out_filename)
    end if

    do icase = 1, num_cases
       call reset_case(cea(icase))

       if (cea(icase)%legacy_mode) then
          !! TEMPORARY WORK AROUND TO REPRODUCE KNOWN BUG !!
          if (icase >= 2) then
             cea(icase)%Mu(:) = cea(icase-1)%Mu(:)
          end if
          !!!!!!!!!!!!!!!!! TO BE DELETED !!!!!!!!!!!!!!!!!!
       end if

       if (cea(icase)%legacy_mode) call write_input_log(cea(icase)%io_log, IOOUT)

       if (cea(icase)%invalid_case) cycle

       if (cea(icase)%Rkt) then
          call ROCKET(cea(icase))
       else if (cea(icase)%Tp .or. cea(icase)%Hp .or. cea(icase)%Sp) then
          call THERMP(cea(icase))
       else if (cea(icase)%Detn) then
          call DETON(cea(icase))
       else if (cea(icase)%Shock) then
          call SHCK(cea(icase))
       end if
    end do

    inquire(IOOUT, opened = is_opened)
    if (is_opened) then
       close(IOOUT)
       IOOUT = 0
    end if

    if (present(plt_filename)) then
       call write_plt_file(cea(1:num_cases), plt_filename)
    end if

    return
  end subroutine run_all_cases


  subroutine run_case(cea, out_filename)
    use mod_types
    use mod_legacy_io, only: REACT, open_legacy_output, write_input_log, write_plt_file

    class(CEA_Core_Problem), intent(inout):: cea
    character(*), intent(in), optional:: out_filename
    logical:: is_opened

    if (cea%legacy_mode) then
       inquire(cea%io_log, opened = is_opened)
       if (cea%io_log == 0 .or. .not. is_opened) open(newunit = cea%io_log, status = 'scratch', form = 'formatted')
    end if

    if (cea%Nlm == 0) call REACT(cea)

    if (.not. associated(cea%points)) call allocate_points(cea)

    call read_libraries(cea)

    if (cea%legacy_mode) then
       call open_legacy_output(IOOUT, out_filename)
    end if

    call reset_case(cea)

    if (cea%legacy_mode) call write_input_log(cea%io_log, IOOUT)

    if (cea%invalid_case) return

    if (cea%Rkt) then
       call ROCKET(cea)
    else if (cea%Tp .or. cea%Hp .or. cea%Sp) then
       call THERMP(cea)
    else if (cea%Detn) then
       call DETON(cea)
    else if (cea%Shock) then
       call SHCK(cea)
    end if

    inquire(IOOUT, opened = is_opened)
    if (is_opened) then
       close(IOOUT)
       IOOUT = 0
    end if

    return
  end subroutine run_case


  subroutine CPHS(cea)
    !***********************************************************************
    ! CALCULATES THERMODYNAMIC PROPERTIES FOR INDIVIDUAL SPECIES
    !***********************************************************************
    use mod_types
    use mod_functions
    implicit none

    type(CEA_Core_Problem), intent(inout):: cea

    integer:: ij, j, jj, k

    k = 1
    if (cea%Tt > cea%Tg(2)) k = 2
    if (cea%Tt > cea%Tg(3)) k = 3

    do concurrent (j = 1:cea%Ng)
       cea%S(j) = calc_entropy(cea%Coef(:, j, k), cea%Tt, cea%Tln, .true.)
       cea%H0(j) = calc_enthalpy(cea%Coef(:, j, k), cea%Tt, cea%Tln, .true.)
    end do

    if (.not. cea%Tp .or. cea%Convg) then
       do concurrent (j = 1:cea%Ng)
          cea%Cp(j) = calc_specific_heat(cea%Coef(:, j, k), cea%Tt)
       end do
    end if

    if (cea%Npr /= 0 .and. k /= 3 .and. cea%Ng /= cea%Ngc) then
       do concurrent (ij = 1:cea%Npr)
          j = cea%Jcond(ij)
          jj = cea%Jcond(ij) - cea%Ng

          cea%S(j) = calc_entropy(cea%Cft(:, jj), cea%Tt, cea%Tln, .true.)
          cea%H0(j) = calc_enthalpy(cea%Cft(:, jj), cea%Tt, cea%Tln, .true.)
          cea%Cp(j) = calc_specific_heat(cea%Cft(:, jj), cea%Tt)
       end do
    end if

    return
  end subroutine CPHS


  subroutine ALLCON(cea)
    !***********************************************************************
    ! CALCULATES THERMODYNAMIC PROPERTIES FOR INDIVIDUAL SPECIES
    !***********************************************************************
    use mod_types
    use mod_functions
    implicit none

    type(CEA_Core_Problem), intent(inout):: cea

    integer:: j, jj

    do concurrent (jj = 1:cea%Nc)
       j = jj + cea%Ng
       cea%S(j) = calc_entropy(cea%Cft(:, jj), cea%Tt, cea%Tln, .true.)
       cea%H0(j) = calc_enthalpy(cea%Cft(:, jj), cea%Tt, cea%Tln, .true.)
       cea%Cp(j) = calc_specific_heat(cea%Cft(:, jj), cea%Tt)
    end do

    return
  end subroutine ALLCON



  subroutine DETON(cea)
    !***********************************************************************
    ! CHAPMAN-JOUGUET DETONATIONS.
    !***********************************************************************
    use mod_types
    use mod_legacy_io
    implicit none

    type(CEA_Core_Problem), intent(inout):: cea

    ! LOCAL VARIABLES
    character(15):: ft1  = 'T1, K'
    character(15):: fh1  = 'H1, CAL/G'
    character(15):: fhs1 = 'H1, KJ/KG'
    character(15):: fm1  = 'M1, (1/n) '
    character(15):: fg1  = 'GAMMA1'
    character(15):: fpp1 = 'P/P1'
    character(15):: ftt1 = 'T/T1'
    character(15):: fmm1 = 'M/M1'
    character(15):: frr1 = 'RHO/RHO1'
    character(15):: fdv  = 'DET VEL,M/SEC'
    character(3):: unit
    integer:: i, ii, iof, itr, mdv, mgam, mh, mmach, mp, mson, mt
    integer:: ip, it
    real(8):: a11, a12, a21, a22, alam, alpha, amm, b1, b2, d, gam, &
         p1, pp1, rk, rr1, T1, tem, tt1, ud, x1, x2

    type(CEA_Point), pointer:: p !< current point
    type(CEA_Point), pointer:: p_tmp

    integer:: i_col_start
    real(8), dimension(cea%max_points):: cpl, gm1, h1, pub, rrho, tub, V, Pcp

    cea%Eql = .true.

    if (cea%T(1) == 0) then
       cea%T(1) = cea%Rtemp(1)
       cea%Nt = 1
    end if

    iofLoop: do iof = 1, cea%Nof
       cea%Tt = cea%T(1)

       cea%Oxfl = cea%Oxf(iof)

       call NEWOF(cea)

       ! BEGIN T LOOP.
       do it = 1, cea%Nt
          T1 = cea%T(it)

          ! BEGIN P LOOP.
          do ip = 1, cea%Np
             cea%ipt = (it - 1) * cea%Np + ip

             i_col_start = int((cea%ipt - 1) / Ncol) * Ncol + 1

             p => cea%points(cea%iOF, cea%ipt)

             if (cea%ipt > 1) then
                if (mod(cea%Isv, Ncol) == 1) cea%Isv = -cea%Isv
                call SETEN(cea)
             end if

             p1 = cea%P(ip)
             cea%Tt = T1
             cea%Pp = p1

             call HCALC(cea)

             if (cea%Tt == 0) return
             if (cea%legacy_mode .and. cea%Detdbg) call OUT1(cea, IOOUT)

             h1(cea%ipt) = cea%Hsub0 * cea%R

             tub(cea%ipt) = T1
             pub(cea%ipt) = p1
             cpl(cea%ipt) = cea%Cpmix * cea%R
             cea%Tt = 3800
             pp1 = 15
             cea%Pp = pp1 * p1

             ! CALCULATE ENTHALPY FOR INITIAL ESTIMATE OF T2(TT AFTER EQLBRM)
             cea%Hsub0 = h1(cea%ipt) / cea%R + 0.75 * T1 * pp1 / cea%Wmix
             cea%Tp = .false.
             cea%Hp = .true.

             call EQLBRM(cea)

             cea%Hsub0 = h1(cea%ipt) / cea%R
             cea%Hp = .false.

             if (cea%Tt /= 0) then
                gam = p%Gammas
                tt1 = cea%Tt / T1
                ii = 0
                tem = tt1 - 0.75 * pp1 / (p%Cpr * cea%Wmix)
                amm = p%Wm / cea%Wmix

                if (cea%legacy_mode .and. cea%Detdbg) write(IOOUT, '(/" T EST.=", F8.2/11X, "P/P1", 17X, "T/T1")') cea%Tt

                ! LOOP FOR IMPROVING T2/T1 AND P2/P1 INITIAL ESTIMATE.
                do ii = 1, 3
                   alpha = amm / tt1
                   pp1 = (1 + gam) * (1 + sqrt(1 - 4 * gam * alpha / (1 + gam)**2)) / (2 * gam * alpha)
                   rk = pp1 * alpha
                   tt1 = tem + 0.5 * pp1 * gam * (rk**2 - 1) / (cea%Wmix * p%Cpr * rk)
                   if (cea%legacy_mode .and. cea%Detdbg) write(IOOUT, '(i5, 2e20.8)') ii, pp1, tt1
                end do

                cea%Tp = .true.
                cea%Tt = T1 * tt1
                rr1 = pp1 * amm / tt1

                ! BEGIN MAIN ITERATION LOOP.
                do itr = 1, 8
                   cea%Pp = p1 * pp1

                   call EQLBRM(cea)

                   if (cea%ipt == 0) exit iofLoop

                   if (cea%Tt == 0) exit

                   gam = p%Gammas
                   amm = p%Wm / cea%Wmix
                   rr1 = pp1 * amm / tt1
                   a11 = 1 / pp1 + gam * rr1 * p%Dlvpt
                   a12 = gam * rr1 * p%Dlvtp
                   a21 = 0.5 * gam * (rr1**2 - 1 - p%Dlvpt * (1 + rr1**2)) + p%Dlvtp - 1
                   a22 = -0.5 * gam * p%Dlvtp * (rr1**2 + 1) - p%Wm * p%Cpr
                   b1 = 1 / pp1 - 1 + gam * (rr1 - 1)
                   b2 = p%Wm * (p%Hsum - h1(cea%ipt) / cea%R) / cea%Tt - 0.5 * gam * (rr1**2 - 1)
                   d = a11 * a22 - a12 * a21
                   x1 = (a22 * b1 - a12 * b2) / d
                   x2 = (a11 * b2 - a21 * b1) / d
                   alam = 1
                   tem = x1

                   if (tem < 0) tem = -tem
                   if (x2 > tem) tem = x2
                   if (-x2 > tem) tem = -x2
                   if (tem > 0.4054652) alam = 0.4054652 / tem

                   pp1 = pp1 * exp(x1 * alam)
                   tt1 = tt1 * exp(x2 * alam)

                   cea%Tt = T1 * tt1
                   ud = rr1 * sqrt(R0 * gam * cea%Tt/p%Wm)

                   if (cea%legacy_mode .and. cea%Detdbg) &
                        write(IOOUT, '(/" ITER =", i2, 5x, "P/P1 =", e15.8, /7x, "T/T1 =", e15.8, 5x, &
                        & "RHO/RHO1 =", e15.8, /7x, "DEL LN P/P1 =", e15.8, 5x, &
                        & "DEL LN T/T1 =", e15.8)') itr, pp1, tt1, rr1, x1, x2

                   ! CONVERGENCE TEST
                   if (tem <= 0.5E-04) exit
                end do

                if (cea%Tt /= 0) then
                   if (itr < 8) then
                      rrho(cea%ipt) = rr1
                      if (cpl(cea%ipt) == 0) then
                         gm1(cea%ipt) = 0
                         p%Vmoc = 0
                      else
                         gm1(cea%ipt) = cpl(cea%ipt) / (cpl(cea%ipt) - cea%R / cea%Wmix)
                         p%Vmoc = ud / sqrt(R0 * gm1(cea%ipt) * T1 / cea%Wmix)
                      end if
                   else
                      if (cea%legacy_mode) write(IOOUT, '(/" CONSERVATION EQNS NOT SATISFIED IN 8 ITERATIONS (DETON)")')
                      cea%Tt = 0
                   end if

                   if (cea%Trnspt) call TRANP(cea)

                   cea%Isv = 0

                   if (ip /= cea%Np .or. it /= cea%Nt .and. cea%Tt /= 0) then
                      cea%Isv = cea%ipt
                      if (mod(cea%ipt, Ncol) /= 0) cycle
                   end if
                end if

                ! OUTPUT
                if (cea%legacy_mode) then
                   write(IOOUT, '(//, 21X, "DETONATION PROPERTIES OF AN IDEAL REACTING GAS")')
                   call OUT1(cea, IOOUT)
                end if

                ! SET MXX ARRAY FOR PLOTTING PARAMETERS
                mp    = 0
                mt    = 0
                mgam  = 0
                mh    = 0
                mdv   = 0
                mson  = 0
                mmach = 0

                do i = 1, cea%Nplt
                   if (index(cea%Pltvar(i)(2:), '1') /= 0) then
                      if (cea%Pltvar(i)(:3) == 'son') then
                         mson = i
                      else if (cea%Pltvar(i)(:3) == 'gam') then
                         mgam = i
                      else if (cea%Pltvar(i)(:1) == 'h') then
                         mh = i
                      else if (cea%Pltvar(i)(:1) == 't') then
                         mt = i
                      else if (cea%Pltvar(i)(:1) == 'p') then
                         mp = i
                      end if
                   else if (index(cea%Pltvar(i), 'vel') /= 0) then
                      mdv = i
                   else if (index(cea%Pltvar(i), 'mach') /= 0) then
                      mmach = i
                   end if
                end do

                if (cea%legacy_mode) write(IOOUT, '(/" UNBURNED GAS"/)')

                cea%fmt(4) = '13'
                cea%fmt(5) = ' '
                cea%fmt(7) = '4,'

                do i = i_col_start, cea%ipt
                   if (cea%SIunit) then
                      V(i) = pub(i)
                      unit = 'BAR'
                   else
                      V(i) = pub(i) / 1.01325d0
                      unit = 'ATM'
                   end if
                   if (mp > 0) cea%Pltout(i, mp) = V(i)
                end do

                if (cea%legacy_mode) then
                   write(IOOUT, cea%fmt) 'P1, ' // unit // '        ', (V(i), i = i_col_start, cea%ipt)

                   cea%fmt(7) = '2,'
                   write(IOOUT, cea%fmt) ft1, (tub(i), i = i_col_start, cea%ipt)

                   if (.not. cea%SIunit) write(IOOUT, cea%fmt) fh1, (h1(i), i = i_col_start, cea%ipt)
                   if (cea%SIunit) write(IOOUT, cea%fmt) fhs1, (h1(i), i = i_col_start, cea%ipt)
                end if

                do concurrent (i = i_col_start:cea%ipt)
                   p => cea%points(cea%iOF, i)
                   V(i) = cea%Wmix
                   p%Sonvel = sqrt(R0 * gm1(i) * tub(i) / cea%Wmix)
                end do

                if (cea%legacy_mode) then
                   cea%fmt(7) = '3,'
                   write(IOOUT, cea%fmt) fm1, (V(i), i = i_col_start, cea%ipt)
                   cea%fmt(7) = '4,'
                   write(IOOUT, cea%fmt) fg1, (gm1(i), i = i_col_start, cea%ipt)
                   cea%fmt(7) = '1,'
                   write(IOOUT, cea%fmt) 'SON VEL1,M/SEC ', (cea%points(cea%iOF, i)%Sonvel, i = i_col_start, cea%ipt)
                end if

                if (cea%Nplt > 0) then
                   do i = i_col_start, cea%ipt
                      p => cea%points(cea%iOF, i)
                      if (mt > 0)   cea%Pltout(i, mt) = tub(i)
                      if (mgam > 0) cea%Pltout(i, mgam) = gm1(i)
                      if (mh > 0)   cea%Pltout(i, mh) = h1(i)
                      if (mson > 0) cea%Pltout(i, mson) = p%Sonvel
                   end do
                end if

                if (cea%legacy_mode) then
                   write(IOOUT, '(/" BURNED GAS"/)')

                   cea%fmt(4) = cea%fmt(6)
                   call OUT2(cea, cea%ipt, IOOUT)

                   if (cea%Trnspt) call OUT4(cea, cea%ipt, IOOUT)

                   write(IOOUT, '(/" DETONATION PARAMETERS"/)')

                   cea%fmt(7) = '3,'
                end if

                do i = i_col_start, cea%ipt
                   p_tmp => cea%points(cea%iOF, i)
                   V(i) = p_tmp%Ppp / pub(i)
                   Pcp(i) = p_tmp%Ttt / tub(i)
                   p_tmp%Sonvel = p_tmp%Sonvel * rrho(i)
                   if (mmach > 0) cea%Pltout(i, mmach) = p_tmp%Vmoc
                   if (mdv > 0)   cea%Pltout(i, mdv) = p_tmp%Sonvel
                end do

                if (cea%legacy_mode) then
                   cea%fmt(4) = '13'
                   cea%fmt(7) = '3,'
                   write(IOOUT, cea%fmt) fpp1, V(i_col_start:cea%ipt)
                   write(IOOUT, cea%fmt) ftt1, Pcp(i_col_start:cea%ipt)
                end if

                do concurrent (i = i_col_start:cea%ipt)
                   p_tmp => cea%points(cea%iOF, i)
                   V(i) = p_tmp%Wm / cea%Wmix
                end do

                if (cea%legacy_mode) then
                   cea%fmt(7) = '4,'
                   write(IOOUT, cea%fmt) fmm1, (V(i), i = i_col_start, cea%ipt)
                   write(IOOUT, cea%fmt) frr1, (rrho(i), i = i_col_start, cea%ipt)
                   write(IOOUT, cea%fmt) 'DET MACH NUMBER', (cea%points(cea%iOF, i)%Vmoc, i = i_col_start, cea%ipt)

                   cea%fmt(7) = '1,'
                   write(IOOUT, cea%fmt) fdv, (cea%points(cea%iOF, i)%Sonvel, i = i_col_start, cea%ipt)
                end if

                cea%Eql = .true.

                if (cea%legacy_mode) call OUT3(cea, cea%ipt, IOOUT)

                cea%Iplt = min(cea%ipt, 500)

                if (cea%Isv == 0 .and. iof == cea%Nof) exit iofLoop
                if (cea%Np == 1 .and. cea%Nt == 1) cycle iofLoop

                if (cea%legacy_mode) write(IOOUT, '(///)')
             end if
          end do
       end do

       cea%Iplt = min(cea%ipt - 1, 500)

    end do iofLoop

    cea%Tp = .false.

    return
  end subroutine DETON


  subroutine EQLBRM(cea)
    !***********************************************************************
    ! CALCULATE EQUILIBRIUM COMPOSITION AND PROPERTIES.
    !***********************************************************************
    use mod_types
    use mod_general
    implicit none

    type(CEA_Core_Problem), intent(inout):: cea

    ! LOCAL VARIABLES
    character(12):: ae, cmp(maxEl)
    character(16):: amb
    logical:: cpcalc, i2many, reduce
    integer:: i, il, ilamb, ilamb1, inc, ipr, iq2, iter, ixsing, iz, j, ja, jb, &
         jbx, jc, jcondi, jcons, jdelg, jkg, jneg, jsw, k, kc, kg, kk, &
         kmat, kneg, l, lc, lcs(maxEl), le, lelim, lk, ll, lncvg, ls, lsing, &
         lz, maxitn, ncvg, njc, nn, numb
    real(8):: aa, ambda, ambda1, bigen, bigneg, delg, dlnt, dpie, esize
    real(8):: gap, gasfrc, pie, siz9, sizeg
    real(8):: sum0, sum1, szgj, tem, tsize, ween, xi, xln, xsize
    real(8):: smalno = 1e-6, smnol = -13.815511
    logical:: mask(cea%Ng)
    integer:: Npt

    type(CEA_Point), pointer:: p !< current point
    p => cea%points(cea%iOF, cea%ipt)

    Npt = mod(cea%ipt - 1, Ncol) + 1
    if (cea%Rkt) then
       if (cea%Fac) then
          Npt = Npt + 3 * ((cea%ipt - 1) / Ncol)
       else
          Npt = Npt + 2 * ((cea%ipt - 1) / Ncol)
       end if
    end if

    ixsing = 0
    lsing = 0
    jsw = 0
    jdelg = 0
    maxitn = 50
    ncvg = 0
    lncvg = 3 * cea%Nlm
    reduce = .false.
    siz9 = cea%Size - 9.2103404d0
    tsize = cea%Size
    xsize = cea%Size + 6.90775528d0

    if (cea%Trace /= 0.) then
       maxitn = maxitn + cea%Ngc / 2
       xsize = -log(cea%Trace)
       if (xsize < cea%Size) xsize = cea%Size + 0.1
    end if

    if (xsize > 80) xsize = 80

    esize = min(80d0, xsize + 6.90775528d0)
    jcons = 0
    pie = 0
    i2many = .false.
    cea%Pderiv = .false.
    cea%Convg = .false.
    numb = 0
    cpcalc = .true.

    if (cea%Tp) cpcalc = .false.

    if (cea%Tt /= 0) then
       if (cea%Npr == 0 .or. (cea%Tt /= cea%T(1) .and. .not. cea%Tp)) go to 400
       k = 1
    else
       cea%Tt = 3800
       if (cea%Npr == 0) go to 400
       k = 1
    end if

100 j = cea%Jcond(k)
    jc = j - cea%Ng
    kg = -cea%Ifz(jc)

    do i = 1, 9
       kg = kg + 1
       kc = jc + kg

       if (cea%Tt <= cea%Temp(2, kc)) then
          if (kg /= 0) then
             cea%Jcond(k) = j + kg
             p%En(j+kg) = p%En(j)
             p%En(j) = 0

             if (cea%legacy_mode .and. cea%Prod(j) /= cea%Prod(j+kg) .and. .not. cea%Short) &
                  write(IOOUT, '(" PHASE CHANGE, REPLACE ", A16, "WITH ", A16)') cea%Prod(j), cea%Prod(j+kg)
          end if
          go to 300
       else if (kc >= cea%Nc .or. cea%Ifz(kc+1) <= cea%Ifz(kc)) then
          exit
       end if
    end do

    if (.not. cea%Tp) then
       cea%Tt = cea%Temp(2, kc) - 10
       k = 1
       go to 100
    end if

    if (cea%legacy_mode) write(IOOUT, '(" REMOVE ", A16)') cea%Prod(j)

    p%En(j) = 0
    cea%Enln(j) = 0
    cea%Deln(j) = 0

    do i = k, cea%Npr
       cea%Jcond(i) = cea%Jcond(i+1)
    end do

    cea%Npr = cea%Npr - 1

300 k = k + 1

    if (k <= cea%Npr) go to 100

400 cea%Tln = log(cea%Tt)

    if (cea%Vol) cea%Pp = R0 * cea%Enn * cea%Tt / cea%Vv

    call CPHS(cea)

    cea%Tm = log(cea%Pp / cea%Enn)
    le = cea%Nlm

    if (cea%Lsave /= 0 .and. cea%Nlm /= cea%Lsave) then
       tem = exp(-tsize)
       do i = cea%Lsave + 1, cea%Nlm
          do concurrent (j = 1:cea%Ng, cea%A(i, j) /= 0)
             p%En(j) = tem
             cea%Enln(j) = -tsize
          end do
       end do
    end if

    ls = cea%Nlm
    lelim = 0
    lz = ls

    if (cea%Ions) lz = ls - 1

    if (cea%legacy_mode .and. Npt == 1 .and. .not. cea%Shock .and. .not. cea%Short) then
       write(IOOUT, '(/" POINT ITN", 6X, "T", 10X, 4(A4, 8X)/(18X, 5(A4, 8X)))') (cea%Elmt(i), i = 1, cea%Nlm)
    end if

    if (p%Debug) then
       do i = 1, cea%Nlm
          cmp(i) = cea%Elmt(i)
       end do
    end if

    ! BEGIN ITERATION
500 if (cpcalc) then
       cea%Cpsum = sum(p%En(1:cea%Ng) * cea%Cp(1:cea%Ng))

       if (cea%Npr /= 0) then
          cea%Cpsum = cea%Cpsum + sum(p%En(cea%Jcond(1:cea%Npr)) * cea%Cp(cea%Jcond(1:cea%Npr)))
          cpcalc = .false.
       end if
    end if

    numb = numb + 1

    call MATRIX(cea)

    iq2 = cea%Iq1 + 1

    if (cea%Convg) cea%Imat = cea%Imat - 1

    if (cea%legacy_mode .and. p%Debug) then
       if (.not. cea%Convg) then
          write(IOOUT, '(/" ITERATION", I3, 6X, "MATRIX ")') numb
       else
          if (.not. cea%Pderiv) write(IOOUT, '(/" T DERIV MATRIX")')
          if (cea%Pderiv) write(IOOUT, '(/" P DERIV MATRIX")')
       end if

       kmat = cea%Imat + 1

       do i = 1, cea%Imat
          write(IOOUT, '(3X, 5E15.6)') (cea%G(i, k), k = 1, kmat)
       end do
    end if

    call gauss_elimination(cea%G(1:cea%Imat, 1:cea%Imat+1), cea%X(1:cea%Imat), cea%Msing)

    if (cea%Msing == 0) then
       if (cea%legacy_mode .and. p%Debug) then
          write(IOOUT, '(/" SOLUTION VECTOR", /, 6x, 5A15/8X, 5A15)') (cmp(k), k = 1, le)
          write(IOOUT, '(3X, 5E15.6)') (cea%X(i), i = 1, cea%Imat)
       end if

       if (.not. cea%Convg) then
          ! OBTAIN CORRECTIONS TO THE ESTIMATES
          if (cea%Vol) cea%X(iq2) = cea%X(cea%Iq1)
          if (cea%Tp) cea%X(iq2) = 0
          dlnt = cea%X(iq2)
          sum0 = cea%X(cea%Iq1)

          if (cea%Vol) then
             cea%X(cea%Iq1) = 0
             sum0 = -dlnt
          end if

          outerLoop0: do j = 1, cea%Ng
             if (lelim /= 0) then
                cea%Deln(j) = 0
                do i = lelim, ls
                   if (cea%A(i, j) /= 0) cycle outerLoop0
                end do
             end if

             cea%Deln(j) = -cea%Mu(j) + cea%H0(j) * dlnt + sum0

!!$           cea%Deln(j) = cea%Deln(j) + sum(A(1:cea%Nlm, j) * X(1:cea%Nlm)) ???
             do k = 1, cea%Nlm
                cea%Deln(j) = cea%Deln(j) + cea%A(k, j) * cea%X(k)
             end do

             if (pie /= 0) cea%Deln(j) = cea%Deln(j) + cea%A(ls, j) * pie
          end do outerLoop0

          if (cea%Npr /= 0) then
             do concurrent (k = 1:cea%Npr)
                cea%Deln(cea%Jcond(k)) = cea%X(cea%Nlm+k)
             end do
          end if

          ! CALCULATE CONTROL FACTOR, AMBDA
          ambda = 1
          ambda1 = 1
          ilamb = 0
          ilamb1 = 0
          sum0 = 5 * max(abs(cea%X(cea%Iq1)), abs(dlnt))

          do j = 1, cea%Ng
             if (cea%Deln(j) > 0) then
                if ((cea%Enln(j) - cea%Ennl + cea%Size) <= 0) then

                   sum1 = abs(cea%Deln(j)-cea%X(cea%Iq1))

                   if (sum1 >= siz9) then
                      sum1 = abs(-9.2103404d0 - cea%Enln(j) + cea%Ennl) / sum1

                      if (sum1 < ambda1) then
                         ambda1 = sum1
                         ilamb1 = j
                      end if
                   end if

                else if (cea%Deln(j) > sum0) then
                   sum0 = cea%Deln(j)
                   ilamb = j
                end if
             end if
          end do

          if (sum0 > 2) ambda = 2 / sum0

          if (ambda1 <= ambda) then
             ambda = ambda1
             ilamb = ilamb1
          end if

          if (cea%legacy_mode .and. p%Debug) then
             ! INTERMEDIATE OUTPUT
             write(IOOUT, '(/" T=", E15.8, " ENN=", E15.8, " ENNL=", E15.8, " PP=", E15.8, &
                  & /" LN P/N=", E15.8, " AMBDA=", E15.8)') cea%Tt, cea%Enn, cea%Ennl, cea%Pp, cea%Tm, ambda

             if (ambda /= 1) then
                amb = 'ENN'
                if (abs(cea%X(iq2)) > abs(cea%X(cea%Iq1))) amb = 'TEMP'
                if (ilamb /= 0) amb = cea%Prod(ilamb)
                write(IOOUT, '(/" AMBDA SET BY ", A16)') amb
             end if

             if (cea%Vol) write(IOOUT, '(" VOLUME=", E15.8, "CC/G")') cea%Vv * 0.001d0

             write(IOOUT, '(/24X, "Nj", 12X, "LN Nj", 8X, "DEL LN Nj", 6X, "H0j/RT", /, 41X, &
                  & "S0j/R", 10X, " G0j/RT", 8X, " Gj/RT")')

             do j = 1, cea%Ngc
                write(IOOUT, '(1X, A16, 4E15.6, /35x, 3E15.6)') &
                     cea%Prod(j), p%En(j), cea%Enln(j), cea%Deln(j), cea%H0(j), cea%S(j), cea%H0(j) - cea%S(j), cea%Mu(j)
             end do
          end if

          ! APPLY CORRECTIONS TO ESTIMATES
          p%Totn = 0

          do concurrent (j = 1:cea%Ng)
             cea%Enln(j) = cea%Enln(j) + ambda * cea%Deln(j)
          end do

          p%En(1:cea%Ng) = 0

          if (lelim == 0) then
             do concurrent (j = 1:cea%Ng, (cea%Enln(j) - cea%Ennl + tsize) > 0)
                p%En(j) = exp(cea%Enln(j))
             end do
          else
             do concurrent (j = 1:cea%Ng, all(cea%A(lelim:ls, j) == 0))
                p%En(j) = exp(cea%Enln(j))
             end do
          end if

          p%Totn = p%Totn + sum(p%En(1:cea%Ng))

          if (cea%Ions .and. cea%Elmt(cea%Nlm) == 'E') then
             mask = .false.
             do concurrent (j = 1:cea%Ng, cea%A(ls, j) /= 0 .and. p%En(j) == 0 .and. (cea%Enln(j) - cea%Ennl + esize) > 0)
                p%En(j) = exp(cea%Enln(j))
                mask(j) = .true.
             end do
             p%Totn = p%Totn + sum(p%En(1:cea%Ng), mask=mask)
          end if

          cea%Sumn = p%Totn

          if (cea%Npr /= 0) then
             do concurrent (k = 1:cea%Npr)
                p%En(cea%Jcond(k)) = p%En(cea%Jcond(k)) + ambda * cea%Deln(cea%Jcond(k))
             end do
             p%Totn = p%Totn + sum(p%En(cea%Jcond(1:cea%Npr)))
          end if

          if (.not. cea%Tp) then
             cea%Tln = cea%Tln + ambda * dlnt
             cea%Tt = exp(cea%Tln)
             cpcalc = .true.
             call CPHS(cea)
          end if

          if (cea%Vol) then
             cea%Enn = cea%Sumn
             cea%Ennl = log(cea%Enn)
             if (cea%Vol) cea%Pp = R0 * cea%Tt * cea%Enn / cea%Vv
          else
             cea%Ennl = cea%Ennl + ambda * cea%X(cea%Iq1)
             cea%Enn = exp(cea%Ennl)
          end if

          cea%Tm = log(cea%Pp / cea%Enn)

          if (cea%Elmt(cea%Nlm) == 'E') then
             ! CHECK ON REMOVING IONS
             if (all(cea%A(cea%Nlm, 1:cea%Ngc) == 0 .or. p%En(1:cea%Ngc) <= 0)) then
                pie = cea%X(cea%Nlm)
                lelim = cea%Nlm
                cea%Nlm = cea%Nlm - 1
                go to 500
             end if
          end if

          ! TEST FOR CONVERGENCE
          if (numb > maxitn) then
             if (cea%legacy_mode) write(IOOUT, '(/, I4, " ITERATIONS DID NOT SATISFY CONVERGENCE", /, 15x, &
                  & " REQUIREMENTS FOR THE POINT", I5, " (EQLBRM)")') maxitn, Npt

             if (cea%Nc == 0 .or. i2many) go to 1500

             i2many = .true.

             if (.not. cea%Hp .or. Npt /= 1 .or. cea%Tt > 100.) then
                if (cea%Npr /= 1 .or. cea%Enn > 1.E-4) go to 1500
                ! HIGH TEMPERATURE, INCLUDED CONDENSED CONDITION
                if (cea%legacy_mode) write(IOOUT, '(/" TRY REMOVING CONDENSED SPECIES (EQLBRM)")')

                cea%Enn = 0.1
                cea%Ennl = -2.3025851
                cea%Sumn = cea%Enn
                xi = cea%Enn / cea%Ng
                xln = log(xi)

                do concurrent (j = 1:cea%Ng)
                   p%En(j) = xi
                   cea%Enln(j) = xln
                end do

                j = cea%Jcond(1)
                k = 1
                go to 1000

             else
                if (cea%legacy_mode) write(IOOUT, '(/" LOW TEMPERATURE IMPLIES A CONDENSED SPECIES SHOULD HA", &
                     & "VE BEEN INSERTED,", &
                     & /" RESTART WITH insert DATASET (EQLBRM)")')
                go to 1500
             end if

          else
             if (abs(cea%X(cea%Iq1) * cea%Enn / p%Totn) > 0.5E-5 .or. &
                  any(abs(cea%Deln(1:cea%Ng)) * p%En(1:cea%Ng) / p%Totn > 0.5d-5) .or. &
                  abs(dlnt) > 1.d-04) go to 500

             if (cea%Npr /= 0) then
                do k = 1, cea%Npr
                   j = cea%Jcond(k)
                   if (abs(cea%Deln(j) / p%Totn) > 0.5d-5) go to 500
                   if (p%En(j) < 0) go to 700
                end do
             end if

             le = cea%Nlm

             if (any(abs(cea%B0(1:cea%Nlm)) >= 1.d-06 .and. &
                  abs(cea%B0(1:cea%Nlm) - sum(spread(p%En(1:cea%Ngc), 1, cea%Nlm) * cea%A(1:cea%Nlm, 1:cea%Ngc), DIM=2)) &
                  > cea%Bcheck)) go to 500

             if (cea%Trace /= 0.) then
                tsize = xsize
                tem = 1

                if (numb /= 1) then
                   lk = lz
                   if (cea%Nlm < lz) lk = cea%Nlm
                   do i = 1, lk
                      if (i /= lsing) then
                         tem = 0
                         if (cea%X(i) /= 0) then
                            tem = abs((cea%pisave(i) - cea%X(i)) / cea%X(i))
                            if (tem > 0.001) exit
                         end if
                      end if
                   end do
                end if

                do concurrent (i = 1:cea%Nlm)
                   cea%pisave(i) = cea%X(i)
                end do

                if (tem > 0.001) go to 500

                if (cea%Ions) then
                   ! CHECK ON ELECTRON BALANCE
                   iter = 1

                   if (pie /= 0) then
                      cea%X(cea%Nlm+1) = pie
                   end if

566                sum1 = 0
                   sum0 = 0
                   pie = cea%X(le)

                   do j = 1, cea%Ng
                      if (cea%A(ls, j) /= 0) then
                         p%En(j) = 0
                         tem = 0

                         if (cea%Enln(j) > -87) tem = exp(cea%Enln(j))

                         if ((cea%Enln(j)-cea%Ennl+tsize) > 0 .and. cea%Elmt(cea%Nlm) == 'E') then
                            pie = 0
                            p%En(j) = tem
                         end if

                         aa = cea%A(ls, j) * tem
                         sum0 = sum0 + aa
                         sum1 = sum1 + aa * cea%A(ls, j)
                      end if
                   end do

                   if (sum1 /= 0) then
                      dpie = -sum0 / sum1

                      do concurrent (j = 1:cea%Ng, cea%A(ls, j) /= 0)
                         cea%Enln(j) = cea%Enln(j) + cea%A(ls, j) * dpie
                      end do

                      if (cea%legacy_mode .and. p%Debug) &
                           write(IOOUT, '(/" ELECTRON BALANCE ITER NO. =", i4, "  DELTA PI =", e14.7)') iter, dpie

                      if (abs(dpie) > 0.0001) then
                         cea%X(le) = cea%X(le) + dpie
                         iter = iter + 1

                         if (iter <= 80) go to 566
                         if (cea%legacy_mode) write(IOOUT, '(/" DID NOT CONVERGE ON ELECTRON BALANCE (EQLBRM)")')
                         go to 1500

                      else if (cea%Elmt(cea%Nlm) == 'E' .and. pie /= 0) then
                         cea%Nlm = cea%Nlm - 1
                         cea%newcom = .true.
                      end if
                   end if
                end if
             end if
          end if

       else if (.not. cea%Pderiv) then
          ! TEMPERATURE DERIVATIVES--CONVG=T, PDERIV=F
          p%Dlvtp = 1. - cea%X(cea%Iq1)
          p%Cpr = cea%G(iq2, iq2)

          p%Cpr = p%Cpr - sum(cea%G(iq2, 1:cea%Iq1) * cea%X(1:cea%Iq1))

          ! PRESSURE DERIVATIVE--CONVG=T, PDERIV=T
          cea%Pderiv = .true.
          go to 500

       else
          p%Dlvpt = -1 + cea%X(cea%Iq1)
          if (cea%Jliq == 0) then
             p%Gammas = -1 / (p%Dlvpt + (p%Dlvtp**2) * cea%Enn / p%Cpr)
          else
             p%En(cea%Jsol) = cea%ensol
             p%Hsum = p%Hsum + p%En(cea%Jliq) * (cea%H0(cea%Jliq) - cea%H0(cea%Jsol))
             p%Gammas = -1. / p%Dlvpt
             cea%Npr = cea%Npr + 1
             cea%Jcond(cea%Npr) = cea%Jliq
          end if
          go to 1400
       end if

       ! SINGULAR MATRIX
    else
       if (cea%Convg) then
          if (cea%legacy_mode) write(IOOUT, '(/" DERIVATIVE MATRIX SINGULAR (EQLBRM)")')
          p%Dlvpt = -1
          p%Dlvtp = 1
          p%Cpr = cea%Cpsum
          p%Gammas = -1 / (p%Dlvpt + p%Dlvtp**2 * cea%Enn / p%Cpr)
          go to 1400

       else
          if (cea%legacy_mode) write(IOOUT, '(/" SINGULAR MATRIX, ITERATION", I3, "  VARIABLE", I3, "(EQLBRM)")') numb, cea%Msing
          lsing = cea%Msing
          ixsing = ixsing + 1
          if (ixsing <= 8) then
             xsize = 80
             tsize = xsize
             if (cea%Msing > cea%Nlm .and. numb < 1 .and. cea%Npr > 1 .and. jdelg > 0) then
                ween = 1000
                j = 0

                do i = 1, cea%Npr
                   jcondi = cea%Jcond(i)
                   if (jcondi /= jdelg) then
                      do ll = 1, cea%Nlm
                         if (cea%A(ll, jdelg) /= 0 .and. cea%A(ll, jcondi) /= 0) then
                            if (p%En(jcondi) <= ween) then
                               ween = p%En(jcondi)
                               j = jcondi
                               k = i
                            end if
                            exit
                         end if
                      end do
                   end if
                end do

                if (j > 0) then
                   if (cea%legacy_mode) write(IOOUT, '(/" TRY REMOVING CONDENSED SPECIES (EQLBRM)")')
                   go to 1000
                end if

             else if (.not. cea%Hp .or. Npt /= 1 .or. cea%Nc == 0 .or. cea%Tt > 100) then
                if (ixsing >= 3) then
                   if (cea%Msing < cea%Iq1) then
                      if (reduce .and. cea%Msing <= cea%Nlm) then
                         if (cea%Nlm < lelim) go to 1500
                         if (cea%legacy_mode) write(IOOUT, '(/" WARNING!! POINT", I3, &
                              & " USES A REDUCED SET OF COMPONENTS", / &
                              & " SPECIES CONTAINING THE ELIMINATED COMPONENT ARE OMITTED.", &
                              & / &
                              & " IT MAY BE NECESSARY TO RERUN WITH INSERTED CONDENSED SPECIES", &
                              & /" CONTAINING COMPONENT ", A8, "(EQLBRM)")') Npt, cea%Elmt(cea%Nlm)
                         cea%Nlm = cea%Nlm - 1
                         go to 500

                      else if (cea%Msing <= cea%Nlm) then
                         ! FIND NEW COMPONENTS
                         if (.not. cea%Ions) go to 1100
                         if (cea%Elmt(cea%Nlm) /= 'E') go to 1100

                         do concurrent (j = 1:cea%Ng, cea%A(cea%Nlm, j) /= 0)
                            p%En(j) = 0
                         end do

                         pie = cea%X(cea%Nlm)
                         cea%Nlm = cea%Nlm - 1
                         if (cea%Msing > cea%Nlm) go to 500
                         go to 1100
                      else
                         ! REMOVE CONDENSED SPECIES TO CORRECT SINGULARITY
                         k = cea%Msing - cea%Nlm
                         j = cea%Jcond(k)

                         if (j /= jcons) then
                            jcons = j
                            go to 1000
                         end if
                      end if
                   end if
                end if

                do concurrent (j = 1:cea%Ng, &
                     .not. (cea%Ions .and. cea%Elmt(cea%Nlm) /= 'E' .and. cea%A(ls, j) /= 0) .and. p%En(j) == 0)
                   p%En(j) = smalno
                   cea%Enln(j) = smnol
                end do

                go to 500
             else
                if (cea%legacy_mode) write(IOOUT, '(/" LOW TEMPERATURE IMPLIES A CONDENSED SPECIES SHOULD HA", &
                     & "VE BEEN INSERTED,", &
                     & /" RESTART WITH insert DATASET (EQLBRM)")')
             end if
          end if
       end if

       go to 1500
    end if

    ! CALCULATE ENTROPY, CHECK ON DELTA S FOR SP PROBLEMS
600 p%Ssum = sum(p%En(1:cea%Ng) * (cea%S(1:cea%Ng) - cea%Enln(1:cea%Ng) - cea%Tm))

    if (cea%Npr > 0) then
       p%Ssum = p%Ssum + sum(p%En(cea%Jcond(1:cea%Npr)) * cea%S(cea%Jcond(1:cea%Npr)))
    end if

    if (.not. cea%Sp) then
       cea%Convg = .true.
    else
       tem = p%Ssum - cea%S0
       if (abs(tem) > 0.0005) go to 500
       if (cea%legacy_mode .and. p%Debug) write(IOOUT, '(/" DELTA S/R =", e15.8)') tem
       cea%Convg = .true.
    end if

    ! CONVERGENCE TESTS ARE SATISFIED, TEST CONDENSED SPECIES.
700 ncvg = ncvg + 1

    if (ncvg > lncvg) then
       ! ERROR, SET TT=0
       if (cea%legacy_mode) write(IOOUT, '(/, I3, " CONVERGENCES FAILED TO ESTABLISH SET OF CONDENSED", " SPECIES (EQLBRM)")') lncvg
       go to 1500
    else
       if (.not. cea%Shock) then
          if (cea%legacy_mode .and. .not. cea%Short) then
             if (cea%newcom) write(IOOUT, '(/" POINT ITN", 6x, "T", 10x, 4a12/(18x, 5a12))') (cmp(k), k = 1, le)
             write(IOOUT, '(i4, i5, 5f12.3, /(12x, 5f12.3))') Npt, numb, cea%Tt, (cea%X(il), il = 1, le)
          end if

          if (.not. cea%Tp .and. cea%Npr == 0 .and. cea%Tt <= cea%Tg(1) * 0.2d0) then
             if (cea%legacy_mode) write(IOOUT, '(/" LOW TEMPERATURE IMPLIES A CONDENSED SPECIES SHOULD HA", "VE BEEN INSERTED,", &
                  & /" RESTART WITH insert DATASET (EQLBRM)")')
             go to 1500
          end if

          cea%newcom = .false.
       end if

       if (cea%Npr /= 0) then
          bigneg = 0
          jneg = 0

          do k = 1, cea%Npr
             j = cea%Jcond(k)
             if (p%En(j) * cea%Cp(j) <= bigneg) then
                bigneg = p%En(j) * cea%Cp(j)
                jneg = j
                kneg = k
             end if
          end do

          if (jneg /= 0) then
             j = jneg
             k = kneg
             if (j == cea%Jsol .or. j == cea%Jliq) then
                cea%Jsol = 0
                cea%Jliq = 0
             end if
             go to 1000
          end if
       end if

       if (cea%Ngc /= cea%Ng .or. cea%Tp) then
          cea%Ng = cea%Ngc

          call CPHS(cea)

          cea%Ng = cea%Ngp1 - 1
          cpcalc = .true.

          if (cea%Ngc == cea%Ng) go to 750

          call ALLCON(cea)

          if (cea%Npr /= 0 .and. .not. cea%Tp) then
             gap = 50

             outerLoop1: do ipr = 1, cea%Npr
                j = cea%Jcond(ipr)
                if (j /= cea%Jsol .and. j /= cea%Jliq) then
                   inc = j - cea%Ng
                   kg = -cea%Ifz(inc)

                   do iz = 1, 20
                      kg = kg + 1
                      kc = inc + kg

                      if (cea%Tt <= cea%Temp(2, kc)) then
                         if (kg /= 0) then
                            jkg = j + kg
                            if (abs(kg) > 1 .or. cea%Prod(j) == cea%Prod(jkg)) go to 740
                            if (jkg == jsw) go to 720
                            if (cea%Tt < cea%Temp(1, inc) - gap .or. cea%Tt > cea%Temp(2, inc) + gap) go to 740
                            go to 720
                         end if
                         cycle outerLoop1
                      else if (cea%Ifz(kc+1) <= cea%Ifz(kc)) then
                         cycle outerLoop1
                      end if
                   end do
                   if (cea%Tt > cea%Temp(2, kc) * 1.2d0) go to 1000
                end if
             end do outerLoop1
          end if

          sizeg = 0
          szgj = 0

          do inc = 1, cea%Nc
             j = inc + cea%Ng

             if (cea%legacy_mode .and. p%Debug) then
                write(IOOUT, '(/1x, a15, 2f10.3, 3x, e15.7)') cea%Prod(j), cea%Temp(1, inc), cea%Temp(2, inc), p%En(j)
             end if

             if (p%En(j) <= 0) then
                if (cea%Tt > cea%Temp(1, inc) .or. cea%Temp(1, inc) == cea%Tg(1)) then
                   if (cea%Tt <= cea%Temp(2, inc)) then
                      delg = (cea%H0(j) - cea%S(j) - sum(cea%A(1:cea%Nlm, j) * cea%X(1:cea%Nlm))) / cea%Mw(j)

                      if (delg < sizeg .and. delg < 0) then
                         if (j /= jcons) then
                            sizeg = delg
                            jdelg = j
                         else
                            szgj = delg
                         end if
                         ipr = ipr - 1
                      end if

                      if (cea%legacy_mode .and. p%Debug) then
                         write(IOOUT, '(" [G0j-SUM(Aij*PIi)]/Mj =", E15.7, 9X, "MAX NEG DELTA G =", E15.7)') delg, sizeg
                      end if
                   end if
                end if
             end if
          end do

          if (sizeg == 0 .and. szgj == 0) go to 750

          if (sizeg /= 0) then
             j = jdelg
             go to 800
          else
             if (cea%legacy_mode) &
                  write(IOOUT, '(/" REINSERTION OF ", A16, " LIKELY TO CAUSE SINGULARITY, ", "(EQLBRM)")') cea%Prod(jcons)
             go to 1500
          end if

720       kk = max(0, kg)
          cea%Tt = cea%Temp(kk+1, inc)
          cea%Tln = log(cea%Tt)
          cea%Jsol = min(j, jkg)
          cea%Jliq = cea%Jsol + 1
          p%En(jkg) = 0.5d0 * p%En(j)
          p%En(j) = p%En(jkg)
          j = jkg
          go to 800

          ! WRONG PHASE INCLUDED FOR T INTERVAL, SWITCH EN
740       p%En(jkg) = p%En(j)
          cea%Jcond(ipr) = jkg
          p%En(j) = 0
          jsw = j

          if (cea%legacy_mode .and. cea%Prod(j) /= cea%Prod(jkg) .and. .not. cea%Short) &
               write(IOOUT, '(" PHASE CHANGE, REPLACE ", A16, "WITH ", A16)') cea%Prod(j), cea%Prod(jkg)

          j = jkg
          go to 900
       end if

       ! CONVERGED WITH NO CONDENSED CHANGES.  IF BOTH SOLID & LIQ PRESENT, 
       ! TEMPORARILY REMOVE LIQ TO PREVENT SINGULAR DERIVATIVE MATRIX.
750    cea%Sumn = cea%Enn
       if (cea%Jsol /= 0) then
          cea%ensol = p%En(cea%Jsol)
          p%En(cea%Jsol) = p%En(cea%Jsol) + p%En(cea%Jliq)
          p%Dlvtp = 0
          p%Cpr = 0
          p%Gammas = 0
          cea%Pderiv = .true.

          do k = 1, cea%Npr
             if (cea%Jcond(k) == cea%Jliq) exit
          end do

          cea%Jcond(k:cea%Npr) = cea%Jcond(k+1:cea%Npr+1)

          cea%Npr = cea%Npr - 1
       end if

       go to 500
    end if

    ! ADD CONDENSED SPECIES
800 cea%Npr = cea%Npr + 1

    cea%Jcond(2:cea%Npr) = cea%Jcond(1:cea%Npr-1)
    cea%Jcond(1) = j

    if (cea%legacy_mode .and. .not. cea%Short) write(IOOUT, '(" ADD ", a16)') cea%Prod(j)

900 inc = j - cea%Ng
    cea%Convg = .false.
    if (cea%Tp) cpcalc = .false.
    numb = -1
    go to 500

    ! REMOVE CONDENSED SPECIES
1000 p%En(j) = 0
    cea%Deln(j) = 0
    cea%Enln(j) = 0

    cea%Jcond(k:cea%Npr) = cea%Jcond(k+1:cea%Npr+1)

    if (cea%legacy_mode .and. .not. cea%Short) write(IOOUT, '(" REMOVE ", A16)') cea%Prod(j)

    cea%Npr = cea%Npr - 1
    do i = 1, cea%Nlm
       if (cmp(i) == cea%Prod(j)) then
          numb = -1
          cea%Convg = .false.
          if (cea%Tp) cpcalc = .false.
          go to 1100
       end if
    end do

    go to 900

1100 cea%newcom = .false.
    nn = cea%Nlm

    if (cea%Elmt(cea%Nlm) == 'E') nn = cea%Nlm - 1

    ! FIND ORDER OF SPECIES FOR COMPONENTS - BIGGEST TO SMALLEST
    njc = 0
    lcs(1:nn) = 0

1200 bigen = -1d-35

    do j = 1, cea%Ng
       if (p%En(j) > bigen) then
          if (.not. cea%Ions .or. cea%A(ls, j) == 0) then
             bigen = p%En(j)
             jbx = j
          end if
       end if
    end do

    if (bigen > 0.) then
       do lc = 1, nn
          if (jbx == 0) jbx = cea%Jx(lc)

          if (cea%A(lc, jbx) > smalno) then
             if (njc /= 0) then
                do i = 1, njc
                   l = lcs(i)
                   if (l == lc) cycle
                   if (l == 0) exit
                   j = cea%Jcm(l)
                   if (all(cea%A(1:nn, jbx) == cea%A(1:nn, j))) cycle
                end do
             end if

             do i = 1, nn
                if (i /= lc .and. abs(cea%A(lc, jbx) * cea%A(i, cea%jx(i)) - cea%A(lc, cea%jx(i)) * cea%A(i, jbx)) <= smalno) then
                   cycle
                end if
             end do

             njc = njc + 1
             if (jbx /= cea%Jcm(lc)) cea%newcom = .true.
             cea%Jcm(lc) = jbx
             lcs(njc) = lc
             exit
          end if
       end do

       p%En(jbx) = -p%En(jbx)
       if (njc < nn) go to 1200
    end if

    do concurrent (j = 1:cea%Ng)
       p%En(j) = abs(p%En(j))
    end do

    if (cea%newcom) then
       ! SWITCH COMPONENTS
       do lc = 1, nn
          jb = cea%Jcm(lc)

          if (cea%A(lc, jb) == 0) then
             jb = cea%Jx(lc)
             cea%Jcm(lc) = jb
          end if

          tem = cea%A(lc, jb)

          if (tem /= 0) then
             cea%pisave(lc) = cea%H0(jb) - cea%S(jb)

             if (jb <= cea%Ng) cea%pisave(lc) = cea%pisave(lc) + cea%Enln(jb) + cea%Tm
             cmp(lc) = trim(cea%Prod(jb))

             ! CALCULATE NEW COEFFICIENTS
             if (tem /= 1) then
                cea%B0(lc) = cea%B0(lc) / tem
                cea%B0p(lc, 1) = cea%B0p(lc, 1) / tem
                cea%B0p(lc, 2) = cea%B0p(lc, 2) / tem

                do concurrent (j = 1:cea%Nspx)
                   cea%A(lc, j) = cea%A(lc, j) / tem
                end do
             end if

             do i = 1, nn
                if (cea%A(i, jb) /= 0. .and. i /= lc) then
                   tem = cea%A(i, jb)

                   do concurrent (j = 1:cea%Nspx)
                      cea%A(i, j) = cea%A(i, j) - cea%A(lc, j) * tem
                   end do

                   do concurrent (j = 1:cea%Nspx, abs(cea%A(i, j)) < 1.E-5)
                      cea%A(i, j) = 0
                   end do

                   cea%B0(i) = cea%B0(i) - cea%B0(lc) * tem
                   cea%B0p(i, 1) = cea%B0p(i, 1) - cea%B0p(lc, 1) * tem
                   cea%B0p(i, 2) = cea%B0p(i, 2) - cea%B0p(lc, 2) * tem
                end if
             end do
          end if
       end do

       if (cea%legacy_mode .and. p%Debug) then
          write(IOOUT, '(/" NEW COMPONENTS")')
          write(IOOUT, '(/2x, 6A12)') (cmp(k), k = 1, nn)
       end if
    end if

    if (cea%Msing /= 0) then
       ! SWITCH ORDER OF MSING AND NLM COMPONENTS
       reduce = .true.
       lelim = cea%Nlm
       lsing = cea%Nlm

       if (cea%Msing /= cea%Nlm) then

          do j = 1, cea%Nspx
             aa = cea%A(cea%Msing, j)
             cea%A(cea%Msing, j) = cea%A(cea%Nlm, j)
             cea%A(cea%Nlm, j) = aa
          end do

          ja = cea%Jcm(cea%Msing)
          cea%Jcm(cea%Msing) = cea%Jcm(cea%Nlm)
          cea%Jcm(cea%Nlm) = ja
          ae = cmp(cea%Msing)
          cmp(cea%Msing) = cmp(cea%Nlm)
          cmp(cea%Nlm) = ae
          ae = cea%Elmt(cea%Msing)
          cea%Elmt(cea%Msing) = cea%Elmt(cea%Nlm)
          cea%Elmt(cea%Nlm) = trim(ae)
          ja = cea%Jx(cea%Msing)
          cea%Jx(cea%Msing) = cea%Jx(cea%Nlm)
          cea%Jx(cea%Nlm) = ja
          aa = cea%Atwt(cea%Msing)
          cea%Atwt(cea%Msing) = cea%Atwt(cea%Nlm)
          cea%Atwt(cea%Nlm) = aa
          aa = cea%B0(cea%Msing)
          cea%B0(cea%Msing) = cea%B0(cea%Nlm)
          cea%B0(cea%Nlm) = aa
          aa = cea%pisave(cea%Msing)
          cea%pisave(cea%Msing) = cea%pisave(cea%Nlm)
          cea%pisave(cea%Nlm) = aa

          do i = 1, 2
             aa = cea%B0p(cea%Msing, i)
             cea%B0p(cea%Msing, i) = cea%B0p(cea%Nlm, i)
             cea%B0p(cea%Nlm, i) = aa
          end do
       end if
    else if (.not. cea%newcom .and. cea%Trace == 0.) then
       go to 600
    end if

    cea%Msing = 0
    tsize = xsize
    go to 500

1400 p%Ttt = cea%Tt
    p%Ppp = cea%Pp
    p%Vlm = R0 * cea%Enn * cea%Tt / cea%Pp
    p%Hsum = p%Hsum * cea%Tt
    p%Wm = 1. / cea%Enn
    gasfrc = cea%Enn / p%Totn

    if (cea%legacy_mode .and. gasfrc < 0.0001) write(IOOUT, '(/" WARNING!  RESULTS MAY BE WRONG FOR POINT", i3, " DUE TO", &
         & /" LOW MOLE FRACTION OF GASES (", e15.8, ") (EQLBRM)")') Npt, gasfrc

    if (cea%Trace /= 0) then
       if (lelim == 0) then
          do concurrent (j = 1:cea%Ng, cea%Enln(j) > -87)
             p%En(j) = exp(cea%Enln(j))
          end do
       else
          do concurrent (j = 1:cea%Ng, all(cea%A(lelim:ls, j) == 0))
             p%En(j) = exp(cea%Enln(j))
          end do
       end if
    end if

    if (cea%legacy_mode .and. p%Debug) write(IOOUT, '(/" POINT=", i3, 3x, "P=", e13.6, 3x, "T=", e13.6, /3x, "H/R=", &
         & e13.6, 3x, "S/R=", e13.6, /3x, "M=", e13.6, 3x, "CP/R=", e13.6, 3x, &
         & "DLVPT=", e13.6, /3x, "DLVTP=", e13.6, 3x, "GAMMA(S)=", e13.6, 3x, "V=", e13.6)') &
         Npt, cea%Pp, cea%Tt, p%Hsum, p%Ssum, p%Wm, p%Cpr, p%Dlvpt, p%Dlvtp, p%Gammas, p%Vlm

    if (cea%Tt >= cea%Tg(1) .and. cea%Tt <= cea%Tg(4)) go to 1600

    if (cea%Shock) go to 1600

    if (cea%legacy_mode) write(IOOUT, '(" THE TEMPERATURE=", e12.4, " IS OUT OF RANGE FOR POINT", i5, "(EQLBRM)")') cea%Tt, Npt

    if (cea%Tt >= cea%Tg(1) * 0.8d0 .and. cea%Tt <= cea%Tg(4) * 1.1d0) go to 1600

    cea%ipt = cea%ipt + 1

1500 cea%Tt = 0
    cea%ipt = cea%ipt - 1
    if (cea%legacy_mode) write(IOOUT, '(/" CALCULATIONS STOPPED AFTER POINT", I3, "(EQLBRM)")') Npt

1600 cea%Lsave = cea%Nlm
    cea%Nlm = ls

    if (cea%Npr > 0) cea%Gonly = .false.

    return
  end subroutine EQLBRM



  subroutine FROZEN(cea)
    !***********************************************************************
    ! CALCULATE PROPERTIES WITH FROZEN COMPOSITION AT ASSIGNED ENTROPY
    ! AND PRESSURE.  CALLED FROM ROCKET.
    !***********************************************************************
    use mod_types
    implicit none

    type(CEA_Core_Problem), intent(inout):: cea

    ! LOCAL VARIABLES
    integer:: iter, j
    real(8):: dlnt, dlpm

    type(CEA_Point), pointer:: p !< current point
    type(CEA_Point), pointer:: pfz

    p => cea%points(cea%iOF, cea%ipt)
    pfz => cea%points(cea%iOF, cea%Nfz)

    cea%Convg = .false.
    cea%Tln = log(cea%Tt)
    dlpm = log(cea%Pp * pfz%Wm)

    do concurrent (j = 1:cea%Ng, pfz%En(j) /= 0)
       cea%Deln(j) = -(log(pfz%En(j)) + dlpm)
    end do

    do iter = 1, 8
       call CPHS(cea)

       cea%Cpsum = sum(pfz%En(1:cea%Ng) * cea%Cp(1:cea%Ng))
       p%Ssum = sum(pfz%En(1:cea%Ng) * (cea%S(1:cea%Ng) + cea%Deln(1:cea%Ng)))

       if (cea%Npr /= 0) then
          cea%Cpsum = cea%Cpsum + sum(pfz%En(cea%Jcond(1:cea%Npr)) * cea%Cp(cea%Jcond(1:cea%Npr)))
          p%Ssum = p%Ssum + sum(pfz%En(cea%Jcond(1:cea%Npr)) * cea%S(cea%Jcond(1:cea%Npr)))
       end if

       if (cea%Convg) then
          p%Hsum = sum(pfz%En(1:cea%Ngc) * cea%H0(1:cea%Ngc)) * cea%Tt

          p%Ttt = cea%Tt
          p%Gammas = cea%Cpsum / (cea%Cpsum - 1 / pfz%Wm)
          p%Vlm = R0 * cea%Tt / (pfz%Wm * cea%Pp)
          p%Wm = pfz%Wm
          p%Dlvpt = -1
          p%Dlvtp = 1
          p%Totn = pfz%Totn
          p%Ppp = cea%Pp
          p%Cpr = cea%Cpsum

          if (cea%Tt >= cea%Tg(1) * 0.8d0) then
             if (all(pfz%En(cea%Ngp1:cea%Ngc) == 0)) return
             if (all(cea%Temp(1, cea%Ngp1-cea%Ng:cea%Ngc-cea%Ng) - 50 <= cea%Tt &
                  .and. cea%Tt <= cea%Temp(2, cea%Ngp1-cea%Ng:cea%Ngc-cea%Ng) + 50)) return
          end if

          cea%Tt = 0
          cea%ipt = cea%ipt - 1
          return

       else
          dlnt = (pfz%Ssum - p%Ssum) / cea%Cpsum
          cea%Tln = cea%Tln + dlnt
          if (abs(dlnt) < 0.5d-4) cea%Convg = .true.
          cea%Tt = exp(cea%Tln)
       end if
    end do

    if (cea%legacy_mode) write(IOOUT, '(/" FROZEN DID NOT CONVERGE IN 8 ITERATIONS (FROZEN)")')

    cea%Tt = 0
    cea%ipt = cea%Nfz - 1

    return
  end subroutine FROZEN


  subroutine HCALC(cea)
    !***********************************************************************
    ! CALCULATE PROPERTIES FOR TOTAL REACTANT USING THERMO DATA FOR
    ! ONE OR MORE REACTANTS. USED ONLY FOR SHOCK AND DETON PROBLEMS.
    !***********************************************************************
    use mod_types
    use mod_functions
    implicit none

    type(CEA_Core_Problem), intent(inout):: cea

    ! LOCAL VARIABLES
    character(6):: date(maxNgc)
    character(2):: el(5)
    character(15):: sub
    integer:: i, icf, ifaz, itot, j, k, l, m, n, nall, nint, ntgas, ntot
    real(8):: bb(5), enj, er, sj, t1, t2, tem, thermo(9, 3), Tsave
    integer:: io_thermo

    type(CEA_Point), pointer:: p !< current point
    p => cea%points(cea%iOF, cea%ipt)

    Tsave = cea%Tt

    cea%Tm = 0
    if (cea%Pp > 0) cea%Tm = log(cea%Pp * cea%Wmix)

    p%Ssum = 0
    cea%Hpp(1) = 0
    cea%Hpp(2) = 0
    cea%Hsub0 = 0
    cea%Cpmix = 0
    tem = (1 + cea%Oxfl)

    ! LOOP ON REACTANTS.
    ! if oxidant, k = 1
    ! if fuel,    k = 2
    cea%Nspr = cea%Nspx
    do n = 1, cea%Nreac
       k = 2
       if (cea%Fox(n)(:1) == 'O' .or. cea%Fox(n)(:1) == 'o') k = 1
       if (cea%Tt == 0) cea%Tt = cea%Rtemp(n)

       j = cea%Jray(n)

       if (j == 0) then
          ! SEARCH FOR REACTANT IN STORED THERMO SPECIES. STORE INDEX IN JRAY(N).
          ifaz = 0
          do j = 1, cea%Ngc
             if (cea%Rname(n) == cea%Prod(j) .or. '*' // cea%Rname(n) == cea%Prod(j)) then
                cea%Jray(n) = j
                if (j > cea%Ng) then
                   if (cea%legacy_mode) write(IOOUT, '(/" REACTANTS MUST BE GASEOUS FOR THIS PROBLEM (HCALC)")')
                   cea%Tt = 0
                   cea%Cpmix = 0
                   return
                end if
                go to 50
             end if
          end do

          ! SEARCH THERMO.LIB FOR SPECIES.
          open(newunit = io_thermo, file = cea%filename_thermo_lib, status = 'old', form = 'unformatted', action = 'read')

          read(io_thermo) cea%Tg, ntgas, ntot, nall

          cea%Nspr = cea%Nspr + 1
          do itot = 1, nall
             if (itot <= ntot) then
                icf = 3
                if (itot > ntgas) icf = 1
                read(io_thermo) sub, nint, date(cea%Nspr), (el(j), bb(j), j = 1, 5), ifaz, &
                     T1, T2, cea%Mw(cea%Nspr), ((thermo(l, m), l = 1, 9), m = 1, icf)
             else
                read(io_thermo) sub, nint, date(cea%Nspr), (el(j), bb(j), j = 1, 5), ifaz, T1, T2, cea%Mw(cea%Nspr), er
                if (nint /= 0) then
                   read(io_thermo) ((thermo(i, j), i = 1, 9), j = 1, nint)
                   icf = nint
                end if
             end if

             if (sub == cea%Rname(n) .or. sub == '*' // cea%Rname(n)) then
                if (ifaz <= 0 .and. nint > 0) then
                   do j = 1, 5
                      if (bb(j) == 0) exit
                      cea%Nfla(n) = j
                      cea%Ratom(n, j) = el(j)
                      cea%Rnum(n, j) = bb(j)
                   end do

                   cea%Jray(n) = cea%Nspr
                   j = cea%Nspr

                   do concurrent (l = 1:icf, m = 1:9)
                      cea%Coef(m, j, l) = thermo(m, l)
                   end do

                   go to 50

                else
                   if (cea%legacy_mode) then
                      if (ifaz > 0) write(IOOUT, '(/" REACTANTS MUST BE GASEOUS FOR THIS PROBLEM (HCALC)")')
                      if (nint == 0) write(IOOUT, '(/" COEFFICIENTS FOR ", a15, " ARE NOT AVAILABLE (HCALC)")') cea%Rname(n)
                   end if
                   cea%Tt = 0
                   cea%Cpmix = 0
                   close(io_thermo)
                   return
                end if
             end if
          end do

          close(io_thermo)

          cea%Nspr = cea%Nspr - 1
          if (cea%legacy_mode) write(IOOUT, '(/" ERROR IN DATA FOR ", a15, " CHECK NAME AND TEMPERATURE", &
               & " RANGE IN", /, " thermo.inp (HCALC)")') cea%Rname(n)

          cea%Energy(n) = ' '
          cea%Tt = 0
          cea%Cpmix = 0

          return
       end if

       ! CALCULATE EN FOR REACTANT AND CALCULATE PROPERTIES.
50     if (cea%Moles) then
          enj = cea%Pecwt(n) / cea%Wp(k)
       else
          enj = cea%Pecwt(n) / cea%Rmw(n)
       end if
       enj = enj / tem
       if (k == 1) enj = enj * cea%Oxfl

       cea%Tln = log(cea%Tt)
       p%En(j) = enj
       l = 1

       if (ifaz <= 0) then
          if (cea%Tt > cea%Tg(2)) l = 2
          if (cea%Tt > cea%Tg(3) .and. ifaz < 0) l = 3
       end if

       cea%S(j) = calc_entropy(cea%Coef(:, j, l), cea%Tt, cea%Tln, .false.)

       cea%H0(j) = calc_enthalpy(cea%Coef(:, j, l), cea%Tt, cea%Tln, .false.)

       cea%Cp(j) = calc_specific_heat(cea%Coef(:, j, l), cea%Tt)

       if (abs(cea%H0(j)) < 0.01) cea%H0(j) = 0

       ! ADD CONTRIBUTION TO CP, H, AND S OF TOTAL REACTANT.
       cea%Cpmix = cea%Cpmix + cea%Cp(j) * enj

       ! FOR CONDENSED SPECIES:  SJ = S(J)
       sj = cea%S(j) - log(enj) - cea%Tm
       p%Ssum = p%Ssum + enj * sj
       er = cea%H0(j) * enj * cea%Tt
       cea%Hsub0 = cea%Hsub0 + er
       cea%Hpp(k) = cea%Hpp(k) + er
    end do

    if (Tsave /= 0) cea%Tt = Tsave

    return
  end subroutine HCALC


  subroutine MATRIX(cea)
    !***********************************************************************
    ! SET UP ITERATION OR DERIVATIVE MATRIX.
    !***********************************************************************
    use mod_types
    implicit none

    type(CEA_Core_Problem), intent(inout):: cea

    ! LOCAL VARIABLES
    integer:: i, iq, iq2, iq3, isym, j, k, kk, kmat
    real(8):: energyl, f, h, ss, sss, term, term1

    type(CEA_Point), pointer:: p !< current point

    p => cea%points(cea%iOF, cea%ipt)

    iq = cea%Nlm + cea%Npr
    cea%Iq1 = iq + 1
    iq2 = cea%Iq1 + 1
    iq3 = iq2 + 1
    kmat = iq3
    if (.not. cea%Convg .and. cea%Tp) kmat = iq2
    cea%Imat = kmat - 1
    cea%G = 0
    sss = 0
    p%Hsum = 0

    ! BEGIN SET-UP OF ITERATION OR DERIVATIVE MATRIX
    do j = 1, cea%Ng
       cea%Mu(j) = cea%H0(j) - cea%S(j) + cea%Enln(j) + cea%Tm
       if (p%En(j) /= 0) then
          h = cea%H0(j) * p%En(j)
          f = cea%Mu(j) * p%En(j)
          ss = h - f
          term1 = h
          if (kmat == iq2) term1 = f

          do i = 1, cea%Nlm
             if (cea%A(i, j) /= 0) then
                term = cea%A(i, j) * p%En(j)
                do concurrent (k = i:cea%Nlm)
                   cea%G(i, k) = cea%G(i, k) + cea%A(k, j) * term
                end do
                cea%G(i, cea%Iq1) = cea%G(i, cea%Iq1) + term
                cea%G(i, iq2) = cea%G(i, iq2) + cea%A(i, j) * term1
                if (.not. (cea%Convg .or. cea%Tp)) then
                   cea%G(i, iq3) = cea%G(i, iq3) + cea%A(i, j) * f
                   if (cea%Sp) cea%G(iq2, i) = cea%G(iq2, i) + cea%A(i, j) * ss
                end if
             end if
          end do

          if (kmat /= iq2) then
             if (cea%Convg .or. cea%Hp) then
                cea%G(iq2, iq2) = cea%G(iq2, iq2) + cea%H0(j) * h
                if (.not. cea%Convg) then
                   cea%G(iq2, iq3) = cea%G(iq2, iq3) + cea%H0(j) * f
                   cea%G(cea%Iq1, iq3) = cea%G(cea%Iq1, iq3) + f
                end if
             else
                cea%G(iq2, cea%Iq1) = cea%G(iq2, cea%Iq1) + ss
                cea%G(iq2, iq2) = cea%G(iq2, iq2) + cea%H0(j) * ss
                cea%G(iq2, iq3) = cea%G(iq2, iq3) + cea%Mu(j) * ss
                cea%G(cea%Iq1, iq3) = cea%G(cea%Iq1, iq3) + f
             end if
          end if
          cea%G(cea%Iq1, iq2) = cea%G(cea%Iq1, iq2) + term1
       end if
    end do

    ! CONDENSED SPECIES
    if (cea%Npr /= 0) then
       do k = 1, cea%Npr
          j = cea%Jcond(k)
          kk = cea%Nlm + k
          cea%Mu(j) = cea%H0(j) - cea%S(j)

          do concurrent (i = 1:cea%Nlm)
             cea%G(i, kk) = cea%A(i, j)
             cea%G(i, kmat) = cea%G(i, kmat) - cea%A(i, j) * p%En(j)
          end do

          cea%G(kk, iq2) = cea%H0(j)
          cea%G(kk, kmat) = cea%Mu(j)
          p%Hsum = p%Hsum + cea%H0(j) * p%En(j)

          if (cea%Sp) then
             sss = sss + cea%S(j) * p%En(j)
             cea%G(iq2, kk) = cea%S(j)
          end if
       end do
    end if

    sss = sss + cea%G(iq2, cea%Iq1)
    p%Hsum = p%Hsum + cea%G(cea%Iq1, iq2)
    cea%G(cea%Iq1, cea%Iq1) = cea%Sumn - cea%Enn

    ! REFLECT SYMMETRIC PORTIONS OF THE MATRIX
    isym = cea%Iq1
    if (cea%Hp .or. cea%Convg) isym = iq2

    do i = 1, isym
       !DIR$ IVDEP
       do concurrent (j = i:isym)
          cea%G(j, i) = cea%G(i, j)
       end do
    end do

    ! COMPLETE THE RIGHT HAND SIDE
    if (.not. cea%Convg) then
       do concurrent (i = 1:cea%Nlm)
          cea%G(i, kmat) = cea%G(i, kmat) + cea%B0(i) - cea%G(i, cea%Iq1)
       end do
       cea%G(cea%Iq1, kmat) = cea%G(cea%Iq1, kmat) + cea%Enn - cea%Sumn

       ! COMPLETE ENERGY ROW AND TEMPERATURE COLUMN
       if (kmat /= iq2) then
          if (cea%Sp) energyl = cea%S0 + cea%Enn - cea%Sumn - sss
          if (cea%Hp) energyl = cea%Hsub0/cea%Tt - p%Hsum
          cea%G(iq2, iq3) = cea%G(iq2, iq3) + energyl
          cea%G(iq2, iq2) = cea%G(iq2, iq2) + cea%Cpsum
       end if

    else
       if (cea%Pderiv) then
          ! PDERIV = .true.-- SET UP MATRIX TO SOLVE FOR DLVPT
          cea%G(cea%Iq1, iq2) = cea%Enn
          do concurrent (i = 1:iq)
             cea%G(i, iq2) = cea%G(i, cea%Iq1)
          end do
       end if
       cea%G(iq2, iq2) = cea%G(iq2, iq2) + cea%Cpsum
    end if

    if (cea%Vol .and. .not. cea%Convg) then
       ! CONSTANT VOLUME MATRIX
       if (kmat == iq2) then
          do concurrent (i = 1:iq)
             cea%G(i, cea%Iq1) = cea%G(i, iq2)
          end do
       else

          !DIR$ IVDEP
          do concurrent (i = 1:iq)
             cea%G(cea%Iq1, i) = cea%G(iq2, i) - cea%G(cea%Iq1, i)
             cea%G(i, cea%Iq1) = cea%G(i, iq2) - cea%G(i, cea%Iq1)
             cea%G(i, iq2) = cea%G(i, iq3)
          end do

          cea%G(cea%Iq1, cea%Iq1) = cea%G(iq2, iq2) - cea%G(cea%Iq1, iq2) - cea%G(iq2, cea%Iq1)
          cea%G(cea%Iq1, iq2) = cea%G(iq2, iq3) - cea%G(cea%Iq1, iq3)
          if (cea%Hp) cea%G(cea%Iq1, iq2) = cea%G(cea%Iq1, iq2) + cea%Enn
       end if

       kmat = cea%Imat
       cea%Imat = cea%Imat - 1
    end if

    return
  end subroutine MATRIX



  subroutine NEWOF(cea)
    !***********************************************************************
    ! CALCULATE NEW VALUES OF B0 AND HSUB0 FOR NEW OF RATIO
    !***********************************************************************
    use mod_types
    implicit none

    type(CEA_Core_Problem), intent(inout):: cea

    ! LOCAL VARIABLES
    integer:: i, j
    real(8):: assval, bigb, bratio, dbi, smalb, tem, v1, v2

    type(CEA_Point), pointer:: p !< current point

    cea%iOF = cea%iOF + 1
    cea%ipt = 1

    p => cea%points(cea%iOF, cea%ipt)

    if (cea%legacy_mode .and. .not. cea%Short) write(IOOUT, '(/" O/F = ", f10.6)') cea%Oxfl
    cea%Eqrat = 0
    tem = cea%Oxfl + 1
    v2 = (cea%Oxfl * cea%Vmin(1) + cea%Vmin(2)) / tem
    v1 = (cea%Oxfl * cea%Vpls(1) + cea%Vpls(2)) / tem
    if (v2 /= 0.) cea%Eqrat = abs(v1/v2)
    do i = 1, cea%Nlm
       cea%B0(i) = (cea%Oxfl * cea%B0p(i, 1) + cea%B0p(i, 2)) / tem
       dbi = abs(cea%B0(i))
       if (i == 1) then
          bigb = dbi
          smalb = dbi
       else if (dbi /= 0.) then
          if (dbi < smalb) smalb = dbi
          if (dbi > bigb) bigb = dbi
       end if
    end do
    cea%Bcheck = bigb * .000001d0

    ! CALCUALTE MOLECULAR WEIGHT OF TOTAL REACTANT, WMIX.
    if (all(cea%Am /= 0)) then
       cea%Wmix = (cea%Oxfl + 1) * cea%Am(1) * cea%Am(2) / (cea%Am(1) + cea%Oxfl * cea%Am(2))
    else
       cea%Wmix = cea%Am(2)
       if (cea%Am(2) == 0.0) cea%Wmix = cea%Am(1)
    end if

    ! IF ASSIGNED U OR H NOT GIVEN IN PROB DATA, INITIAL HSUB0 = 1.d30
    if (cea%Size == 0) assval = cea%Hsub0
    if (assval >= 1.d30) cea%Hsub0 = (cea%Oxfl * cea%Hpp(1) + cea%Hpp(2)) / tem

    ! NOTE THAT "BRATIO" IS "BRATIO" IN SEC 3.2 IN RP-1311.
    bratio = smalb / bigb
    cea%Size = 18.420681d0
    if (bratio < 1.d-5) cea%Size = log(1000/bratio)
    cea%Jsol = 0
    cea%Jliq = 0

    if (cea%legacy_mode .and. .not. cea%Short) then
       write(IOOUT, '(/, 23x, "EFFECTIVE FUEL", 5x, "EFFECTIVE OXIDANT", 8x, "MIXTURE")')
       if (cea%Vol) write(IOOUT, '(" INTERNAL ENERGY", 11x, "u(2)/R", 14x, "u(1)/R", 14x, "u0/R")')
       if (.not. cea%Vol) write(IOOUT, '(" ENTHALPY", 18x, "h(2)/R", 14x, "h(1)/R", 15x, "h0/R")')
       write(IOOUT, '(" (KG-MOL)(K)/KG", 4x, e18.8, 2e20.8)') cea%Hpp(2), cea%Hpp(1), cea%Hsub0
       write(IOOUT, '(/" KG-FORM.WT./KG", 13x, "bi(2)", 15x, "bi(1)", 15x, "b0i")')
    end if

    do i = 1, cea%Nlm
       j = cea%Jcm(i)
       if (cea%legacy_mode .and. .not. cea%Short) &
            write(IOOUT, '(1x, a16, 3e20.8)') cea%Prod(j), cea%B0p(i, 2), cea%B0p(i, 1), cea%B0(i)
    end do

    return
  end subroutine NEWOF


  subroutine ROCKET(cea)
    !***********************************************************************
    ! EXECUTIVE ROUTINE FOR ROCKET PROBLEMS.
    !***********************************************************************
    use mod_types
    use mod_legacy_io, only: RKTOUT
    implicit none

    type(CEA_Core_Problem), intent(inout):: cea

    ! LOCAL VARIABLES
    integer:: i, i01, i12, iof, iplt1, iplte, ipp, isub, isup1, isupsv, itnum, &
         itrot, nipp, niter, nn, npr1, nptth
    logical:: done, seql, thi
    real(8):: a1l = -1.26505, b1 = 1.0257, c1 = -1.2318, pa = 1e5
    real(8):: acatsv, aeatl, appl, aratio, asq, check, cprf, dd, dh, &
         dlnp, dlnpe, dlt, dp, eln, mat, msq, Pp_old, pcpa, pcplt, pinf, pinj, &
         pinjas, pjrat, ppa, pr, pracat, prat, pratsv, pvg, test, tmelt, usq
    integer:: ip, it

    integer:: iar
    type(CEA_Point), pointer:: p !< current point
    type(CEA_Point), pointer:: p1, p2, p4, p12
    type(CEA_Point), pointer:: pfz
    type(CEA_Point), pointer:: pth

    iplte = cea%Iplt
    isup1 = 1

    do concurrent (iof = 1:cea%Nof)
       cea%points(iof, 1)%App = 1
    end do

    cea%Iopt = 0
    cea%Npp = cea%Npp + 2
    nn = cea%Npp
    i01 = 0
    i12 = 1
    nipp = 1
    nptth = 2
    if (cea%Fac) then
       cea%Eql = .true.
       cea%Npp = cea%Npp + 1
       if (cea%AcAt /= 0) then
          cea%Iopt = 1
       else if (cea%mdotByAc /= 0) then
          cea%Iopt = 2
       else
          if (cea%legacy_mode) write(IOOUT, '(/" FATAL ERROR!! EITHER mdot OR ac/at MISSING FOR fac PROBLEM (ROCKET)")')
          cea%Tt = 0
          cea%Iplt = max(cea%Iplt, iplte)
          return
       end if
       i01 = 1
       i12 = 2
       nipp = 2
       nptth = 3
       do i = cea%Nsub, 1, - 1
          cea%Subar(i+1) = cea%Subar(i)
       end do
       cea%Nsub = cea%Nsub + 1
       if (cea%Iopt /= 1) then
          if (cea%AcAt == 0) cea%AcAt = 2
       end if
       cea%Subar(1) = cea%AcAt
    else if (.not. cea%Eql .and. cea%Nfz > 1 .and. cea%Nsub > 0) then
       cea%Nsub = 0
       if (cea%legacy_mode) write(IOOUT, '(/" WARNING!!  FOR FROZEN PERFORMANCE, SUBSONIC AREA ", /, &
            & " RATIOS WERE OMITTED SINCE nfz IS GREATER THAN 1 (ROCKET)")')
    end if
    nn = nn + cea%Nsub + cea%Nsup
    if (cea%Nfz > 2 .and. nn > Ncol-2) then
       if (cea%legacy_mode) write(IOOUT, '(/" WARNING!!  nfz NOT ALLOWED TO BE > 2 IF THE TOTAL", /, &
            & " NUMBER OF POINTS IS >", i3, " (ROCKET)")') Ncol - 2
       cea%Nfz = 1
       cea%Froz = .false.
    end if
    seql = cea%Eql
    cea%Tt = cea%Tcest
    cea%Pp = cea%P(1)

    do concurrent (iof = 1:cea%Nof)
       cea%points(iof, i12)%App = 1
    end do

    ! LOOP FOR EACH O/F
    iofLoop: do iof = 1, cea%Nof
       cea%Oxfl = cea%Oxf(iof)
       if (cea%T(1) /= 0.) then
          cea%Tp = .true.
       else
          cea%Hp = .true.
       end if
       cea%Sp = .false.
       call NEWOF(cea)

       cea%ipt = 1

       p1 => cea%points(cea%iOF, 1)
       p2 => cea%points(cea%iOF, 2)
       p4 => cea%points(cea%iOF, 4)
       p12 => cea%points(cea%iOF, i12)
       pfz => cea%points(cea%iOF, cea%Nfz)
       pth => cea%points(cea%iOF, nptth)

       do it = 1, cea%Nt
          cea%Tt = cea%T(it)

          ! LOOP FOR CHAMBER PRESSURES
          do ip = 1, cea%Np
             itnum = 0
             cea%Area = .false.
             if (cea%T(1) == 0) cea%Hp = .true.
             if (cea%T(1) /= 0) cea%Tp = .true.
             cea%Sp = .false.
             cea%Eql = .true.
             isub = 1
             cea%Isup = 1
             cea%Pp = cea%P(ip)
             pinf = cea%Pp
             ipp = 1
             itrot = 3
             isupsv = 1
             niter = 1
             cea%Page1 = .true.
             iplt1 = iplte
             cea%Iplt = iplte
             done = .false.

             ! LOOP FOR OUTPUT COLUMNS
250          iar = cea%ipt
             p => cea%points(cea%iOF, cea%ipt)
             if (cea%Eql) then
                call EQLBRM(cea)
                if (cea%ipt == cea%Nfz) cprf = cea%Cpsum
             else
                call FROZEN(cea)
             end if
             ! TT = 0 IF NO CONVERGENCE
             if (cea%Tt /= 0.) then
                ! TEST FOR FINITE AREA COMBUSTOR
                if (.not. cea%Fac) go to 400
                pinjas = cea%P(ip) * pa
                pinj = pinjas
                if (cea%ipt == 1 .and. cea%Trnspt) call TRANP(cea)
                if (cea%ipt == 2) pinf = p%Ppp
                if (cea%ipt /= 1) go to 400

                ! INITIAL ESTIMATE FOR PC (AND ACAT IF NOT ASSIGNED)
                do i = 1, 4
                   prat = (b1 + c1 * cea%AcAt) / (1 + a1l * cea%AcAt)
                   ppa = pinj * prat
                   if (cea%Iopt == 1) exit
                   cea%AcAt = ppa / (cea%mdotByAc * 2350)
                   if (cea%AcAt >= 1) then
                      pratsv = prat
                      if (cea%legacy_mode .and. cea%Debugf) then
                         if (i <= 1) write(IOOUT, '(/"  ITERATION", 9x, "PC", 7x, "CONTRACTION RATIO")')
                         write(IOOUT, '(5x, i2, 7x, f12.2, 3x, f12.6)') i, ppa, cea%AcAt
                      end if
                   else
                      if (cea%legacy_mode) write(IOOUT, '(/" INPUT VALUE OF mdot/a =", f12.3, " IS TOO LARGE."/ &
                           & " GIVES CONTRACTION RATIO ESTIMATE LESS THAN 1 (ROCKET)")') cea%mdotByAc
                      cea%Tt = 0
                      exit iofLoop
                   end if
                end do

                if (i > 4) cea%Subar(1) = cea%AcAt

                cea%Pp = ppa / pa
                p1%App = cea%Pp / p1%Ppp
                go to 1100

             else
                if (cea%ipt < 1) exit iofLoop
                if (.not. cea%Area) go to 600
                cea%ipt = iar - 1
                cea%Isup = cea%Nsup + 2
                cea%Isv = 0
                itnum = 0
                go to 950
             end if

             ! INITIALIZE FOR THROAT
400          if (ipp > nipp) then
                p => cea%points(cea%iOF, cea%ipt)
                usq = 2 * (cea%points(cea%iOF, 1)%Hsum - p%Hsum) * R0
                if (ipp > nptth) go to 600
                ! THROAT
                if (.not. thi) then
                   cea%Vv = pth%Vlm
                   pvg = cea%Pp * cea%Vv * pth%Gammas
                   if (pvg == 0) then
                      if (cea%legacy_mode) write(IOOUT, '(/" WARNING!!  DIFFICULTY IN LOCATING THROAT (ROCKET)")')
                      go to 550
                   else
                      msq = usq / pvg
                      if (cea%legacy_mode .and. p1%Debug .or. p2%Debug) &
                           write(IOOUT, '(/" USQ=", e15.8, 5x, "PVG=", e15.8)') usq, pvg
                      dh = abs(msq - 1)
                      if (dh <= 0.4d-4) go to 550
                      if (itrot > 0) then
                         Pp_old = cea%Pp
                         if (cea%Jsol /= 0) then
                            tmelt = cea%Tt
                            cea%Pp = cea%Pp * (1 + msq * pth%Gammas) / (pth%Gammas + 1)
                         else if (tmelt == 0) then
                            cea%Pp = cea%Pp * (1 + msq * pth%Gammas) / (pth%Gammas + 1)
                         else
                            if (cea%legacy_mode) write(IOOUT, '(/" WARNING!!  DISCONTINUITY AT THE THROAT (ROCKET)")')
                            dlt = log(tmelt / cea%Tt)
                            dd = dlt * pth%Cpr / (cea%Enn * pth%Dlvtp)
                            cea%Pp = cea%Pp * EXP(dd)
                            pth%App = cea%P(ip) / cea%Pp
                            if (cea%Fac) pth%App = pinf / cea%Pp
                            if (cea%legacy_mode .and. cea%Eql .and. .not. cea%Short) write(IOOUT, '(" Pinf/Pt =", F9.6)') pth%App
                            thi = .true.
                            go to 250
                         end if
                         go to 500
                      else if (itrot < 0) then
                         if (itrot < -19) then
                            if (cea%legacy_mode) write(IOOUT, '(/" WARNING!!  DIFFICULTY IN LOCATING THROAT (ROCKET)")')
                            go to 550
                         else
                            if (cea%Npr /= npr1) go to 550
                            cea%Pp = cea%Pp - dp
                            go to 500
                         end if
                      else if (cea%Npr == npr1) then
                         if (cea%legacy_mode) write(IOOUT, '(/" WARNING!!  DIFFICULTY IN LOCATING THROAT (ROCKET)")')
                         go to 550
                      else
                         dp = abs(cea%Pp - Pp_old) / 20
                         cea%Pp = max(cea%Pp, Pp_old)
                         if (cea%legacy_mode) write(IOOUT, '(/" WARNING!!  DISCONTINUITY AT THE THROAT (ROCKET)")')
                         cea%Pp = cea%Pp - dp
                         go to 500
                      end if
                   end if
                else
                   pth%Gammas = 0
                   go to 550
                end if
             else
                if (.not. cea%Fac .and. cea%Trnspt) call TRANP(cea)
                if (cea%ipt == cea%Nfz) cea%Eql = seql
                cea%Tp = .false.
                cea%Hp = .false.
                cea%Sp = .true.
                cea%S0 = p12%Ssum
             end if

450          tmelt = 0
             itrot = 3
             thi = .false.
             pth%App = ((p12%Gammas + 1) / 2)**(p12%Gammas / (p12%Gammas - 1))
             if (cea%legacy_mode .and. cea%Eql .and. .not. cea%Short) write(IOOUT, '(" Pinf/Pt =", f9.6)') pth%App
             cea%Pp = Pinf / pth%App
             cea%Isv = -i12
             go to 1200

500          npr1 = cea%Npr
             pth%App = cea%P(ip) / cea%Pp
             if (cea%Fac) pth%App = Pinf / cea%Pp
             if (cea%legacy_mode .and. cea%Eql .and. .not. cea%Short) write(IOOUT, '(" Pinf/Pt =", f9.6)') pth%App
             itrot = itrot - 1
             go to 250

550          cea%Awt = cea%Enn * cea%Tt / (cea%Pp * sqrt(usq))
             pcplt = log(pth%App)

600          cea%Isv = 0
             p => cea%points(cea%iOF, cea%ipt)
             p%AeAt = cea%Enn * p%Ttt / (cea%Pp * sqrt(usq) * cea%Awt)
             if (cea%Tt == 0) go to 1150
             if (cea%Area) go to 750
             if (cea%Trnspt .and. (.not. cea%Fac .or. done .or. cea%ipt > 2)) call TRANP(cea)
             if (cea%ipt == cea%Nfz) cea%Eql = seql
             if (cea%Fac) then
                if (cea%ipt == nptth) then
                   cea%Area = .true.
                   go to 750
                else if (cea%ipt == 2 .and. done) then
                   cea%ipt = 3
                   !  The following statement was corrected 1/30/2004.  Only fac parameters 
                   !    after combustion were affected--generally extra or missing points.
                   !  (remove) if (ipp <= cea%Npp) ipp = ipp - 1
                   if (ipp < cea%Npp .or. cea%Npp == 4) ipp = ipp - 1
                end if
             end if

650          if (ipp < cea%Npp) go to 1100

700          if (cea%Nsub == i01 .and. cea%Nsup == 0) go to 1150
             cea%Area = .true.

             ! PCP ESTIMATES FOR AREA RATIOS
750          p => cea%points(cea%iOF, cea%ipt)
             if (itnum == 0) then
                dlnp = 1
                itnum = 1
                aratio = cea%Subar(isub)
                if ((.not. cea%Fac .or. done) .and. cea%Nsub <= i01) aratio = cea%Supar(cea%Isup)
                if (.not. cea%Eql .and. cea%Nfz >= 3) then
                   if (aratio <= pfz%Aeat) then
                      if (cea%legacy_mode) write(IOOUT, '(/, " WARNING!! FOR FROZEN PERFORMANCE, POINTS WERE OMITTED", &
                           & " WHERE THE ASSIGNED", /, " SUPERSONIC AREA RATIOS WERE ", &
                           & "LESS THAN THE VALUE AT POINT nfz =", i3, " (ROCKET)")') cea%Nfz
                      go to 1050
                   end if
                end if
                if (aratio  <  1) then
                   if (cea%legacy_mode) write(IOOUT, '(/" AN ASSIGNED AREA RATIO IS < 1 (ROCKET)")')
                   go to 1050
                end if
                eln = log(aratio)

                if ((done .or. .not. cea%Fac) .and. cea%Nsub <= i01) then
                   if (cea%Nfz == ipp) isupsv = cea%Isup
                   if (cea%Supar(cea%Isup) < 2) then
                      appl = sqrt(eln*(1.535d0 + 3.294d0 * eln)) + pcplt
                   else
                      if (cea%Isup > isup1) then
                         if (cea%Supar(cea%Isup-1) >= 2) go to 850
                      end if
                      appl = pth%Gammas + eln * 1.4
                   end if
                   go to 1100
                end if

                appl = pcplt / (cea%Subar(isub) + (10.587 * eln**2 + 9.454) * eln)
                if (Aratio < 1.09) appl = 0.9 * appl
                if (Aratio > 10) appl = appl / Aratio
                if (isub > 1 .or. mod(cea%ipt, Ncol) == 0) go to 1100
                go to 1200

                ! TEST FOR CONVERGENCE ON AREA RATIO.
             else if (p%Gammas > 0) then
                check = 0.00004
                p => cea%points(cea%iOF, cea%ipt)
                if (cea%legacy_mode .and. p%Debug) &
                     write(IOOUT, '(/" ITER=", i2, 2x, "ASSIGNED AE/AT=", f14.7, 3x, "AE/AT=", f14.7, &
                     & /, 2x, "PC/P=", f14.7, 2x, "DELTA LN PCP=", f14.7)') &
                     itnum, aratio, p%AeAt, p%App, dlnp
                if (abs(p%AeAt - Aratio) / Aratio <= check) go to 900
                if (abs(dlnp) < 0.00004) go to 900
                aeatl = log(p%AeAt)
                itnum = itnum + 1
                if (itnum > 10) then
                   if (cea%legacy_mode) &
                        write(IOOUT, '(/" WARNING!!  DID NOT CONVERGE FOR AREA RATIO =", F10.5, " (ROCKET)")') aratio
                   go to 900
                else
                   p => cea%points(cea%iOF, cea%ipt)
                   ! IMPROVED PCP ESTIMATES.
                   asq = p%Gammas * cea%Enn * R0 * cea%Tt
                   dlnpe = p%Gammas * usq / (usq - asq)
                   go to 850
                end if
             else
                if (cea%legacy_mode) write(IOOUT, '(/" WARNING!!  AREA RATIO CALCULATION CANNOT BE DONE ", &
                     & "BECAUSE GAMMAs", /, " CALCULATION IMPOSSIBLE. (ROCKET)")')
                cea%ipt = cea%ipt - 1
                if (cea%Nsub <= 0) isup1 = 100
                if (cea%Nsub < 0.) cea%Nsup = cea%Isup - 1
                if (cea%Nsub > 0) cea%Nsub = isub - 1
                go to 1000
             end if

850          dlnp = dlnpe * eln - dlnpe * aeatl
             appl = appl + dlnp
             if (itnum == 1) go to 1100
             if (appl < 0.) appl = 0.000001
             p => cea%points(cea%iOF, cea%ipt)
             p%App = EXP(appl)
             cea%Pp = Pinf / p%App
             go to 250

             ! CONVERGENCE HAS BEEN REACHED FOR ASSIGNED AREA RATIO
900          p = cea%points(cea%iOF, cea%ipt)
             p%AeAt = Aratio
             if (cea%Fac .and. .not. done) then
                if (cea%Iopt == 1) then
                   ! OPTION 1 FOR FINITE AREA COMBUSTOR. INPUT IS ASSIGNED INJECTOR
                   ! PRESSURE AND CONTRACTION RATIO. IMPROVED ESTIMATE FOR PC
                   cea%Area = .false.
                   itnum = 0
                   ppa = p%Ppp * pa
                   pinj = ppa + 1.d05 * usq / p%Vlm
                   test = (pinj - pinjas) / pinjas
                   pcpa = pinf * pa
                   if (cea%legacy_mode .and. cea%Debugf) then
                      write(IOOUT, '(" ITER", 3x, "TEST", 3x, "ASSIGNED PINJ", 1x, "CALC PINJ", 5x, &
                           & "PC", 7x, "P AT ACAT", 3x, "PREV ACAT", 2x, "ACAT")')
                      write(IOOUT, '(i3, f10.6, 1x, 4f12.2, 2f9.5)') niter, test, pinjas, pinj, pcpa, ppa, acatsv, cea%AcAt
                   end if
                else if (cea%Iopt == 2) then
                   ! OPTION 2 FOR FINITE AREA COMBUSTOR. INPUT IS ASSIGNED INJECTOR
                   ! PRESSURE AND MASS FLOW PER UNIT AREA. IMPROVED ESTIMATE FOR PC
                   ! AND ACAT
                   acatsv = cea%AcAt
                   pratsv = prat
                   cea%Area = .false.
                   itnum = 0
                   ppa = p4%Ppp * pa
                   pinj = ppa + 1.d05 * usq / p4%Vlm
                   mat = pa / (cea%Awt * R0)
                   cea%AcAt = mat / cea%mdotByAc
                   prat = (b1 + c1 * cea%AcAt) / (1 + a1l * cea%AcAt)
                   test = (pinj - pinjas) / pinjas
                   pcpa = pinf * pa
                   if (cea%legacy_mode .and. cea%Debugf) then
                      write(IOOUT, '(" ITER", 3x, "TEST", 3x, "ASSIGNED PINJ", 1x, "CALC PINJ", 5x, &
                           & "PC", 7x, "P AT ACAT", 3x, "PREV ACAT", 2x, "ACAT")')
                      write(IOOUT, '(i3, f10.6, 1x, 4f12.2, 2f9.5)') niter, test, pinjas, pinj, pcpa, ppa, acatsv, cea%AcAt
                   end if
                end if

                if (abs(test) < 0.00002) then
                   done = .true.
                   p1%App = p2%Ppp / p1%Ppp
                   cea%Area = .false.
                   if (cea%Nsub > 1) isub = 2
                   cea%Isv = 4
                   cea%ipt = 2
                   ipp = min(4, cea%Npp)
                   call SETEN(cea)
                   p2%Cpr = p4%Cpr
                   p2%Dlvpt = p4%Dlvpt
                   p2%Dlvtp = p4%Dlvtp
                   p2%Gammas = p4%Gammas
                   p2%Hsum = p4%Hsum
                   p2%Ppp = p4%Ppp
                   p2%App = p1%Ppp/pinf
                   p2%Ssum = p4%Ssum
                   p2%Totn = p4%Totn
                   p2%Ttt = p4%Ttt
                   p2%Vlm = p4%Vlm
                   p2%Wm = p4%Wm
                   if (cea%legacy_mode .and. .not. cea%Short) write(IOOUT, '(" END OF CHAMBER ITERATIONS")')
                   go to 600
                end if

                if (cea%Iopt == 1) then
                   prat = pinjas / pinj
                   cea%Pp = pinf * prat
                else if (cea%Iopt == 2) then
                   pjrat = pinj / pinjas
                   cea%Pp = pinf
                   do i = 1, 2
                      pracat = pratsv / prat
                      pr = pjrat * pracat
                      cea%Pp = cea%Pp / pr
                      pcpa = cea%Pp * pa
                      cea%AcAt = cea%AcAt / pr
                      cea%Subar(1) = cea%AcAt
                      pratsv = prat
                      pjrat = 1
                      prat = (b1 + c1 * cea%AcAt) / (1 + a1l * cea%AcAt)
                      if (cea%legacy_mode .and. cea%Debugf) &
                           write(IOOUT, '(" NEW PC = ", f10.2, 2x, "NEW ACAT = ", f9.6, 2x, "PJRAT =", &
                           & f10.7, " PRACAT =", f10.7)') pcpa, cea%AcAt, pjrat, pracat
                   end do
                end if

                cea%Hp = .true.
                cea%Sp = .false.
                niter = niter + 1
                cea%Isv = 0
                cea%ipt = 2
                ipp = 2
                call SETEN(cea)
                go to 250
             end if

950          if (cea%Trnspt) call TRANP(cea)
             if (cea%ipt == cea%Nfz) cea%Eql = seql

1000         itnum = 0
             if (cea%Nsub > i01) then
                isub = isub + 1
                if (isub <= cea%Nsub) go to 750
                isub = 1
                cea%Nsub = -cea%Nsub
                if (cea%Isup <= cea%Nsup) go to 750
                cea%Area = .false.
                go to 1150
             end if

1050         cea%Isup = cea%Isup + 1
             itnum = 0
             if (cea%Isup <= cea%Nsup) go to 750
             cea%Isup = isupsv
             cea%Area = .false.
             go to 1150

             ! TEST FOR OUTPUT -- SCHEDULES COMPLETE OR NPT=Ncol
1100         cea%Isv = cea%ipt
             if (cea%ipt <= nptth .or. mod(cea%ipt - nptth, Ncol - nptth) > 0) go to 1200

1150         if (.not. cea%Eql) then
                if (cea%Nfz <= 1) then
                   pfz%Cpr = cprf
                   pfz%Gammas = cprf / (cprf - 1 / pfz%Wm)
                end if
             end if

             call RKTOUT(cea, it)

             cea%Iplt = cea%ipt + nptth * ((cea%ipt - 1) / Ncol)
             if (.not. cea%Page1) then
                cea%Iplt = cea%Iplt - 2
                if (cea%Iopt /= 0) cea%Iplt = cea%Iplt - 1
                cea%Iplt = min(cea%Iplt, 500)
             else
                cea%Page1 = .false.
             end if
             iplte = max(iplte, cea%Iplt)
             dlnp = 1
             if (cea%Tt == 0) cea%Area = .false.
             if (cea%legacy_mode .and. .not. cea%Eql .and. cea%Tt == 0.) &
                  write(IOOUT, '(/" WARNING!!  CALCULATIONS WERE STOPPED BECAUSE NEXT ", &
                  & "POINT IS MORE", /, " THAN 50 K BELOW THE TEMPERATURE", &
                  & " RANGE OF A CONDENSED SPECIES (ROCKET)")')
             if (cea%Isv == 0) then
                ! PCP, SUBAR, AND SUPAR SCHEDULES COMPLETED
                if (cea%Nsub < 0) cea%Nsub = -cea%Nsub
                if (.not. cea%Froz .or. .not. cea%Eql) go to 1300
                ! SET UP FOR FROZEN.
                if (cea%Eql) cea%Iplt = iplt1
                cea%Eql = .false.
                cea%Page1 = .true.
                call SETEN(cea)
                cea%Tt = pfz%Ttt
                ipp = cea%Nfz
                if (cea%Nfz == cea%ipt) go to 1150
                cea%ipt = cea%Nfz
                cea%Enn = 1 / pfz%Wm
                if (cea%Nfz == 1) go to 450
                if (cea%Nsub > 0) then
                   cea%Nsub = -cea%Nsub
                   if (cea%legacy_mode) write(IOOUT, '(/" WARNING!!  FOR FROZEN PERFORMANCE, SUBSONIC AREA ", /, &
                        & " RATIOS WERE OMITTED SINCE nfz IS GREATER THAN 1 (ROCKET)")')
                end if

                if (pfz%App < pth%App) then
                   if (cea%legacy_mode) write(IOOUT, '(/" WARNING!!  FREEZING IS NOT ALLOWED AT A SUBSONIC ", &
                        & "PRESSURE RATIO FOR nfz GREATER"/" THAN 1. FROZEN ", &
                        & "PERFORMANCE CALCULATIONS WERE OMITTED (ROCKET)")')
                   go to 1300
                else
                   if (cea%Nfz < cea%Npp) go to 1200
                   go to 700
                end if
             else
                if (cea%legacy_mode .and. cea%Eql) write(IOOUT, '(////)')
             end if

             ! SET INDICES AND ESTIMATES FOR NEXT POINT.
1200         cea%ipt = cea%ipt + 1
             if (cea%Eql .or. (cea%Isv == -i12 .and. .not. seql)) then
                ! THE FOLLOWING STATEMENT WAS ADDED TO TAKE CARE OF A SITUATION
                ! WHERE EQLBRM WENT SINGULAR WHEN STARTING FROM ESTIMATES WHERE
                ! BOTH SOLID AND LIQUID WERE INCLUDED.  JULY 27, 1990.
                if (cea%Jliq /= 0 .and. cea%Isv > 0) cea%Isv = 0
                call SETEN(cea)
             end if

             p => cea%points(cea%iOF, cea%ipt)

             do
                ipp = ipp + 1
                if (cea%ipt > nptth) then
                   if (cea%Area) then
                      p%App = EXP(appl)
                   else
                      p%App = cea%Pcp(ipp - nptth)
                      if (cea%Fac) p%App = p%App * Pinf / p1%Ppp
                      if (.not. cea%Eql .and. p%App < pfz%App) then
                         if (cea%legacy_mode) &
                              write(IOOUT, '(/, " WARNING!! FOR FROZEN PERFORMANCE, POINTS WERE OMITTED", &
                              & " WHERE THE ASSIGNED", /, &
                              & " PRESSURE RATIOS WERE LESS THAN ", &
                              & "THE VALUE AT POINT nfz =", i3, " (ROCKET)")') cea%Nfz
                         cycle
                      end if
                   end if
                   cea%Pp = pinf / p%App
                   if (cea%Fac) then
                      if (cea%Area) then
                         if (isub <= cea%Nsub .and. isub > i01 .and. aratio >= cea%points(cea%iOF, 2)%AeAt) then
                            if (cea%legacy_mode) &
                                 write(IOOUT, '(/" WARNING!!  ASSIGNED subae/at =", f10.5, " IS NOT ", &
                                 & "PERMITTED TO BE GREATER"/" THAN ac/at =", f9.5, &
                                 & ".  POINT OMITTED (ROCKET)")') aratio, cea%points(cea%iOF, 2)%AeAt
                            cea%ipt = cea%ipt - 1
                            go to 1000
                         end if
                      else if (cea%ipt > nptth .and. cea%Pcp(ipp-3) < p1%Ppp / p2%Ppp) then
                         if (cea%legacy_mode) &
                              write(IOOUT, '(/" WARNING!!  ASSIGNED pip =", F10.5, &
                              & " IS NOT PERMITTED"/" TO BE LESS THAN  Pinj/Pc =", f9.5, &
                              & ". POINT OMITTED", " (ROCKET)")') cea%Pcp(ipp-3), p1%Ppp / p2%Ppp
                         cea%ipt = cea%ipt - 1
                         go to 650
                      end if
                   end if
                end if
                exit
             end do

             go to 250

1300         cea%ipt = 1

             ! CHECK FOR COMPLETED SCHEDULES -
             ! 1) CHAMBER PRESSURES(IP = NP)
             ! 2) CHAMBER TEMPERATURES(IT = NT)
             ! 3) O/F VALUES(IOF = NOF)
             if (ip == cea%Np .and. it == cea%Nt .and. iof == cea%Nof) exit iofLoop
             if (cea%legacy_mode) write(IOOUT, '(////)')
             call SETEN(cea)
             cea%Tt = p12%Ttt
          end do
       end do
    end do iofLoop

    cea%Iplt = max(cea%Iplt, iplte)

    return
  end subroutine ROCKET


  subroutine read_libraries(cea)
    use mod_types
    implicit none

    class(CEA_Core_Problem), intent(inout):: cea
    integer:: i, j
    real(8):: xi, xln

    if (cea%invalid_case) return

    cea%Nonly = cea%Nonly_in
    cea%Prod(:) = cea%Prod_in(:)

    if (cea%Ions) then
       if (cea%Elmt(cea%Nlm) /= 'E') then
          cea%Nlm = cea%Nlm + 1
          cea%Elmt(cea%Nlm) = 'E'

          cea%B0p(cea%Nlm, 1) = 0
          cea%B0p(cea%Nlm, 2) = 0
       end if
    else if (cea%Elmt(cea%Nlm) == 'E') then
       cea%Nlm = cea%Nlm - 1
    end if

    cea%Jray(1:cea%Nreac) = 0

    call SEARCH(cea)

    if (cea%Ngc == 0) then
       cea%invalid_case = .true.
       return
    end if

    if (cea%Trnspt) call READTR(cea)

    ! INITIAL ESTIMATES
    cea%Npr = 0
    cea%Gonly = .true.
    cea%Enn = 0.1d0
    cea%Ennl = -2.3025851
    cea%Sumn = cea%Enn
    xi = cea%Ng
    if (xi == 0.) xi = 1
    xi = cea%Enn/xi
    xln = log(xi)

    do concurrent (i = 1:cea%Nof)
       cea%points(i, 1)%En(cea%Ng+1:cea%Ng+cea%Nc) = 0
       cea%points(i, 1)%En(1:cea%Ng) = xi
    end do

    cea%Enln(cea%Ng+1:cea%Ng+cea%Nc) = 0
    cea%Enln(1:cea%Ng) = xln

    if (cea%Nc /= 0 .and. cea%Nsert /= 0) then
       innerLoop: do i = 1, cea%Nsert
          do j = cea%Ngc, cea%Ngp1, - 1
             if (cea%Prod(j) == cea%ensert(i)) then
                cea%Npr = cea%Npr + 1
                cea%Jcond(cea%Npr) = j
                if (.not. cea%Short) write(cea%io_log, '(1X, A16, "INSERTED")') cea%Prod(j)
                cycle innerLoop
             end if
          end do
          write(cea%io_log, '(/" WARNING!!!", A16, "NOT FOUND FOR INSERTION")') cea%ensert(i)
       end do innerLoop
    end if

    return
  end subroutine read_libraries


  subroutine SEARCH(cea)
    !***********************************************************************
    ! SEARCH THERMO.LIB FOR THERMO DATA FOR SPECIES TO BE CONSIDERED.
    !***********************************************************************
    use mod_types
    implicit none

    type(CEA_Core_Problem), intent(inout):: cea

    ! LOCAL VARIABLES
    character(6):: date(maxNgc)
    character(2):: el(5)
    character(15):: sub
    integer:: i, j, k, ii
    integer:: i5, ifaz, itot, nall, ne, nint, ntgas, ntot
    real(8):: b(5), t1, t2, thermo(9, 3)
    integer:: io_thermo
    logical:: is_opened

    cea%Nc = 0
    ne = 0
    do i = 1, cea%Nlm
       cea%Jx(i) = 0
    end do
    do j = 1, maxNgc
       cea%S(j) = 0
       cea%H0(j) = 0
       cea%Deln(j) = 0
       do i = 1, cea%Nlm
          cea%A(i, j) = 0
       end do
    end do
    ! READ TEMPERATURE RANGES FOR COEFFICIENTS OF GASEOUS SPECIES.
    ! SOME DEFINITIONS:
    !   NTGAS = NUMBER OF GASEOUS SPECIES IN THERMO.LIB.
    !   NTOT =  NTGAS PLUS NUMBER OF TEMPERATURE INTERVALS FOR CONDENSED.
    !   NALL =  NTOT PLUS THE NUMBER OF REACTANT SPECIES IN THERMO.LIB.
    !   NG =    NUMBER OF GASES WITH STORED COEFFICIENTS.
    !   NC =    NUMBER OF CONDENSED INTERVALS WITH STORED COEFFICIENTS.
    !   NGC =    NG + NC
    !   THDATE = DATE READ FROM THERMO.INP FILE
    open(newunit = io_thermo, file = cea%filename_thermo_lib, status = 'old', form = 'unformatted', action = 'read')
    read(io_thermo) cea%Tg, ntgas, ntot, nall, cea%Thdate
    cea%Ngc = 1
    cea%Nc = 1

    ! BEGIN LOOP FOR READING SPECIES DATA FROM THERMO.LIB.
    outerLoop: do itot = 1, ntot
       if (itot > ntgas) then
          read(io_thermo) sub, nint, date(cea%Ngc), (el(j), b(j), j = 1, 5), cea%Ifz(cea%Nc), &
               cea%Temp(1, cea%Nc), cea%Temp(2, cea%Nc), cea%Mw(cea%Ngc), (cea%Cft(k, cea%Nc), k = 1, 9)
       else
          read(io_thermo) sub, nint, date(cea%Ngc), (el(j), b(j), j = 1, 5), ifaz, T1, T2, cea%Mw(cea%Ngc), thermo
       end if
       if (cea%Nonly /= 0) then
          i = 1
20        if (cea%Prod(i) /= sub .and. '*' // cea%Prod(i) /= sub) then
             i = i + 1
             if (i <= cea%Nonly) go to 20
             cycle outerLoop
          else
             if (sub == cea%Prod(cea%Ngc-1)) then
                cea%Nonly = cea%Nonly + 1
                do k = cea%Nonly, i+1, -1
                   cea%Prod(k) = cea%Prod(k-1)
                end do
             else
                cea%Prod(i) = cea%Prod(cea%Ngc)
             end if
             cea%Prod(cea%Ngc) = sub
          end if
       else if (cea%Nomit /= 0) then
          do i = 1, cea%Nomit
             if (cea%Omit(i) == sub .or. '*' // cea%Omit(i) == sub) cycle outerLoop
          end do
       end if
       innerLoop: do k = 1, 5
          if (b(k) == 0) exit
          do i = 1, cea%Nlm
             if (cea%Elmt(i) == el(k)) then
                cea%A(i, cea%Ngc) = b(k)
                cycle innerLoop
             end if
          end do
          do j = 1, cea%Nlm
             cea%A(j, cea%Ngc) = 0
          end do
          cycle outerLoop
       end do innerLoop
       cea%Prod(cea%Ngc) = sub
       if (itot > ntgas) then
          cea%Nc = cea%Nc + 1
          if (cea%Nc > maxNc) go to 400
       else
          cea%Ng = cea%Ngc
          if (cea%Ng > maxNg) go to 400
          do i = 1, 3
             do j = 1, 9
                cea%Coef(j, cea%Ng, i) = thermo(j, i)
             end do
          end do
          ! IF SPECIES IS AN ATOMIC GAS, STORE INDEX IN JX
          if (b(2) == 0 .and. b(1) == 1) then
             do i = 1, cea%Nlm
                if (cea%Elmt(i) == el(1)) then
                   ne = ne + 1
                   cea%Jx(i) = cea%Ngc
                   cea%Jcm(i) = cea%Ngc
                   go to 150
                end if
             end do
          end if
       end if
150    cea%Ngc = cea%Ngc + 1
       if (cea%Ngc > maxNgc) go to 400
    end do outerLoop

    close(io_thermo)

    ! FINISHED READING THERMO DATA FROM I/O UNIT io_thermo.
    cea%Ifz(cea%Nc) = 0
    cea%Nc = cea%Nc - 1
    cea%Ngc = cea%Ngc - 1
    cea%Ngp1 = cea%Ng + 1
    if (cea%Ngc < cea%Nonly) then
       do k = cea%Ngc+1, cea%Nonly
          write(cea%io_log, '(/" WARNING!!  ", a15, " NOT A PRODUCT IN thermo.lib FILE (SEARCH)")') cea%Prod(k)
       end do
    end if
    ! FIND MISSING ELEMENTS (IF ANY) FOR COMPONENTS
    cea%Nspx = cea%Ngc
    if (ne < cea%Nlm) then
       do i = 1, cea%Nlm
          if (cea%Nspx > maxNgc) go to 400
          if (cea%Jx(i) == 0) then
             cea%Nspx = cea%Nspx + 1
             do k = 1, cea%Nlm
                cea%A(k, cea%Nspx) = 0
             end do
             cea%A(i, cea%Nspx) = 1
             cea%Prod(cea%Nspx) = cea%Elmt(i)
             do k = 1, 100
                if (cea%Elmt(i) == atomic_symbol(k)) then
                   cea%Mw(cea%Nspx) = atomic_mass(k)
                   cea%Atwt(i) = atomic_mass(k)
                   cea%Cp(cea%Nspx) = 2.5d0
                   go to 210
                end if
             end do
210          cea%Jx(i) = cea%Nspx
             cea%Jcm(i) = cea%Nspx
          end if
       end do
    end if
    ! ARE ALL ELEMENTS IN PRODUCT SPECIES?
    outerLoop2: do i = 1, cea%Nlm
       do j = 1, cea%Ngc
          if (cea%A(i, j) /= 0) cycle outerLoop2
          ii = i
       end do
       write(cea%io_log, '(/" PRODUCT SPECIES CONTAINING THE ELEMENT", a3, " MISSING", &
            & //, 13x, "FATAL ERROR (SEARCH)")') cea%Elmt(ii)
       cea%Ngc = 0
       return
    end do outerLoop2
    ! WRITE POSSIBLE PRODUCT LIST
    if (.not. cea%Short) then
       write(cea%io_log, '(/2x, "SPECIES BEING CONSIDERED IN THIS SYSTEM", &
            & /" (CONDENSED PHASE MAY HAVE NAME LISTED SEVERAL TIMES)", &
            & /"  LAST thermo.inp UPDATE: ", a10, /)') cea%Thdate
       do i = 1, cea%Ngc, 3
          i5 = i + 2
          if (cea%Ngc < i5) i5 = cea%Ngc
          write(cea%io_log, '(3(2X, A6, 2X, A15))') (date(j), cea%Prod(j), j = i, i5)
       end do
    end if
    return

400 write(cea%io_log, '(/" INSUFFICIENT STORAGE FOR PRODUCTS-SEE RP-1311,", /"   PART 2, PAGE 39. (SEARCH)")')
    cea%Ngc = 0

    inquire(io_thermo, opened = is_opened)
    if (is_opened) close(io_thermo)

    return
  end subroutine SEARCH


  subroutine READTR(cea)
    ! SEARCH FOR TRANSPORT PROPERTIES FOR THIS CHEMICAL SYSTEM
    use mod_types
    implicit none

    type(CEA_Core_Problem), intent(inout):: cea

    character(16):: bin(2, 40), pure(6), spece(2)
    integer:: i, j, k, jj(2), jk, ir, lineb, npure, nrec
    real(8):: trdata(36)
    integer:: io_transport

    if (allocated(cea%transport_properties)) return ! Already done.

    open(newunit = io_transport, file = trim(cea%filename_trans_lib), status = 'old', form = 'unformatted', action = 'read')

    cea%Ntape = 0
    npure = 0
    lineb = 1

    if (.not. cea%Short) write(cea%io_log, '(/" SPECIES WITH TRANSPORT PROPERTIES"//8X, "PURE SPECIES"/)')

    read(io_transport) nrec

    allocate(cea%transport_properties(nrec))

    do ir = 1, nrec
       read(io_transport) spece, trdata
       k = 1
       do j = 1, cea%Ng
          if (spece(k) == cea%Prod(j) .or. '*' // spece(k) == cea%Prod(j)) then
             jj(k) = j
             if (k == 2) then
                ! STORE NAMES FOR BINARIES IN BIN ARRAY.
                do k = 1, 2
                   bin(k, lineb) = spece(k)
                end do
                lineb = lineb + 1
                exit
             else
                jj(2) = j
                if (spece(2) == ' ') then
                   ! WRITE NAMES FOR PURE SPECIES.
                   npure = npure + 1
                   pure(npure) = spece(1)
                   exit
                else
                   k = 2
                   cycle
                end if
             end if
          end if
       end do

       if (j <= cea%Ng) then
          cea%Ntape = cea%Ntape + 1
          cea%transport_properties(cea%Ntape)%j = jj
          cea%transport_properties(cea%Ntape)%data = reshape(trdata, [6, 3, 2])
       end if

       if (npure /= 0 .and. (npure >= 6 .or. ir >= nrec)) then
          if (.not. cea%Short) write(cea%io_log, '(4(2x, A16))') (pure(jk), jk = 1, npure)
          npure = 0
       end if
    end do

    lineb = lineb - 1

    if (.not. cea%Short) then
       write(cea%io_log, '(/"     BINARY INTERACTIONS"/)')
       do j = 1, lineb
          write(cea%io_log, '(5X, 2A16)') (bin(i, j), i = 1, 2)
       end do
    end if

    write(cea%io_log, *)

    close(io_transport)

    return
  end subroutine READTR



  subroutine SETEN(cea)
    !***********************************************************************
    ! USE COMPOSITIONS FROM PREVIOUS POINT AS INITIAL ESTIMATES FOR
    ! CURRENT POINT NPT.  IF -
    !  ISV>0  USE COMPOSITIONS FROM POINT ISV.
    !  ISV<0  SAVE COMPOSITIONS FROM POINT -ISV FOR POSSIBLE LATER USE.
    !         ALSO USE COMPOSITIONS FROM POINT -ISV FOR NPT.
    !  ISV=0  USE COMPOSITIONS SAVED WHEN ISV<0.
    !***********************************************************************
    use mod_types
    implicit none

    type(CEA_Core_Problem), intent(inout):: cea

    ! LOCAL VARIABLES
    integer:: j
    type(CEA_Point), pointer:: p
    type(CEA_Point), pointer:: psv

    p => cea%points(cea%iOF, cea%ipt)

    if (cea%Isv > 0) then
       ! USE COMPOSITIONS FROM PREVIOUS POINT
       psv => cea%points(cea%iOF, cea%Isv)
       do j = 1, cea%Ngc
          p%En(j) = psv%En(j)
       end do
    else if (cea%Isv < 0) then
       ! FIRST T--SAVE COMPOSITIONS FOR FUTURE POINTS WITH THIS T
       cea%Isv = -cea%Isv
       psv => cea%points(cea%iOF, cea%Isv)
       cea%Tsave = psv%Ttt
       cea%Ensave = cea%Enn
       cea%Enlsav = cea%Ennl
       cea%lsav = cea%lsave
       do j = 1, cea%Ng
          cea%Sln(j) = cea%Enln(j)
       end do
       do j = 1, cea%Ng
          p%En(j) = psv%En(j)
       end do
       cea%Npr = 0
       do j = cea%Ngp1, cea%Ngc
          cea%Sln(j) = psv%En(j)
          p%En(j) = cea%Sln(j)
          if (cea%Jliq == j) then
             p%En(cea%Jsol) = psv%En(cea%Jsol) + psv%En(cea%Jliq)
             p%En(cea%Jliq) = 0
             cea%Jsol = 0
             cea%Jliq = 0
             cea%Tsave = cea%Tsave - 5
             cea%Tt = cea%Tsave
             cea%Sln(j) = 0
          else if (p%En(j) > 0) then
             cea%Npr = cea%Npr + 1
             cea%Jcond(cea%Npr) = j
          end if
       end do
    else ! if (cea%Isv == 0) then
       ! NEXT POINT FIRST T IN SCHEDULE, USE PREVIOUS COMPOSITIONS FOR THIS T
       cea%Jsol = 0
       cea%Jliq = 0
       cea%Enn = cea%Ensave
       cea%Ennl = cea%Enlsav
       cea%lsave = cea%lsav
       cea%Npr = 0
       do j = cea%Ngp1, cea%Ngc
          p%En(j) = cea%Sln(j)
          if (p%En(j) > 0.d0) then
             cea%Npr = cea%Npr + 1
             cea%Jcond(cea%Npr) = j
          end if
       end do
       do j = 1, cea%Ng
          p%En(j) = 0
          cea%Enln(j) = cea%Sln(j)
          if (cea%Sln(j) /= 0) then
             if ((cea%Enln(j) - cea%Ennl + 18.5) > 0) p%En(j) = exp(cea%Enln(j))
          end if
       end do
       if (.not. cea%Tp) cea%Tt = cea%Tsave
       cea%Sumn = cea%Enn
    end if
  end subroutine SETEN



  subroutine SHCK(cea)
    !***********************************************************************
    ! PRIMARY ROUTINE FOR SHOCK PROBLEMS.
    !***********************************************************************
    use mod_types
    use mod_legacy_io
    implicit none

    type(CEA_Core_Problem), intent(inout):: cea

    ! LOCAL VARIABLES
    character(1):: cr12, cr52
    integer:: i, iof, it1, it2, itr, j, n
    logical:: refl, seql, srefl
    real(8):: ax, axx, b2, cormax, gg, hs, m2m1(Ncol), mis(13), mu12rt, p1, p21, &
         p21l, p2p1(Ncol), pmn, rho12, rho52, rrho(Ncol), sg(78), T1, T21, &
         t21l, t2t1(Ncol), ttmax, u1u2(Ncol), uis(13), utwo(Ncol), uu, wmx, ww

    type(CEA_Point), pointer:: p, p_prev

    if (cea%Trace == 0) then
       if (cea%legacy_mode) then
          cea%Trace = 5e-9
       else
          cea%Trace = 5d-9
       end if
    end if

    cea%Tp = .true.
    cea%Cpmix = 0
    srefl = .false.
    if (cea%legacy_mode .and. .not. cea%Short) then
       write(IOOUT, '(/"   *** INPUT FOR SHOCK PROBLEMS ***")')
       write(IOOUT, '(/" INCDEQ =", L2, "   REFLEQ =", L2, "   INCDFZ =", L2, "    REFLFZ =", L2)') &
            cea%Incdeq, cea%Refleq, cea%Incdfz, cea%Reflfz
    end if
    if (cea%Refleq .or. cea%Reflfz) srefl = .true.
    seql = cea%Incdeq
    if (cea%T(1) == 0) cea%T(1) = cea%Rtemp(1)
    do i = 1, cea%Nsk
       uis(i) = cea%points(1, i)%U1
       mis(i) = cea%points(1, i)%Mach1
       if (cea%points(1, i)%Mach1 == 0 .and. cea%points(1, i)%U1 == 0) exit
    end do
    if (cea%Nsk > Ncol) then
       write(IOOUT, '(/" WARNING!!  ONLY ", I2, " u1 OR mach1 VALUES ALLOWED (SHCK)")') Ncol
       cea%Nsk = Ncol
    end if
    if (cea%legacy_mode .and. .not. cea%Short) then
       write(IOOUT, '(/1p, " U1 =   ", 5E13.6, /(8X, 5E13.6))') (cea%points(1, i)%U1, i = 1, cea%Nsk)
       write(IOOUT, '(/1p, " MACH1 =", 5E13.6, /(8X, 5E13.6))') (cea%points(1, i)%Mach1, i = 1, cea%Nsk)
    end if

    do iof = 1, cea%Nof
       cea%Oxfl = cea%Oxf(iof)
       call NEWOF(cea)
       cea%Incdeq = seql

       do
          refl = .false.
          it2 = 2
          it1 = 1
          cea%Pp = cea%P(1)
          cea%Tt = cea%T(1)

          if (.not. cea%Incdeq) then
             ! FROZEN
             do i = 1, cea%Nsk
                p => cea%points(iof, i)
                p%Dlvtp = 1
                p%Dlvpt = -1
             end do
          end if

          do i = 1, cea%Nsk
             p => cea%points(iof, i)

             p%Ppp = cea%P(i)
             p%Ttt = cea%T(i)

             if (i > 1) then
                p_prev => cea%points(iof, i-1)
                if (p%Ppp == 0) p%Ppp = p_prev%Ppp
                if (p%Ttt == 0) p%Ttt = p_prev%Ttt
                p%Ssum = p_prev%Ssum
                p%Hsum = p_prev%Hsum
             end if

             if (i == 1 .or. p%Ttt /= cea%Tt .or. p%Ppp /= cea%Pp) then
                cea%Pp = p%Ppp
                cea%Tt = p%Ttt
                if (cea%Tt >= cea%Tg(1) * 0.8d0) then
                   cea%ipt = i
                   call HCALC(cea)
                   p%Hsum = cea%Hsub0
                else
                   if (cea%legacy_mode) &
                        write(IOOUT, '(/" TEMPERATURE=", E12.4, " IS OUT OF EXTENDED RANGE ", "FOR POINT", I5, " (SHCK)")') &
                        cea%Tt, i
                   return
                end if
             end if

             if (cea%Cpmix /= 0) cea%Gamma1 = cea%Cpmix / (cea%Cpmix - 1/cea%Wmix)
             cea%A1 = sqrt(R0 * cea%Gamma1 * cea%Tt / cea%Wmix)
             if (p%U1 == 0) p%U1 = cea%A1 * p%Mach1
             if (p%Mach1 == 0) p%Mach1 = p%U1 / cea%A1
             p%Wm = cea%Wmix
             p%Cpr = cea%Cpmix
             p%Gammas = cea%Gamma1
             p%Vlm = R0 * cea%Tt / (cea%Wmix * cea%Pp)
          end do

          if (cea%legacy_mode) then
             ! OUTPUT--1ST CONDITION
             write(IOOUT, '(////25X, "SHOCK WAVE PARAMETERS ASSUMING")')
             if (.not. cea%Incdeq) then
                write(IOOUT, '(/, 17X, " FROZEN COMPOSITION FOR INCIDENT SHOCKED CONDITI1ONS"//)')
             else
                write(IOOUT, '(/, 16X, " EQUILIBRIUM COMPOSITION FOR INCIDENT SHOCKED CONDITIONS"//)')
             end if
             cea%Eql = .false.
             call OUT1(cea, IOOUT)
             write(IOOUT, '(/" INITIAL GAS (1)")')
             cea%fmt(4) = '13'
             cea%fmt(5) = ' '
             cea%fmt(7) = '4,'
             write(IOOUT, cea%fmt) 'MACH NUMBER1   ', (cea%points(iof, j)%Mach1, j = 1, cea%Nsk)
             cea%fmt(7) = '2,'
             write(IOOUT, cea%fmt) 'U1, M/SEC      ', (cea%points(iof, j)%U1, j = 1, cea%Nsk)
             call OUT2(cea, cea%Nsk, IOOUT)
          end if

          ! BEGIN CALCULATIONS FOR 2ND CONDITION
          if (cea%Incdeq) cea%Eql = .true.

          cea%ipt = 1
          do
             p => cea%points(iof, cea%ipt)

             cea%Gamma1 = p%Gammas
             uu = p%U1
             wmx = p%Wm
             p1 = p%Ppp
             T1 = p%Ttt
             hs = p%Hsum
             if (refl) uu = u1u2(cea%ipt)
             mu12rt = wmx * uu**2 / (R0 * T1)
             if (refl) then
                ! REFLECTED--SUBSCRIPTS 2=1, 5=2, P52=P21
                T21 = 2
                b2 = (-1 - mu12rt - T21) / 2
                p21 = -b2 + sqrt(b2**2 - T21)
             else
                p21 = (2 * cea%Gamma1 * p%Mach1**2 - cea%Gamma1 + 1) / (cea%Gamma1 + 1)
                ! THE FOLLOWING IMPROVED FORMULATION FOR THE INITIAL ESTIMATE FOR THE
                ! 2ND CONDITION WAS MADE AND TESTED BY S. GORDON 7/10/89.
                if (.not. cea%Eql) then
                   T21 = p21 * (2 / p%Mach1**2 + cea%Gamma1 - 1) / (cea%Gamma1 + 1)
                else
                   cea%Pp = p21 * p1
                   cea%Tp = .false.
                   cea%Hp = .true.
                   cea%Hsub0 = hs + uu**2 / (2 * R0)

                   call EQLBRM(cea)

                   T21 = p%Ttt / T1
                   cea%Hp = .false.
                   cea%Tp = .true.
                end if
             end if
             p21l = log(p21)
             ttmax = 1.05 * cea%Tg(4) / T1
             T21 = min(T21, ttmax)
             t21l = log(T21)

             ax = 1
             do itr = 1, 60
                if (cea%legacy_mode .and. cea%Shkdbg) &
                     write(IOOUT, '(/" ITR NO.=", I3, 3X, "P", I1, "/P", I1, " =", F9.4, 3X, "T", I1, &
                     & "/T", I1, " =", F9.4, "   RHO2/RHO1 =", F9.6)') itr, it2, it1, p21, it2, it1, T21, rho52
                cea%Tt = T21 * T1
                cea%Pp = p21 * p1
                if (.not. cea%Eql) then
                   ! FROZEN
                   cea%Tln = log(cea%Tt)
                   if (.not. cea%Incdeq) then
                      call HCALC(cea)
                      if (cea%Tt == 0) exit
                      p%Hsum = cea%Hsub0
                      p%Cpr = cea%Cpmix
                   else
                      call CPHS(cea)
                      p%Cpr = cea%Cpsum
                      p%Hsum = 0
                      do j = 1, cea%Ng
                         p%Hsum = p%Hsum + cea%H0(j) * p%En(j)
                      end do
                      p%Hsum = p%Hsum * cea%Tt
                   end if
                else
                   call EQLBRM(cea)

                   if (cea%Tt == 0) exit
                end if
                rho12 = wmx * T21 / (p%Wm * p21)
                gg = rho12 * mu12rt
                rho52 = 1 / rho12
                if (refl) gg = -mu12rt * rho52 / (rho52 - 1)**2
                cea%G(1, 1) = -gg * p%Dlvpt - p21
                cea%G(1, 2) = -gg * p%Dlvtp
                cea%G(1, 3) = p21 - 1 + gg - mu12rt
                if (refl) cea%G(1, 3) = p21 - 1 + gg * (rho52 - 1)
                gg = gg * T1 / wmx
                if (.not. refl) gg = gg * rho12
                cea%G(2, 1) = -gg * p%Dlvpt + cea%Tt * (p%Dlvtp - 1) / p%Wm
                cea%G(2, 2) = -gg * p%Dlvtp - cea%Tt * p%Cpr
                gg = 1 - rho12**2
                if (refl) gg = (rho52 + 1) / (rho52 - 1)
                cea%G(2, 3) = p%Hsum - hs - uu**2 * gg / (2 * R0)
                cea%X(3) = cea%G(1, 1) * cea%G(2, 2) - cea%G(1, 2) * cea%G(2, 1)
                cea%X(1) = (cea%G(1, 3) * cea%G(2, 2) - cea%G(2, 3) * cea%G(1, 2)) / cea%X(3)
                cea%X(2) = (cea%G(1, 1) * cea%G(2, 3) - cea%G(2, 1) * cea%G(1, 3)) / cea%X(3)
                if (cea%legacy_mode .and. cea%Shkdbg) then
                   write(IOOUT, '(/" G(I,J)  ", 3E15.8)') cea%G(1, 1), cea%G(1, 2), cea%G(1, 3)
                   write(IOOUT, '(/" G(I,J)  ", 3E15.8)') cea%G(2, 1), cea%G(2, 2), cea%G(2, 3)
                   write(IOOUT, '(/" X       ", 2E15.8)') cea%X(1), cea%X(2)
                   write(IOOUT, '(/" HSUM HS UU U2 ", 4E15.8)') p%Hsum, hs, uu, uu * rho12
                end if
                ax = abs(cea%X(1))
                axx = abs(cea%X(2))
                if (axx > ax) ax = axx

                if (ax < 0.00005) exit

                cormax = 0.40546511
                if (itr > 4) cormax = 0.22314355
                if (itr > 12) cormax = 0.09531018
                if (itr > 20) cormax = 0.04879016
                ax = ax / cormax
                if (ax > 1) then
                   cea%X(1) = cea%X(1) / ax
                   cea%X(2) = cea%X(2) / ax
                end if
                p21l = p21l + cea%X(1)
                t21l = t21l + cea%X(2)
                p21 = exp(p21l)
                T21 = exp(t21l)

                if (cea%legacy_mode .and. cea%Shkdbg) &
                     write(IOOUT, '(/" MAX.COR.=", e13.6, " X(1)=", e13.6, " X(2)=", e13.6)') cormax, cea%X(1), cea%X(2)

                if (itr == 1 .and. T21 >= ttmax) then
                   cea%Tt = 0
                   cea%ipt = cea%ipt - 1
                end if
             end do

             p => cea%points(iof, cea%ipt)

             if (.not. cea%Eql .or. cea%Tt /= 0) then
                if (cea%legacy_mode .and. itr > 60) then
                   write(IOOUT, '(/6x, " WARNING!!  NO CONVERGENCE FOR u1=", F8.1, &
                        & /"  ANSWERS NOT RELIABLE, SOLUTION MAY NOT EXIST (SHCK)")') p%U1
                end if

                if ((itr > 60) .or. (ax < 0.00005) .or. ((.not. cea%Eql) .and. (.not. cea%Incdeq) .and. (cea%Tt == 0))) then
                   ! CONVERGED OR TOOK 60 ITERATIONS WITHOUT CONVERGING.
                   ! STORE RESULTS.
                   rrho(cea%ipt) = rho52
                   m2m1(cea%ipt) = p%Wm / wmx
                   p2p1(cea%ipt) = p21
                   t2t1(cea%ipt) = T21
                   utwo(cea%ipt) = uu * rho12
                   u1u2(cea%ipt) = uu - utwo(cea%ipt)
                   if (cea%Tt >= cea%Tg(1) * 0.8d0 .and. cea%Tt <= cea%Tg(4) * 1.1d0) then
                      if (.not. cea%Eql) then
                         ! FROZEN
                         p%Ppp = cea%Pp
                         p%Ttt = cea%Tt
                         p%Gammas = p%Cpr / (p%Cpr - 1 / wmx)
                         p%Vlm = R0 * cea%Tt / (wmx * cea%Pp)
                         if (cea%Incdeq) then
                            p%Ssum = 0
                            do j = 1, cea%Ngc
                               pmn = cea%Pp * wmx * p%En(j)
                               if (p%En(j) > 0) p%Ssum = p%Ssum + p%En(j) * (cea%S(j) - log(pmn))
                            end do
                         end if
                      end if
                      go to 900
                   end if
                end if

                if (cea%legacy_mode) &
                     write(IOOUT, '(/" TEMPERATURE=", E12.4, " IS OUT OF EXTENDED RANGE ", &
                     & "FOR POINT", I5, " (SHCK)")') cea%Tt, cea%ipt
                cea%Tt = 0
             end if

             if (cea%ipt < 1) return

             cea%Nsk = cea%ipt

900          if (cea%Trnspt) then
                call TRANP(cea)
             end if

             cea%Isv = 0
             if (cea%ipt < cea%Nsk) cea%Isv = cea%ipt
             if (cea%ipt == 1) cea%Isv = -1
             cea%ipt = cea%ipt + 1

             if (cea%Eql .and. cea%ipt <= cea%Nsk) then
                call SETEN(cea)
             end if

             if (cea%ipt <= cea%Nsk) cycle
             cea%ipt = cea%Nsk

             if (cea%legacy_mode) then
                if (refl) then
                   if (cea%Eql) then
                      write(IOOUT, '(/" SHOCKED GAS (5)--REFLECTED--EQUILIBRIUM")')
                   else
                      write(IOOUT, '(/" SHOCKED GAS (5)--REFLECTED--FROZEN")')
                   end if
                   cr12 = '2'
                   cr52 = '5'
                else
                   if (cea%Eql) then
                      write(IOOUT, '(/" SHOCKED GAS (2)--INCIDENT--EQUILIBRIUM")')
                   else
                      write(IOOUT, '(/" SHOCKED GAS (2)--INCIDENT--FROZEN")')
                   end if
                   cr12 = '1'
                   cr52 = '2'
                end if
                cea%fmt(7) = '2,'
                write(IOOUT, cea%fmt) 'U' // cr52 // ', M/SEC      ', (utwo(j), j = 1, cea%ipt)

                call OUT2(cea, cea%ipt, IOOUT)

                if (cea%Trnspt) call OUT4(cea, cea%ipt, IOOUT)
                write(IOOUT, *)
                cea%fmt(7) = '3,'
                write(IOOUT, cea%fmt) 'P' // cr52 // '/P' // cr12 // '           ', (p2p1(j), j = 1, cea%ipt)
                write(IOOUT, cea%fmt) 'T' // cr52 // '/T' // cr12 // '           ', (t2t1(j), j = 1, cea%ipt)
                cea%fmt(7) = '4,'
                write(IOOUT, cea%fmt) 'M' // cr52 // '/M' // cr12 // '           ', (m2m1(j), j = 1, cea%ipt)
                write(IOOUT, cea%fmt) 'RHO' // cr52 // '/RHO' // cr12 // '       ', (rrho(j), j = 1, cea%ipt)
                cea%fmt(7) = '2,'
                if (.not. refl) write(IOOUT, cea%fmt) 'V2, M/SEC      ', (u1u2(j), j = 1, cea%ipt)
                if (refl) write(IOOUT, cea%fmt) 'U5+V2,M/SEC    ', (u1u2(j), j = 1, cea%ipt)

                if (.not. cea%Eql) then
                   ! WRITE FROZEN MOLE (OR MASS) FRACTIONS
                   cea%fmt(7) = '5,'
                   if (.not. cea%Incdeq) then
                      if (cea%Massf) then
                         write(IOOUT, '(/1x, A4, " FRACTIONS"/)') 'MASS'
                      else
                         write(IOOUT, '(/1x, A4, " FRACTIONS"/)') 'MOLE'
                         ww = wmx
                      end if
                      do n = 1, cea%Nreac
                         j = cea%Jray(n)
                         if (cea%Massf) ww = cea%Mw(j)
                         write(IOOUT, '(" ", A16, F8.5, 12F9.5)') cea%Prod(j), (cea%points(cea%iOF, i)%En(j) * ww, i = 1, cea%ipt)
                      end do
                   else
                      cea%Eql = .true.
                      call OUT3(cea, cea%ipt, IOOUT)
                      cea%Eql = .false.
                   end if
                else
                   call OUT3(cea, cea%ipt, IOOUT)
                end if
             end if

             cea%Iplt = min(cea%Iplt + cea%ipt, 500)

             if (srefl) then
                if (.not. refl) then
                   refl = .true.
                   it2 = 5
                   it1 = 2
                   cea%Eql = .true.
                   if (cea%Reflfz) then
                      cea%Eql = .false.
                      if (cea%Refleq) then
                         j = 1
                         do i = 1, cea%ipt
                            p => cea%points(iof, i)
                            sg(j)   = u1u2(i)
                            sg(j+1) = p%Wm
                            sg(j+2) = p%Ppp
                            sg(j+3) = p%Ttt
                            sg(j+4) = p%Hsum
                            sg(j+5) = p%Gammas
                            j = j + 6
                         end do
                      end if
                   end if
                   cea%ipt = 1
                   cycle
                else if (.not. cea%Eql .and. cea%Refleq) then
                   j = 1
                   do i = 1, cea%ipt
                      p => cea%points(iof, i)
                      u1u2(i) = sg(j)
                      p%Wm = sg(j+1)
                      p%Ppp = sg(j+2)
                      p%Ttt = sg(j+3)
                      p%Hsum = sg(j+4)
                      p%Gammas = sg(j+5)
                      j = j + 6
                   end do
                   cea%Eql = .true.
                   cea%ipt = 1
                   cycle
                end if
             end if

             exit
          end do

          if (.not. cea%Incdeq .or. .not. cea%Incdfz) exit

          cea%Incdeq = .false.
          cea%Eql = .false.
       end do

       do concurrent (i = 1:cea%Nsk)
          cea%points(iof, i)%U1 = uis(i)
          cea%points(iof, i)%Mach1 = mis(i)
       end do
    end do

    cea%Tp = .false.
    do n = 1, cea%Nreac
       cea%Rtemp(n) = cea%T(1)
    end do

    return
  end subroutine SHCK



  subroutine THERMP(cea)
    !***********************************************************************
    ! ASSIGNED THERMODYNAMIC STATES.  HP, SP, TP, UV, SV, AND TV PROBLEMS.
    !***********************************************************************
    use mod_types
    use mod_legacy_io
    implicit none

    type(CEA_Core_Problem), intent(inout):: cea

    ! LOCAL VARIABLES
    integer:: iof
    integer:: ip, it

    cea%Eql = .true.
    iofLoop: do iof = 1, cea%Nof
       cea%Oxfl = cea%Oxf(iof)
       call NEWOF(cea)
       ! SET ASSIGNED P OR VOLUME
       do ip = 1, cea%Np
          cea%Pp = cea%P(ip)
          ! SET ASSIGNED T
          do it = 1, cea%Nt
             cea%ipt = (ip - 1) * cea%Nt + it

             if (iof > 1 .or. cea%ipt > 1) then
                if (.not. cea%Tp .and. cea%Tt /= 0) cea%T(1) = cea%Tt

                if (cea%ipt == 2) then
                   cea%Isv = -1
                else if (cea%Nt /= 1 .and. (it == 1 .or. cea%Tt == 0.)) then
                   cea%Isv = 0
                else
                   cea%Isv = cea%ipt - 1
                end if
                call SETEN(cea)
             end if

             cea%Vv = cea%V(ip)
             cea%Tt = cea%T(it)
             call EQLBRM(cea)

             if (cea%ipt == 0) return

             if (cea%Trnspt .and. cea%Tt /= 0) call TRANP(cea)

             if (ip /= cea%Np .or. it /= cea%Nt .and. cea%Tt /= 0) then
                if (mod(cea%ipt, Ncol) /= 0) cycle
             end if

             if (cea%legacy_mode) then
                if (.not. cea%Hp) write(IOOUT, '(////15X, "THERMODYNAMIC EQUILIBRIUM PROPERTIES AT ASSIGNED")')
                if (cea%Hp) write(IOOUT, '(////9X, "THERMODYNAMIC EQUILIBRIUM COMBUSTION PROPERTIES AT ASSIGNED")')
                if (.not. cea%Vol) then
                   if (cea%Hp) write(IOOUT, '(/34X, " PRESSURES"/)')
                   if (cea%Tp) write(IOOUT, '(/27X, "TEMPERATURE AND PRESSURE"/)')
                   if (cea%Sp) write(IOOUT, '(/29X, "ENTROPY AND PRESSURE"/)')
                else
                   if (cea%Hp) write(IOOUT, '(/36X, " VOLUME"/)')
                   if (cea%Tp) write(IOOUT, '(/28X, "TEMPERATURE AND VOLUME"/)')
                   if (cea%Sp) write(IOOUT, '(/30X, "ENTROPY AND VOLUME"/)')
                end if
                call OUT1(cea, IOOUT)
                write(IOOUT, '(/" THERMODYNAMIC PROPERTIES"/)')
                call OUT2(cea, cea%ipt, IOOUT)
                if (cea%Trnspt) call OUT4(cea, cea%ipt, IOOUT)
                call OUT3(cea, cea%ipt, IOOUT)
             end if

             cea%Iplt = min((iof - 1) * cea%Nt * cea%Np + (ip - 1) * cea%Nt + it, 500)

             if ((ip == cea%Np .and. it == cea%Nt .or. cea%Tt == 0) .and. iof == cea%Nof) return

             if (cea%legacy_mode) write(IOOUT, '(////)')
          end do
       end do
    end do iofLoop

    return
  end subroutine THERMP



  subroutine TRANIN(cea)
    !***********************************************************************
    ! BRINGS IN AND SORTS OUT INPUT FOR TRANSPORT CALCULATIONS
    !***********************************************************************
    use mod_types
    implicit none

    type(CEA_Core_Problem), intent(inout):: cea

    ! LOCAL VARIABLES
    integer:: i, ii, inds(maxTr), ir, j, jtape(2), k, k1, k2, kt, kvc, l, loop, m, nms
    logical:: change, elc1, elc2, ion1, ion2, setx
    real(8):: coeff, debye, ekt, enel, enmin, ionic, lambda, omega, prop, qc, ratio, &
         stcf(maxTr, maxTr), stcoef(maxTr), te, testen, testot, total, &
         trc(6, 3, 2), wmols(maxTr), wmred, xsel, xss(maxTr)

    type(CEA_Point), pointer:: p !< current point
    type(CEA_Point), pointer:: p1

    p1 => cea%points(cea%iOF, 1)

    xss(:) = 0
    wmols(:) = 0
    inds(:) = 0
    nms = 0
    setx = .false.

    if (.not. cea%Eql) then
       if (.not. cea%Shock) then
          if (.not. setx) then
             setx = .true.
             cea%Nm = nms
             do i = 1, cea%Nm
                cea%Xs(i) = xss(i)
                cea%Wmol(i) = wmols(i)
                cea%Ind(i) = inds(i)
             end do
          end if
          go to 300
       else if (.not. cea%Incdeq) then
          if (cea%ipt <= 1) then
             cea%Nm = cea%Nreac
             do i = 1, cea%Nm
                j = cea%Jray(i)
                cea%Ind(i) = j
                cea%Wmol(i) = cea%Mw(j)
                cea%Xs(i) = p1%En(j) * p1%Wm
             end do
          end if
          go to 300
       end if
    end if

    ! PICK OUT IMPORTANT SPECIES
    p => cea%points(cea%iOF, cea%ipt)
    cea%Nm = 0
    total = 0
    enmin = 1.0d-11 / p%Wm
    testot = 0.999999999d0 / p%Wm
    do i = 1, cea%Lsave
       j = cea%Jcm(i)
       if (p%En(j) <= 0 .and. j <= cea%Ngc) then
          if ((cea%Enln(j) - cea%Ennl + 25.328436d0) > 0) p%En(j) = exp(cea%Enln(j))
       end if
       cea%Nm = cea%Nm + 1
       cea%Ind(cea%Nm) = j
       total = total + p%En(j)
       if (cea%Mw(j) < 1) enel = p%En(j)
       p%En(j) = -p%En(j)
    end do
    testen = 1 / (cea%Ng * p%Wm)

    if (total <= testot) then
       outerLoop1: do loop = 1, cea%Ng
          testen = testen / 10

          do j = 1, cea%Ng
             if (p%En(j) >= testen) then
                if (cea%Nm >= maxTr) then
                   if (cea%legacy_mode) &
                        write(IOOUT, '(/" WARNING!!  MAXIMUM ALLOWED NO. OF SPECIES", I3, " WAS USED IN ", &
                        & /" TRANSPORT PROPERTY CALCULATIONS FOR POINT", I3, "(TRANIN))")') cea%Nm, cea%ipt
                   exit outerLoop1
                else
                   total = total + p%En(j)
                   cea%Nm = cea%Nm + 1
                   cea%Ind(cea%Nm) = j
                   p%En(j) = -p%En(j)
                end if
             end if
          end do

          if (testen <= enmin) exit
       end do outerLoop1
    end if

    ! CALCULATE MOLE FRACTIONS FROM THE EN(J, NPT)
    do j = 1, cea%Ng
       p%En(j) = abs(p%En(j))
    end do
    do i = 1, cea%Nm
       j = cea%Ind(i)
       cea%Wmol(i) = cea%Mw(j)
       cea%Xs(i) = p%En(j) / total
    end do
    if (cea%ipt == cea%Nfz) then
       nms = cea%Nm
       do i = 1, cea%Nm
          xss(i) = cea%Xs(i)
          wmols(i) = cea%Wmol(i)
          inds(i) = cea%Ind(i)
       end do
       setx = .false.
    end if
    ! REWRITE REACTIONS TO ELIMINATE TRACE ELEMENTS
    cea%Nr = cea%Nm - cea%Lsave
    if (cea%Nr /= 0) then
       do k = 1, maxTr
          do m = 1, maxTr
             cea%Stc(k, m) = 0
          end do
       end do
       k = 1
       do i = cea%Lsave + 1, cea%Nm
          cea%Stc(k, i) = -1
          j = cea%Ind(i)
          do m = 1, cea%Lsave
             cea%Stc(k, m) = cea%A(m, j)
          end do
          k = k + 1
       end do
       do i = 1, cea%Nm
          if (cea%Xs(i) < 1.0d-10) then
             m = 1
             change = .false.

             do j = 1, cea%Nr
                coeff = cea%Stc(j, i)
                if (abs(coeff) > 1.0d-05) then
                   if (.not. change) then
                      change = .true.
                      do k = 1, cea%Nm
                         stcoef(k) = cea%Stc(j, k) / coeff
                      end do
                      cycle
                   else
                      do k = 1, cea%Nm
                         cea%Stc(j, k) = cea%Stc(j, k) / coeff - stcoef(k)
                      end do
                   end if
                end if
                do k = 1, cea%Nm
                   stcf(m, k) = cea%Stc(j, k)
                end do
                m = m + 1
             end do

             do ii = 1, cea%Nm
                do j = 1, cea%Nr
                   cea%Stc(j, ii) = stcf(j, ii)
                end do
             end do
             cea%Nr = m - 1
          end if
       end do
    end if

    ! FIND TRANSPORT DATA FOR IMPORTANT INTERACTIONS
300 do i = 1, cea%Nm
       cea%Con(i) = 0
       do j = 1, cea%Nm
          cea%Eta(i, j) = 0
       end do
    end do

    outerLoop2: do ir = 1, cea%Ntape
       jtape = cea%transport_properties(ir)%j
       trc = cea%transport_properties(ir)%data

       innerLoop: do k = 1, 2
          do i = 1, cea%Nm
             j = cea%Ind(i)
             if (j == jtape(k)) then
                l = i
                if (k == 2) then
                   kvc = 1

                   do
                      kt = 1
                      if (trc(2, 1, kvc) /= 0) then
                         if (trc(2, 2, kvc) /= 0) then
                            if (cea%Tt > trc(2, 1, kvc)) kt = 2
                            if (trc(2, 3, kvc) /= 0) then
                               if (cea%Tt > trc(2, 2, kvc)) kt = 3
                            end if
                         end if

                         prop = exp(trc(6, kt, kvc) + (trc(5, kt, kvc) / cea%Tt + trc(4, kt, kvc)) / cea%Tt &
                              + trc(3, kt, kvc) * cea%Tln)

                         if (kvc == 2) then
                            cea%Con(l) = prop
                            cycle outerLoop2
                         else
                            cea%Eta(l, m) = prop
                            if (l /= m) cea%Eta(m, l) = cea%Eta(l, m)
                         end if

                      else if (kvc == 2) then
                         cycle outerLoop2
                      end if

                      kvc = 2
                   end do

                else
                   m = i
                   cycle innerLoop
                end if
             end if
          end do

          exit
       end do innerLoop
    end do outerLoop2

    ! MAKE ESTIMATES FOR MISSING DATA
    !
    ! INCLUDES ION CROSS SECTION ESTIMATES
    ! ESTIMATES FOR  E-ION, ION-ION, E-NEUTRAL, ION-NEUTRAL
    ! DEBYE SHIELDING WITH IONIC CUTOFF DISTANCE
    if (cea%Ions) then
       te = cea%Tt / 1000
       ekt = 4.8032d0**2 / (Boltz * te)
       qc = 100 * ekt**2
       xsel = enel / total
       if (xsel < 1.0d-12) xsel = 1.0d-12
       debye = ((22.5d0 / pi) * (R0 / Avgdr * 100) * (te/xsel)) / ekt**3
       ionic = ((810 / (4*pi)) * (R0 / Avgdr * 100d0) * (te/xsel))**(2/3.) / ekt**2
       lambda = sqrt(debye + ionic)
       lambda = max(lambda, 2.71828183d0)
    end if
    do i = 1, cea%Nm
       k = cea%Ind(i)
       cea%Cprr(i) = cea%Cp(k)
       if (.not. (cea%Ions .and. (abs(cea%A(cea%Nlm, k)) == 1) .and. (cea%Eta(i, i) == 0))) then
          if (cea%Eta(i, i) == 0) then
             omega = log(50 * cea%Wmol(i)**4.6 / cea%Tt**1.4)
             omega = max(omega, 1d0)
             cea%Eta(i, i) = cea%Viscns * sqrt(cea%Wmol(i) * cea%Tt) / omega
          end if
          if (cea%Con(i) == 0) cea%Con(i) = cea%Eta(i, i) * R0 * (0.00375d0 + 0.00132d0 * (cea%Cprr(i) - 2.5d0)) / cea%Wmol(i)
       end if
    end do

    do i = 1, cea%Nm
       do j = i, cea%Nm
          ion1 = .false.
          ion2 = .false.
          elc1 = .false.
          elc2 = .false.
          omega = 0
          if (cea%Eta(i, j) == 0) cea%Eta(i, j) = cea%Eta(j, i)
          if (cea%Eta(j, i) == 0) cea%Eta(j, i) = cea%Eta(i, j)
          if (cea%Eta(i, j) == 0) then
             if (cea%Ions) then
                ! ESTIMATE FOR IONS
                k1 = cea%Ind(i)
                k2 = cea%Ind(j)
                if (abs(cea%A(cea%Nlm, k1)) == 1) ion1 = .true.
                if (abs(cea%A(cea%Nlm, k2)) == 1) ion2 = .true.
                if (cea%Wmol(i) < 1) elc1 = .true.
                if (cea%Wmol(j) < 1) elc2 = .true.
                if (ion1 .and. ion2) omega = 1.36d0 * qc * log(lambda)
                if ((ion1 .and. elc2) .or. (ion2 .and. elc1)) omega = 1.29d0 * qc * log(lambda)
                if ((ion1 .and. .not. ion2) .or. (ion2 .and. .not. ion1)) omega = exp(6.776 - 0.4 * cea%Tln)
                if (omega /= 0) then
                   wmred = sqrt(2 * cea%Tt * cea%Wmol(i) * cea%Wmol(j) / (cea%Wmol(i) + cea%Wmol(j)))
                   cea%Eta(i, j) = cea%Viscns * wmred * pi / omega
                   cea%Eta(j, i) = cea%Eta(i, j)
                   if (i == j) then
                      cea%Cprr(i) = cea%Cp(k1)
                      cea%Con(i) = cea%Eta(i, i) * R0 * (0.00375d0 + 0.00132d0 * (cea%Cprr(i) - 2.5d0)) / cea%Wmol(i)
                   end if
                   cycle
                end if
             end if
             ! ESTIMATE FOR UNLIKE INTERACTIONS FROM RIGID SPHERE ANALOGY
             ratio = sqrt(cea%Wmol(j) / cea%Wmol(i))
             cea%Eta(i, j) = 5.656854d0 * cea%Eta(i, i) * sqrt(cea%Wmol(j) / (cea%Wmol(i) + cea%Wmol(j)))
             cea%Eta(i, j) = cea%Eta(i, j) / (1 + sqrt(ratio * cea%Eta(i, i) / cea%Eta(j, j)))**2
             cea%Eta(j, i) = cea%Eta(i, j)
          end if
       end do
    end do

    return
  end subroutine TRANIN



  subroutine TRANP(cea)
    !***********************************************************************
    ! CALCULATES GAS TRANSPORT PROPERTIES
    !
    !   NUMBER OF GASEOUS SPECIES = NM   (MAXIMUM maxTr)
    !   NUMBER OF CHEMICAL REACTIONS = NR (NM - NLM)
    !   ARRAY OF STOICHIOMETRIC COEFFICIENTS = STC
    !***********************************************************************
    use mod_types
    use mod_general
    implicit none

    type(CEA_Core_Problem), intent(inout):: cea

    ! LOCAL VARIABLES
    integer:: i, i1, j, jj, k, m, mm, nlmm, nmm
    real(8):: cpreac, delh(maxTr), gmat(maxMat, maxMat+1), phi(maxTr, maxTr), &
         psi(maxTr, maxTr), reacon, rtpd(maxTr, maxTr), stx(maxTr), &
         stxij(maxTr, maxTr), sumc, sumv, wtmol, xskm(maxTr, maxTr)

    type(CEA_Point), pointer:: p !< current point
    p => cea%points(cea%iOF, cea%ipt)

    call TRANIN(cea)
    ! CALCULATE VISCOSITY AND FROZEN THERMAL CONDUCTIVITY
    nmm = cea%Nm - 1
    do i = 1, cea%Nm
       rtpd(i, i) = 0
       phi(i, i) = 1
       psi(i, i) = 1
    end do
    p%Confro = 0
    p%Vis = 0
    do i = 1, nmm
       i1 = i + 1
       !DIR$ IVDEP
       do j = i1, cea%Nm
          sumc = 2 / (cea%Eta(i, j) * (cea%Wmol(i) + cea%Wmol(j)))
          phi(i, j) = sumc * cea%Wmol(j) * cea%Eta(i, i)
          phi(j, i) = sumc * cea%Wmol(i) * cea%Eta(j, j)
          sumc = (cea%Wmol(i) + cea%Wmol(j))**2
          psi(i, j) = phi(i, j) * (1 + 2.41d0 * (cea%Wmol(i) - cea%Wmol(j)) * (cea%Wmol(i) - 0.142d0 * cea%Wmol(j)) / sumc)
          psi(j, i) = phi(j, i) * (1 + 2.41d0 * (cea%Wmol(j) - cea%Wmol(i)) * (cea%Wmol(j) - 0.142d0 * cea%Wmol(i)) / sumc)
       end do
    end do
    do i = 1, cea%Nm
       sumc = 0
       sumv = 0
       do j = 1, cea%Nm
          sumc = sumc + psi(i, j) * cea%Xs(j)
          sumv = sumv + phi(i, j) * cea%Xs(j)
       end do
       p%Vis = p%Vis + cea%Eta(i, i) * cea%Xs(i) / sumv
       p%Confro = p%Confro + cea%Con(i) * cea%Xs(i) / sumc
    end do
    if (cea%Eql .and. cea%Nr > 0) then
       ! CALCULATE REACTION HEAT CAPACITY AND THERMAL CONDUCTIVITY
       m = cea%Nr + 1
       do i = 1, cea%Nr
          delh(i) = 0
          do k = 1, cea%Lsave
             j = cea%Jcm(k)
             delh(i) = cea%Stc(i, k) * cea%H0(j) + delh(i)
          end do
          nlmm = cea%Lsave + 1
          do k = nlmm, cea%Nm
             j = cea%Ind(k)
             delh(i) = cea%Stc(i, k) * cea%H0(j) + delh(i)
          end do
          cea%G(i, m) = delh(i)
       end do
       do i = 1, maxTr
          do j = 1, maxTr
             if (abs(cea%Stc(i, j)) < 1.0d-6) cea%Stc(i, j) = 0
          end do
       end do
       jj = cea%Nm - 1
       do k = 1, jj
          mm = k + 1
          do m = mm, cea%Nm
             rtpd(k, m) = cea%Wmol(k) * cea%Wmol(m) / (1.1 * cea%Eta(k, m) * (cea%Wmol(k) + cea%Wmol(m)))
             xskm(k, m) = cea%Xs(k) * cea%Xs(m)
             xskm(m, k) = xskm(k, m)
             rtpd(m, k) = rtpd(k, m)
          end do
       end do
       do i = 1, cea%Nr
          do j = i, cea%Nr
             cea%G(i, j) = 0
             gmat(i, j) = 0
          end do
       end do
       do k = 1, jj
          mm = k + 1
          do m = mm, cea%Nm
             if (cea%Xs(k) >= 1.0d-10 .and. cea%Xs(m) >= 1.0d-10) then
                do j = 1, cea%Nr
                   if ((cea%Stc(j, k) == 0) .and. (cea%Stc(j, m) == 0)) stx(j) = 0
                   if ((cea%Stc(j, k) /= 0) .or. (cea%Stc(j, m) /= 0)) &
                        stx(j) = cea%Xs(m) * cea%Stc(j, k) - cea%Xs(k) * cea%Stc(j, m)
                end do
                do i = 1, cea%Nr
                   do j = i, cea%Nr
                      stxij(i, j) = stx(i) * stx(j) / xskm(k, m)
                      cea%G(i, j) = cea%G(i, j) + stxij(i, j)
                      gmat(i, j) = gmat(i, j) + stxij(i, j) * rtpd(k, m)
                   end do
                end do
             end if
          end do
       end do
       m = 1 + cea%Nr
       do i = 1, cea%Nr
          !DIR$ IVDEP
          do j = i, cea%Nr
             cea%G(j, i) = cea%G(i, j)
          end do
          cea%G(i, m) = delh(i)
       end do
       cea%Imat = cea%Nr

       call gauss_elimination(cea%G(1:cea%Imat, 1:cea%Imat+1), cea%X(1:cea%Imat))

       cpreac = 0
       do i = 1, cea%Nr
          cea%G(i, m) = delh(i)
          cpreac = cpreac + cea%R * delh(i) * cea%X(i)
          !DIR$ IVDEP
          do j = i, cea%Nr
             cea%G(i, j) = gmat(i, j)
             cea%G(j, i) = cea%G(i, j)
          end do
       end do

       call gauss_elimination(cea%G(1:cea%Imat, 1:cea%Imat+1), cea%X(1:cea%Imat))

       reacon = 0
       do i = 1, cea%Nr
          reacon = reacon + cea%R*delh(i)*cea%X(i)
       end do
       reacon = 0.6d0 * reacon
    else
       cpreac = 0
       reacon = 0
    end if
    ! CALCULATE OTHER ANSWERS
    p%Cpfro = 0
    wtmol = 0
    do i = 1, cea%Nm
       p%Cpfro = p%Cpfro + cea%Xs(i) * cea%Cprr(i)
       wtmol = wtmol + cea%Xs(i) * cea%Wmol(i)
    end do
    p%Cpfro = p%Cpfro * cea%R / wtmol
    p%Confro = p%Confro / 1000
    if (.not. cea%SIunit) p%Confro = p%Confro / cal_to_J
    p%Vis = p%Vis / 1000
    p%Prfro = p%Vis * p%Cpfro / p%Confro
    if (cea%Eql) then
       cpreac = cpreac / wtmol
       reacon = reacon / 1000
       p%Cpeql = cpreac + p%Cpfro
       p%Coneql = p%Confro + reacon
       p%Preql = p%Vis * p%Cpeql / p%Coneql
    end if
  end subroutine TRANP

end module mod_cea
