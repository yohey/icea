!***********************************************************************
!                     P R O G R A M      C E A 2
!
!             CHEMICAL EQULIBRIUM WITH APPLICATIONS         5/21/04
!***********************************************************************
program main
  use cea
  implicit none
! LOCAL VARIABLES
  character(15):: ensert(20)
  character(200):: infile, ofile
  character(196):: prefix
  logical:: caseOK, ex, readOK
  integer:: i, inc, iof, j, ln, n
  real(8):: xi, xln

  write(*, '(//" ENTER INPUT FILE NAME WITHOUT .inp EXTENSION."/  &
       & "   THE OUTPUT FILES FOR LISTING AND PLOTTING WILL HAVE", / &
       & " THE SAME NAME WITH EXTENSIONS .out AND .plt RESPECTIVELY" &
       & //)')
  read(*, '(a)') prefix
  ln = index(prefix, ' ') - 1
  infile = prefix(1:ln) // '.inp'
  ofile  = prefix(1:ln) // '.out'
  Pfile  = prefix(1:ln) // '.plt'
  inquire(file=infile, exist=ex)
  if (.not. ex) then
     print *, infile, ' DOES NOT EXIST'
     go to 400
  end if
  open(IOINP, file=infile, status='old', form='formatted')
  open(IOOUT, file=ofile, status='unknown', form='formatted')
  open(IOSCH, status='scratch', form='unformatted')
  open(IOTHM, file='thermo.lib', form='unformatted')
  open(IOTRN, file='trans.lib', form='unformatted')
  write(IOOUT, '(/" *******************************************************************************")')
  write(IOOUT, '(/, 9x, "NASA-GLENN CHEMICAL EQUILIBRIUM PROGRAM CEA2,", &
       & " MAY 21, 2004", /19x, "BY  BONNIE MCBRIDE", &
       & " AND SANFORD GORDON", /5x, &
       & " REFS: NASA RP-1311, PART I, 1994", &
       & " AND NASA RP-1311, PART II, 1996")')
  write(IOOUT, '(/" *******************************************************************************")')
  readOK = .true.
  Newr = .false.
100 Iplt = 0
  Nplt = 0
  call INPUT(readOK, caseOK, ensert)
  if (caseOK .and. readOK) then
     do iof = 1, Nof
        if (Oxf(iof) == 0. .and. B0p(1, 1) /= 0.) then
           do i = 1, Nlm
              if (B0p(i, 1) == 0. .or. B0p(i, 2) == 0.) then
                 write(IOOUT, '(/, "OXIDANT NOT PERMITTED WHEN SPECIFYING 100% FUEL(main)")')
                 go to 200
              end if
           end do
        end if
     end do
     if (Ions) then
        if (Elmt(Nlm) /= 'E') then
           Nlm = Nlm + 1
           Elmt(Nlm) = 'E'
           B0p(Nlm, 1) = 0.
           B0p(Nlm, 2) = 0.
        end if
     else if (Elmt(Nlm) == 'E') then
        Nlm = Nlm - 1
     end if
     do n = 1, Nreac
        Jray(n) = 0
     end do
     call SEARCH
     if (Ngc == 0) go to 300
     Newr = .false.
     if (Trnspt) call READTR
! INITIAL ESTIMATES
     Npr = 0
     Gonly = .true.
     Enn = .1d0
     Ennl = -2.3025851
     Sumn = Enn
     xi = Ng
     if (xi == 0.) xi = 1.
     xi = Enn/xi
     xln = log(xi)
     do inc = 1, Nc
        j = Ng + inc
        En(j, 1) = 0.d0
        Enln(j) = 0.d0
     end do
     do j = 1, Ng
        En(j, 1) = xi
        Enln(j) = xln
     end do
     if (Nc /= 0 .and. Nsert /= 0) then
        do 120 i = 1, Nsert
           do j = Ngc, Ngp1, - 1
              if (Prod(j) == ensert(i)) then
                 Npr = Npr + 1
                 Jcond(Npr) = j
                 if (.not. Short) write(IOOUT, '(1X, A16, "INSERTED")') Prod(j)
                 go to 120
              end if
           end do
           write(IOOUT, '(/" WARNING!!!", A16, "NOT FOUND FOR INSERTION")') ensert(i)
120     continue
     end if
     if (Rkt) then
        call ROCKET
     else if (Tp .or. Hp .or. Sp) then
        call THERMP
     else if (Detn) then
        call DETON
     else if (Shock) then
        call SHCK
     end if
     if (Nplt > 0) then
        open(IOPLT, file=Pfile, form='formatted')
        write(IOPLT, '("#", 2x, 20A12)') (Pltvar(j), j = 1, Nplt)
        do i = 1, Iplt
           write(IOPLT, '(1x, 1p, 20E12.4)') (Pltout(i, j), j = 1, Nplt)
        end do
        write(IOPLT, '("#", 2x, 20A12)') (Pltvar(j), j = 1, Nplt)
     end if
  end if
200 if (readOK) go to 100
300 close(IOINP)
  close(IOOUT)
  close(IOSCH)
  close(IOTHM)
  close(IOTRN)
  close(IOPLT)
400 stop
end program



subroutine CPHS
!***********************************************************************
! CALCULATES THERMODYNAMIC PROPERTIES FOR INDIVIDUAL SPECIES
!***********************************************************************
  use cea
  implicit none
! LOCAL VARIABLES
  real(8):: cx(7) = (/0d0, 0d0, 1d0, 0.5d0, 0.6666666666666667d0, 0.75d0, 0.8d0/)
  real(8):: hcx(7) = (/0d0, 0d0, 1d0, 0d0, 0d0, 0d0, 0d0/)
  real(8), save:: scx(7)
  integer, save:: i, ij, j, jj, k

  k = 1
  if (Tt > Tg(2)) k = 2
  if (Tt > Tg(3)) k = 3
  cx(2) = 1.d0 / Tt
  cx(1) = cx(2)**2
  scx(3) = Tln
  scx(2) = -cx(2)
  hcx(2) = Tln * cx(2)
  hcx(1) = -cx(1)
  scx(1) = hcx(1) * .5d0
  do i = 4, 7
     hcx(i) = cx(i) * Tt
     scx(i) = cx(i-1) * Tt
  end do
  do j = 1, Ng
     H0(j) = 0.d0
     S(j) = 0.d0
  end do
  do i = 7, 4, - 1
     do j = 1, Ng
        S(j) = (S(j) + Coef(j, i, k)) * scx(i)
        H0(j) = (H0(j) + Coef(j, i, k)) * hcx(i)
     end do
  end do
  do i = 1, 3
     do j = 1, Ng
        S(j) = S(j) + Coef(j, i, k) * scx(i)
        H0(j) = H0(j) + Coef(j, i, k) * hcx(i)
     end do
  end do
  do j = 1, Ng
     S(j) = S(j) + Coef(j, 9, k)
     H0(j) = H0(j) + Coef(j, 8, k) * cx(2)
  end do
  if (.not. Tp .or. Convg) then
     do j = 1, Ng
        Cp(j) = 0.d0
     end do
     do i = 7, 4, - 1
        do j = 1, Ng
           Cp(j) = (Cp(j) + Coef(j, i, k)) * Tt
        end do
     end do
     do i = 1, 3
        do j = 1, Ng
           Cp(j) = Cp(j) + Coef(j, i, k) * cx(i)
        end do
     end do
  end if
  if (Npr /= 0 .and. k /= 3 .and. Ng /= Ngc) then
     do ij = 1, Npr
        j = Jcond(ij)
        jj = Jcond(ij) - Ng
        Cp(j) = 0.d0
        H0(j) = 0.d0
        S(j) = 0.d0
        do i = 7, 4, - 1
           S(j) = (S(j) + Cft(jj, i)) * scx(i)
           H0(j) = (H0(j) + Cft(jj, i)) * hcx(i)
           Cp(j) = (Cp(j) + Cft(jj, i)) * Tt
        end do
        do i = 1, 3
           S(j) = S(j) + Cft(jj, i) * scx(i)
           H0(j) = H0(j) + Cft(jj, i) * hcx(i)
           Cp(j) = Cp(j) + Cft(jj, i) * cx(i)
        end do
        S(j) = S(j) + Cft(jj, 9)
        H0(j) = H0(j) + Cft(jj, 8) * cx(2)
     end do
  end if
  go to 99999
  entry ALLCON
  do jj = 1, Nc
     j = jj + Ng
     Cp(j) = 0.d0
     H0(j) = 0.d0
     S(j) = 0.d0
     do i = 7, 4, - 1
        S(j) = (S(j) + Cft(jj, i)) * scx(i)
        H0(j) = (H0(j) + Cft(jj, i)) * hcx(i)
        Cp(j) = (Cp(j) + Cft(jj, i)) * Tt
     end do
     do i = 1, 3
        S(j) = S(j) + Cft(jj, i) * scx(i)
        H0(j) = H0(j) + Cft(jj, i) * hcx(i)
        Cp(j) = Cp(j) + Cft(jj, i) * cx(i)
     end do
     S(j) = S(j) + Cft(jj, 9)
     H0(j) = H0(j) + Cft(jj, 8) * cx(2)
  end do
99999 return
end subroutine CPHS



subroutine DETON
!***********************************************************************
! CHAPMAN-JOUGUET DETONATIONS.
!***********************************************************************
  use cea
  implicit none
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
  character(3), save:: unit
  integer, save:: i, ii, iof, itr, j, mdv, mgam, mh, mmach, mp, mson, mt
  real(8), save:: a11, a12, a21, a22, alam, alpha, amm, b1, b2, cpl(Ncol), d, gam, &
       gm1(Ncol), h1(Ncol), p1, pp1, pub(Ncol), rk, rr1, rrho(Ncol), t1, &
       tem, tt1, tub(Ncol), ud, x1, x2


  iof = 0
  Eql = .true.
  if (T(1) == 0.) then
     T(1) = Rtemp(1)
     Nt = 1
  end if
100 Tt = T(1)
  iof = iof + 1
  Oxfl = Oxf(iof)
  call NEWOF
! BEGIN T LOOP.
  do It = 1, Nt
     t1 = T(It)
! BEGIN P LOOP.
     do Ip = 1, Np
        p1 = P(Ip)
        Tt = t1
        Pp = p1
        call HCALC
        if (Tt == 0.) return
        if (Detdbg) call OUT1
        h1(Npt) = Hsub0 * R
        tub(Npt) = t1
        pub(Npt) = p1
        cpl(Npt) = Cpmix * R
        itr = 0
        Tt = 3800.
        pp1 = 15.
        Pp = pp1 * p1
! CALCULATE ENTHALPY FOR INITIAL ESTIMATE OF T2(TT AFTER EQLBRM)
        Hsub0 = h1(Npt) / R + .75 * t1 * pp1 / Wmix
        Tp = .false.
        Hp = .true.
        call EQLBRM
        Hsub0 = h1(Npt) / R
        Hp = .false.
        if (Tt /= 0.) then
           gam = Gammas(Npt)
           tt1 = Tt / t1
           ii = 0
           tem = tt1 - .75 * pp1 / (Cpr(Npt) * Wmix)
           amm = Wm(Npt) / Wmix
           if (Detdbg) write(IOOUT, '(/" T EST.=", F8.2/11X, "P/P1", 17X, "T/T1")') Tt
! LOOP FOR IMPROVING T2/T1 AND P2/P1 INITIAL ESTIMATE.
           do ii = 1, 3
              alpha = amm / tt1
              pp1 = (1. + gam) * (1. + (1. - 4. * gam * alpha / (1. + gam)**2)**.5) &
                   / (2. * gam * alpha)
              rk = pp1 * alpha
              tt1 = tem + .5 * pp1 * gam * (rk**2 - 1.) / (Wmix * Cpr(Npt) * rk)
              if (Detdbg) write(IOOUT, '(I5, 2E20.8)') ii, pp1, tt1
           end do
           Tp = .true.
           Tt = t1 * tt1
           rr1 = pp1 * amm / tt1
! BEGIN MAIN ITERATION LOOP.
110        itr = itr + 1
           Pp = p1 * pp1
           call EQLBRM
           if (Npt == 0) go to 200
           if (Tt /= 0.) then
              gam = Gammas(Npt)
              amm = Wm(Npt) / Wmix
              rr1 = pp1 * amm / tt1
              a11 = 1. / pp1 + gam * rr1 * Dlvpt(Npt)
              a12 = gam * rr1 * Dlvtp(Npt)
              a21 = .5 * gam * (rr1**2 - 1. - Dlvpt(Npt) * (1. + rr1**2)) + Dlvtp(Npt) - 1.
              a22 = -.5 * gam * Dlvtp(Npt) * (rr1**2 + 1.) - Wm(Npt) * Cpr(Npt)
              b1 = 1. / pp1 - 1. + gam * (rr1 - 1.)
              b2 = Wm(Npt) * (Hsum(Npt) - h1(Npt) / R) / Tt - .5 * gam * (rr1**2 - 1.)
              d = a11 * a22 - a12 * a21
              x1 = (a22 * b1 - a12 * b2) / d
              x2 = (a11 * b2 - a21 * b1) / d
              alam = 1.
              tem = x1
              if (tem < 0.) tem = -tem
              if (x2 > tem) tem = x2
              if (-x2 > tem) tem = -x2
              if (tem > 0.4054652) alam = .4054652 / tem
              pp1 = pp1 * exp(x1 * alam)
              tt1 = tt1 * exp(x2 * alam)
              Tt = t1 * tt1
              ud = rr1 * (Rr * gam * Tt/Wm(Npt))**.5
              if (Detdbg) write(IOOUT, '(/" ITER =", I2, 5X, "P/P1 =", E15.8, /7X, "T/T1 =", E15.8, 5X, &
                   & "RHO/RHO1 =", E15.8, /7X, "DEL LN P/P1 =", E15.8, 5X, &
                   & "DEL LN T/T1 =", E15.8)') itr, pp1, tt1, rr1, x1, x2
! CONVERGENCE TEST
              if (itr < 8 .and. tem > 0.5E-04) go to 110
              if (itr < 8) then
                 rrho(Npt) = rr1
                 if (cpl(Npt) == 0.) then
                    gm1(Npt) = 0.
                    Vmoc(Npt) = 0.
                 else
                    gm1(Npt) = cpl(Npt) / (cpl(Npt) - R / Wmix)
                    Vmoc(Npt) = ud / (Rr * gm1(Npt) * t1 / Wmix)**.5
                 end if
              else
                 write(IOOUT, '(/" CONSERVATION EQNS NOT SATISFIED IN 8 ITERATIONS (DETON)")')
                 Npt = Npt - 1
                 Tt = 0.
              end if
              if (Trnspt) call TRANP
              Isv = 0
              if (Ip /= Np .or. It /= Nt .and. Tt /= 0.) then
                 Isv = Npt
                 if (Npt /= Ncol) go to 120
              end if
           end if
! OUTPUT
           write(IOOUT, '(//, 21X, "DETONATION PROPERTIES OF AN IDEAL REACTING GAS")')
           call OUT1
! SET MXX ARRAY FOR PLOTTING PARAMETERS
           mp    = 0
           mt    = 0
           mgam  = 0
           mh    = 0
           mdv   = 0
           mson  = 0
           mmach = 0
           do i = 1, Nplt
              if (index(Pltvar(i)(2:), '1') /= 0) then
                 if (Pltvar(i)(:3) == 'son') then
                    mson = i
                 else if (Pltvar(i)(:3) == 'gam') then
                    mgam = i
                 else if (Pltvar(i)(:1) == 'h') then
                    mh = i
                 else if (Pltvar(i)(:1) == 't') then
                    mt = i
                 else if (Pltvar(i)(:1) == 'p') then
                    mp = i
                 end if
              else if (index(Pltvar(i), 'vel') /= 0) then
                 mdv = i
              else if (index(Pltvar(i), 'mach') /= 0) then
                 mmach = i
              end if
           end do
           write(IOOUT, '(/" UNBURNED GAS"/)')
           Fmt(4) = '13'
           Fmt(5) = ' '
           Fmt(7) = '4,'
           do i = 1, Npt
              if (SIunit) then
                 V(i) = pub(i)
                 unit = 'BAR'
              else
                 V(i) = pub(i) / 1.01325d0
                 unit = 'ATM'
              end if
              if (mp > 0) Pltout(i+Iplt, mp) = V(i)
           end do
           write(IOOUT, Fmt) 'P1, '//unit//'        ', (V(j), j = 1, Npt)
           Fmt(7) = '2,'
           write(IOOUT, Fmt) ft1, (tub(j), j = 1, Npt)
           if (.not. SIunit) write(IOOUT, Fmt) fh1, (h1(j), j = 1, Npt)
           if (SIunit) write(IOOUT, Fmt) fhs1, (h1(j), j = 1, Npt)
           do i = 1, Npt
              V(i) = Wmix
              Sonvel(i) = (Rr * gm1(i) * tub(i) / Wmix)**.5
           end do
           Fmt(7) = '3,'
           write(IOOUT, Fmt) fm1, (V(j), j = 1, Npt)
           Fmt(7) = '4,'
           write(IOOUT, Fmt) fg1, (gm1(j), j = 1, Npt)
           Fmt(7) = '1,'
           write(IOOUT, Fmt) 'SON VEL1,M/SEC ', (Sonvel(j), j = 1, Npt)
           if (Nplt > 0) then
              do i = 1, Npt
                 if (mt > 0)   Pltout(i+Iplt, mt) = tub(i)
                 if (mgam > 0) Pltout(i+Iplt, mgam) = gm1(i)
                 if (mh > 0)   Pltout(i+Iplt, mh) = h1(i)
                 if (mson > 0) Pltout(i+Iplt, mson) = Sonvel(i)
              end do
           end if
           write(IOOUT, '(/" BURNED GAS"/)')
           Fmt(4) = Fmt(6)
           call OUT2
           if (Trnspt) call OUT4
           write(IOOUT, '(/" DETONATION PARAMETERS"/)')
           Fmt(7) = '3,'
           do i = 1, Npt
              V(i) = Ppp(i) / pub(i)
              Pcp(i) = Ttt(i) / tub(i)
              Sonvel(i) = Sonvel(i) * rrho(i)
              if (mmach > 0) Pltout(i+Iplt, mmach) = Vmoc(i)
              if (mdv > 0)   Pltout(i+Iplt, mdv) = Sonvel(i)
           end do
           write(IOOUT, Fmt) fpp1, (V(j), j = 1, Npt)
           write(IOOUT, Fmt) ftt1, (Pcp(j), j = 1, Npt)
           do i = 1, Npt
              V(i) = Wm(i) / Wmix
           end do
           Fmt(7) = '4,'
           write(IOOUT, Fmt) fmm1, (V(j), j = 1, Npt)
           write(IOOUT, Fmt) frr1, (rrho(j), j = 1, Npt)
           write(IOOUT, Fmt) 'DET MACH NUMBER', (Vmoc(j), j = 1, Npt)
           Fmt(7) = '1,'
           write(IOOUT, Fmt) fdv, (Sonvel(j), j = 1, Npt)
           Eql = .true.
           call OUT3
           Iplt = min(Iplt+Npt, 500)
           if (Isv == 0 .and. iof == Nof) go to 200
           if (Np == 1 .and. Nt == 1) go to 100
           write(IOOUT, '(///)')
           Npt = 0
120        Npt = Npt + 1
           if (Isv == 1) Isv = -1
           call SETEN
        end if
     end do
  end do
!     Iplt = min(Iplt + Npt, 500)
  Iplt = min(Iplt + Npt - 1, 500)
  if (iof < Nof) go to 100
200 Tp = .false.
  return
end subroutine DETON



subroutine EFMT(Fone, Aa, Vx)
!***********************************************************************
! WRITE OUTPUT RECORD WITH NUMERICAL VALUES IN SPECIAL EXPONENT FORM.
!***********************************************************************
  use cea
  implicit none
! DUMMY ARGUMENTS
  character(15):: Aa
  character(4):: Fone
  real(8):: Vx(maxMat)
! LOCAL VARIABLES
  character(4):: fmix(5), frmt(8)
  integer, save:: i, j, j1, ne(Ncol)
  real(8), save:: ee, fe, w(Ncol)

  data frmt /'(1H ', ',A15', ',', '9X,', '13(F', '6.4,', 'I2,', '1X))'/
  data fmix /'I3,', '6.4,', 'I2,', '9X,', '5.3,'/

  frmt(6) = fmix(2)
  frmt(7) = fmix(3)
  j1 = 1
  frmt(4) = '1x,'
  if ( Fone == '9X,') then
     j1 = 2
     frmt(4) = fmix(4)
  end if
  do i = j1, Npt
     if (Vx(i) /= 0.) then
        ee = log10(abs(Vx(i)))
        ne(i) = int(ee)
        fe = ne(i)
        if (ee < -.2181E-05 .and. fe /= ee) ne(i) = ne(i) - 1
        if (abs(ne(i)) >= 10) then
           frmt(6) = fmix(5)
           frmt(7) = fmix(1)
        end if
        w(i) = Vx(i) / 10.**ne(i)
     else
        w(i) = 0.
        ne(i) = 0
     end if
  end do
  write(IOOUT, frmt) Aa, (w(j), ne(j), j = j1, Npt)
end subroutine EFMT



subroutine EQLBRM
!***********************************************************************
! CALCULATE EQUILIBRIUM COMPOSITION AND PROPERTIES.
!***********************************************************************
  use cea
  implicit none
! LOCAL VARIABLES
  character(12), save:: ae, cmp(maxEl)
  character(16), save:: amb
  logical, save:: cpcalc, i2many, newcom, reduce
  integer, save:: i, il, ilamb, ilamb1, inc, ipr, iq2, iter, ix, ixsing, iz, j, ja, jb, &
       jbx, jc, jcondi, jcons, jdelg, jex, jj, jkg, jneg, jsw, k, kc, kg, kk, &
       kmat, kneg, l, lc, lcs(maxEl), le, lelim, lk, ll, lncvg, ls, lsing, &
       lz, maxitn, ncvg, njc, nn, numb
  real(8), save:: aa, ambda, ambda1, bigen, bigneg, delg, dlnt, dpie, ensol, esize, &
       gap, gasfrc, pie, pisave(maxMat-2), siz9, sizeg, &
       sum, sum1, szgj, tem, tmelt, tsize, ween, xi, xln, xsize, xx(maxMat)
  real(8):: smalno = 1e-6, smnol = -13.815511

  ixsing = 0
  lsing = 0
  jsw = 0
  jdelg = 0
  maxitn = 50
  ncvg = 0
  lncvg = 3 * Nlm
  reduce = .false.
  siz9 = Size - 9.2103404d0
  tsize = Size
  xsize = Size + 6.90775528d0
  if (Trace /= 0.) then
     maxitn = maxitn + Ngc / 2
     xsize = -log(Trace)
     if (xsize < Size) xsize = Size + .1
  end if
  if (xsize > 80.) xsize = 80.d0
  esize = min(80.d0, xsize + 6.90775528d0)
  jcons = 0
  pie = 0.
  i2many = .false.
  Pderiv = .false.
  Convg = .false.
  numb = 0
  cpcalc = .true.
  if (Tp) cpcalc = .false.
  if (Tt /= 0.d0) then
     if (Npr == 0 .or. (Tt /= T(1) .and. .not. Tp)) go to 400
     k = 1
  else
     Tt = 3800.d0
     if (Npr == 0) go to 400
     k = 1
  end if
100 j = Jcond(k)
  jc = j - Ng
  kg = -Ifz(jc)
  do i = 1, 9
     kg = kg + 1
     kc = jc + kg
     if (Tt <= Temp(2, kc)) then
        if (kg /= 0) then
           Jcond(k) = j + kg
           En(j+kg, Npt) = En(j, Npt)
           En(j, Npt) = 0.
           if (Prod(j) /= Prod(j+kg) .and. .not. Short)  &
                & write(IOOUT, '(" PHASE CHANGE, REPLACE ", A16, "WITH ", A16)') Prod(j), Prod(j+kg)
        end if
        go to 300
     else if (kc >= Nc .or. Ifz(kc+1) <= Ifz(kc)) then
        go to 200
     end if
  end do
200 if (.not. Tp) then
     Tt = Temp(2, kc) - 10.d0
     k = 1
     go to 100
  end if
  write(IOOUT, '(" REMOVE ", A16)') Prod(j)
  En(j, Npt) = 0.d0
  Enln(j) = 0.d0
  Deln(j) = 0.d0
  do i = k, Npr
     Jcond(i) = Jcond(i+1)
  end do
  Npr = Npr - 1
300 k = k + 1
  if (k <= Npr) go to 100
400 Tln = log(Tt)
  if (Vol) Pp = Rr * Enn * Tt / Vv
  call CPHS
  Tm = log(Pp / Enn)
  le = Nlm
  if (Lsave /= 0 .and. Nlm /= Lsave) then
     tem = exp(-tsize)
     do i = Lsave + 1, Nlm
        do j = 1, Ng
           if (A(i, j) /= 0.) then
              En(j, Npt) = tem
              Enln(j) = -tsize
           end if
        end do
     end do
  end if
  ls = Nlm
  lelim = 0
  lz = ls
  if (Ions) lz = ls - 1
  if (Npt == 1 .and. .not. Shock .and. .not. Short) then
     write(IOOUT, '(/" POINT ITN", 6X, "T", 10X, 4(A4, 8X)/(18X, 5(A4, 8X)))') (Elmt(i), i = 1, Nlm)
  end if
  if (Debug(Npt)) then
     do i = 1, Nlm
        cmp(i) = Elmt(i)
     end do
  end if
! BEGIN ITERATION
500 if (cpcalc) then
     Cpsum = 0.d0
     do j = 1, Ng
        Cpsum = Cpsum + En(j, Npt) * Cp(j)
     end do
     if (Npr /= 0) then
        do k = 1, Npr
           j = Jcond(k)
           Cpsum = Cpsum + En(j, Npt) * Cp(j)
        end do
        cpcalc = .false.
     end if
  end if
  numb = numb + 1
  call MATRIX
  iq2 = Iq1 + 1
  if (Convg) Imat = Imat - 1
  if (Debug(Npt)) then
     if (.not. Convg) then
        write(IOOUT, '(/" ITERATION", I3, 6X, "MATRIX ")') numb
     else
        if (.not. Pderiv) write(IOOUT, '(/" T DERIV MATRIX")')
        if (Pderiv) write(IOOUT, '(/" P DERIV MATRIX")')
     end if
     kmat = Imat + 1
     do i = 1, Imat
        write(IOOUT, '(3X, 5E15.6)') (G(i, k), k = 1, kmat)
     end do
  end if
  Msing = 0
  call GAUSS
  if (Msing == 0) then
     if (Debug(Npt)) then
        write(IOOUT, '(/" SOLUTION VECTOR", /, 6x, 5A15/8X, 5A15)') (cmp(k), k = 1, le)
        write(IOOUT, '(3X, 5E15.6)') (X(i), i = 1, Imat)
     end if
     if (.not. Convg) then
! OBTAIN CORRECTIONS TO THE ESTIMATES
        if (Vol) X(iq2) = X(Iq1)
        if (Tp) X(iq2) = 0.
        dlnt = X(iq2)
        sum = X(Iq1)
        if (Vol) then
           X(Iq1) = 0.
           sum = -dlnt
        end if
        do 520 j = 1, Ng
           if (lelim /= 0) then
              Deln(j) = 0.
              do i = lelim, ls
                 if (A(i, j) /= 0.) go to 520
              end do
           end if
           Deln(j) = -Mu(j) + H0(j) * dlnt + sum
           do k = 1, Nlm
              Deln(j) = Deln(j) + A(k, j) * X(k)
           end do
           if (pie /= 0.) Deln(j) = Deln(j) + A(ls, j) * pie
520     continue
        if (Npr /= 0) then
           do k = 1, Npr
              j = Jcond(k)
              kk = Nlm + k
              Deln(j) = X(kk)
           end do
        end if
! CALCULATE CONTROL FACTOR, AMBDA
        ambda = 1.d0
        ambda1 = 1.d0
        ilamb = 0
        ilamb1 = 0
        sum = max(abs(X(Iq1)), abs(dlnt))
        sum = sum*5.
        do j = 1, Ng
           if (Deln(j) > 0.) then
              if ((Enln(j) - Ennl + Size) <= 0.) then
                 sum1 = abs(Deln(j)-X(Iq1))
                 if (sum1 >= siz9) then
                    sum1 = abs(-9.2103404d0 - Enln(j) + Ennl) / sum1
                    if (sum1 < ambda1) then
                       ambda1 = sum1
                       ilamb1 = j
                    end if
                 end if
              else if (Deln(j) > sum) then
                 sum = Deln(j)
                 ilamb = j
              end if
           end if
        end do
        if (sum > 2.d0) ambda = 2.d0/sum
        if (ambda1 <= ambda) then
           ambda = ambda1
           ilamb = ilamb1
        end if
        if (Debug(Npt)) then
! INTERMEDIATE OUTPUT
           write(IOOUT, '(/" T=", E15.8, " ENN=", E15.8, " ENNL=", E15.8, " PP=", E15.8, &
                & /" LN P/N=", E15.8, " AMBDA=", E15.8)') Tt, Enn, Ennl, Pp, Tm, ambda
           if (ambda /= 1.d0) then
              amb = 'ENN'
              if (abs(X(iq2)) > abs(X(Iq1))) amb = 'TEMP'
              if (ilamb /= 0) amb = Prod(ilamb)
              write(IOOUT, '(/" AMBDA SET BY ", A16)') amb
           end if
           if (Vol) write(IOOUT, '(" VOLUME=", E15.8, "CC/G")') Vv * .001d0
           write(IOOUT, '(/24X, "Nj", 12X, "LN Nj", 8X, "DEL LN Nj", 6X, "H0j/RT", /, 41X, &
                & "S0j/R", 10X, " G0j/RT", 8X, " Gj/RT")')
           do j = 1, Ngc
              write(IOOUT, '(1X, A16, 4E15.6, /35x, 3E15.6)') Prod(j), En(j, Npt), Enln(j), Deln(j), &
                   H0(j), S(j), H0(j) - S(j), Mu(j)
           end do
        end if
! APPLY CORRECTIONS TO ESTIMATES
        Totn(Npt) = 0.d0
        do j = 1, Ng
           Enln(j) = Enln(j) + ambda * Deln(j)
        end do
        do 540 j = 1, Ng
           En(j, Npt) = 0.
           if (lelim /= 0) then
              do i = lelim, ls
                 if (A(i, j) /= 0.) go to 540
              end do
           end if
           if ((Enln(j) - Ennl + tsize) > 0.) then
              En(j, Npt) = exp(Enln(j))
              Totn(Npt) = Totn(Npt) + En(j, Npt)
           end if
540     continue
        if (Ions .and. Elmt(Nlm) == 'E') then
           do j = 1, Ng
              if (A(ls, j) /= 0. .and. En(j, Npt) == 0.) then
                 if ((Enln(j)-Ennl+esize) > 0.) then
                    En(j, Npt) = exp(Enln(j))
                    Totn(Npt) = Totn(Npt) + En(j, Npt)
                 end if
              end if
           end do
        end if
        Sumn = Totn(Npt)
        if (Npr /= 0) then
           do k = 1, Npr
              j = Jcond(k)
              En(j, Npt) = En(j, Npt) + ambda * Deln(j)
              Totn(Npt) = Totn(Npt) + En(j, Npt)
           end do
        end if
        if (.not. Tp) then
           Tln = Tln + ambda * dlnt
           Tt = exp(Tln)
           cpcalc = .true.
           call CPHS
        end if
        if (Vol) then
           Enn = Sumn
           Ennl = log(Enn)
           if (Vol) Pp = Rr * Tt * Enn / Vv
        else
           Ennl = Ennl + ambda * X(Iq1)
           Enn = exp(Ennl)
        end if
        Tm = log(Pp / Enn)
        if (Elmt(Nlm) == 'E') then
! CHECK ON REMOVING IONS
           do j = 1, Ngc
              if (A(Nlm, j) /= 0.) then
                 if (En(j, Npt) > 0.) go to 560
              end if
           end do
           pie = X(Nlm)
           lelim = Nlm
           Nlm = Nlm - 1
           go to 500
        end if
! TEST FOR CONVERGENCE
560     if (numb > maxitn) then
           write(IOOUT, '(/, I4, " ITERATIONS DID NOT SATISFY CONVERGENCE", /, 15x, &
                & " REQUIREMENTS FOR THE POINT", I5, " (EQLBRM)")') maxitn, Npt
           if (Nc == 0 .or. i2many) go to 1500
           i2many = .true.
           if (.not. Hp .or. Npt /= 1 .or. Tt > 100.) then
              if (Npr /= 1 .or. Enn > 1.E-4) go to 1500
! HIGH TEMPERATURE, INCLUDED CONDENSED CONDITION
              write(IOOUT, '(/" TRY REMOVING CONDENSED SPECIES (EQLBRM)")')
              Enn = .1
              Ennl = -2.3025851
              Sumn = Enn
              xi = Ng
              xi = Enn/xi
              xln = log(xi)
              do j = 1, Ng
                 En(j, Npt) = xi
                 Enln(j) = xln
              end do
              j = Jcond(1)
              k = 1
              go to 1000
           else
              write(IOOUT, '(/" LOW TEMPERATURE IMPLIES A CONDENSED SPECIES SHOULD HA", &
                   & "VE BEEN INSERTED,", &
                   & /" RESTART WITH insert DATASET (EQLBRM)")')
              go to 1500
           end if
        else
           sum = (X(Iq1) * Enn / Totn(Npt))
           if (abs(sum) > 0.5E-5) go to 500
           do j = 1, Ng
              if (abs(Deln(j)) * En(j, Npt) / Totn(Npt) > 0.5d-5) go to 500
           end do
           if (abs(dlnt) > 1.d-04) go to 500
           if (Npr /= 0) then
              do k = 1, Npr
                 j = Jcond(k)
                 if (abs(Deln(j)/Totn(Npt)) > 0.5d-5) go to 500
                 if (En(j, Npt) < 0.) go to 700
              end do
           end if
           le = Nlm
           do i = 1, Nlm
              if (abs(B0(i)) >= 1.d-06) then
                 sum = 0.
                 do j = 1, Ngc
                    sum = sum + En(j, Npt) * A(i, j)
                 end do
                 if (abs(B0(i)-sum) > Bcheck) go to 500
              end if
           end do
           if (Trace /= 0.) then
              tsize = xsize
              tem = 1.
              if (numb /= 1) then
                 lk = lz
                 if (Nlm < lz) lk = Nlm
                 do i = 1, lk
                    if (i /= lsing) then
                       tem = 0.
                       if (X(i) /= 0.) then
                          tem = abs((pisave(i) - X(i)) / X(i))
                          if (tem > .001) go to 565
                       end if
                    end if
                 end do
              end if
565           do i = 1, Nlm
                 pisave(i) = X(i)
              end do
              if (tem > .001) go to 500
              if (Ions) then
! CHECK ON ELECTRON BALANCE
                 iter = 1
                 if (pie /= 0.) then
                    le = Nlm + 1
                    X(le) = pie
                 end if
566              sum1 = 0.
                 sum = 0.
                 pie = X(le)
                 do j = 1, Ng
                    if (A(ls, j) /= 0.) then
                       En(j, Npt) = 0.
                       tem = 0.
                       if (Enln(j) > -87.) tem = exp(Enln(j))
                       if ((Enln(j)-Ennl+tsize) > 0. .and. Elmt(Nlm) == 'E') then
                          pie = 0.
                          En(j, Npt) = tem
                       end if
                       aa = A(ls, j) * tem
                       sum = sum + aa
                       sum1 = sum1 + aa * A(ls, j)
                    end if
                 end do
                 if (sum1 /= 0.) then
                    dpie = -sum / sum1
                    do j = 1, Ng
                       if (A(ls, j) /= 0.) Enln(j) = Enln(j) + A(ls, j) * dpie
                    end do
                    if (Debug(Npt)) write(IOOUT, '(/" ELECTRON BALANCE ITER NO. =", I4, "  DELTA PI =", E14.7)') iter, dpie
                    if (abs(dpie) > .0001) then
                       X(le) = X(le) + dpie
                       iter = iter + 1
                       if (iter <= 80) go to 566
                       write(IOOUT, '(/" DID NOT CONVERGE ON ELECTRON BALANCE (EQLBRM)")')
                       go to 1500
                    else if (Elmt(Nlm) == 'E' .and. pie /= 0.) then
                       Nlm = Nlm - 1
                       newcom = .true.
                    end if
                 end if
              end if
           end if
        end if
     else if (.not. Pderiv) then
! TEMPERATURE DERIVATIVES--CONVG=T, PDERIV=F
        Dlvtp(Npt) = 1. - X(Iq1)
        Cpr(Npt) = G(iq2, iq2)
        do j = 1, Iq1
           Cpr(Npt) = Cpr(Npt) - G(iq2, j) * X(j)
        end do
! PRESSURE DERIVATIVE--CONVG=T, PDERIV=T
        Pderiv = .true.
        go to 500
     else
        Dlvpt(Npt) = -1. + X(Iq1)
        if (Jliq == 0) then
           Gammas(Npt) = -1. / (Dlvpt(Npt) + (Dlvtp(Npt)**2) * Enn / Cpr(Npt))
        else
           En(Jsol, Npt) = ensol
           Hsum(Npt) = Hsum(Npt) + En(Jliq, Npt) * (H0(Jliq) - H0(Jsol))
           Gammas(Npt) = -1. / Dlvpt(Npt)
           Npr = Npr + 1
           Jcond(Npr) = Jliq
        end if
        go to 1400
     end if
! SINGULAR MATRIX
  else
     if (Convg) then
        write(IOOUT, '(/" DERIVATIVE MATRIX SINGULAR (EQLBRM)")')
        Dlvpt(Npt) = -1.
        Dlvtp(Npt) = 1.
        Cpr(Npt) = Cpsum
        Gammas(Npt) = -1. / (Dlvpt(Npt) + Dlvtp(Npt)**2 * Enn / Cpr(Npt))
        go to 1400
     else
        write(IOOUT, '(/" SINGULAR MATRIX, ITERATION", I3, "  VARIABLE", I3, "(EQLBRM)")') numb, Msing
        lsing = Msing
        ixsing = ixsing + 1
        if (ixsing <= 8) then
           xsize = 80.
           tsize = xsize
           if (Msing > Nlm .and. numb < 1 .and. Npr > 1 .and. jdelg > 0) then
              ween = 1000.
              j = 0
              do 570 i = 1, Npr
                 jcondi = Jcond(i)
                 if (jcondi /= jdelg) then
                    do ll = 1, Nlm
                       if (A(ll, jdelg) /= 0 .and. A(ll, jcondi) /= 0.) then
                          if (En(jcondi, Npt) <= ween) then
                             ween = En(jcondi, Npt)
                             j = jcondi
                             k = i
                          end if
                          go to 570
                       end if
                    end do
                 end if
570           continue
              if (j > 0) then
                 write(IOOUT, '(/" TRY REMOVING CONDENSED SPECIES (EQLBRM)")')
                 go to 1000
              end if
           else if (.not. Hp .or. Npt /= 1 .or. Nc == 0 .or. Tt > 100.) then
              if (ixsing >= 3) then
                 if (Msing < Iq1) then
                    if (reduce .and. Msing <= Nlm) then
                       if (Nlm < lelim) go to 1500
                       write(IOOUT, '(/" WARNING!! POINT", I3, &
                            & " USES A REDUCED SET OF COMPONENTS", / &
                            & " SPECIES CONTAINING THE ELIMINATED COMPONENT ARE OMITTED.", &
                            & / &
                            & " IT MAY BE NECESSARY TO RERUN WITH INSERTED CONDENSED SPECIES", &
                            & /" CONTAINING COMPONENT ", A8, "(EQLBRM)")') Npt, Elmt(Nlm)
                       Nlm = Nlm - 1
                       go to 500
                    else if (Msing <= Nlm) then
! FIND NEW COMPONENTS
                       if (.not. Ions) go to 1100
                       if (Elmt(Nlm) /= 'E') go to 1100
                       do j = 1, Ng
                          if (A(Nlm, j) /= 0.) En(j, Npt) = 0.d0
                       end do
                       pie = X(Nlm)
                       Nlm = Nlm - 1
                       if (Msing > Nlm) go to 500
                       go to 1100
                    else
! REMOVE CONDENSED SPECIES TO CORRECT SINGULARITY
                       k = Msing - Nlm
                       j = Jcond(k)
                       if (j /= jcons) then
                          jcons = j
                          go to 1000
                       end if
                    end if
                 end if
              end if
              do 575 jj = 1, Ng
                 if (Ions) then
                    if (Elmt(Nlm) /= 'E') then
                       if (A(ls, jj) /= 0.) go to 575
                    end if
                 end if
                 if (En(jj, Npt) == 0.) then
                    En(jj, Npt) = smalno
                    Enln(jj) = smnol
                 end if
575           continue
              go to 500
           else
              write(IOOUT, '(/" LOW TEMPERATURE IMPLIES A CONDENSED SPECIES SHOULD HA", &
                   & "VE BEEN INSERTED,", &
                   & /" RESTART WITH insert DATASET (EQLBRM)")')
           end if
        end if
     end if
     go to 1500
  end if
! CALCULATE ENTROPY, CHECK ON DELTA S FOR SP PROBLEMS
600 Ssum(Npt) = 0.
  do j = 1, Ng
     Ssum(Npt) = Ssum(Npt) + En(j, Npt) * (S(j) - Enln(j) - Tm)
  end do
  if (Npr > 0) then
     do k = 1, Npr
        j = Jcond(k)
        Ssum(Npt) = Ssum(Npt) + En(j, Npt) * S(j)
     end do
  end if
  if (.not. Sp) then
     Convg = .true.
  else
     tem = Ssum(Npt) - S0
     if (abs(tem) > .0005) go to 500
     if (Debug(Npt)) write(IOOUT, '(/" DELTA S/R =", E15.8)') tem
     Convg = .true.
  end if
! CONVERGENCE TESTS ARE SATISFIED, TEST CONDENSED SPECIES.
700 ncvg = ncvg + 1
  if (ncvg > lncvg) then
! ERROR, SET TT=0
     write(IOOUT, '(/, I3, " CONVERGENCES FAILED TO ESTABLISH SET OF CONDENSED", " SPECIES (EQLBRM)")') lncvg
     go to 1500
  else
     if (.not. Shock) then
        do il = 1, le
           xx(il) = X(il)
        end do
        if (.not. Short) then
           if (newcom) write(IOOUT, '(/" POINT ITN", 6X, "T", 10X, 4A12/(18X, 5A12))') (cmp(k), k = 1, le)
           write(IOOUT, '(I4, I5, 5F12.3, /(12X, 5F12.3))') Npt, numb, Tt, (xx(il), il = 1, le)
        end if
        if (.not. Tp .and. Npr == 0 .and. Tt <= Tg(1) * .2d0) then
           write(IOOUT, '(/" LOW TEMPERATURE IMPLIES A CONDENSED SPECIES SHOULD HA", &
                & "VE BEEN INSERTED,", &
                & /" RESTART WITH insert DATASET (EQLBRM)")')
           go to 1500
        end if
        newcom = .false.
     end if
     if (Npr /= 0) then
        bigneg = 0.
        jneg = 0
        do k = 1, Npr
           j = Jcond(k)
           if (En(j, Npt)*Cp(j) <= bigneg) then
              bigneg = En(j, Npt) * Cp(j)
              jneg = j
              kneg = k
           end if
        end do
        if (jneg /= 0) then
           j = jneg
           k = kneg
           if (j == Jsol .or. j == Jliq) then
              Jsol = 0
              Jliq = 0
           end if
           go to 1000
        end if
     end if
     if (Ngc /= Ng .or. Tp) then
        Ng = Ngc
        call CPHS
        Ng = Ngp1 - 1
        cpcalc = .true.
        if (Ngc == Ng) go to 750
        call ALLCON
        if (Npr /= 0 .and. .not. Tp) then
           gap = 50.
           do 710 ipr = 1, Npr
              j = Jcond(ipr)
              if (j /= Jsol .and. j /= Jliq) then
                 inc = j - Ng
                 kg = -Ifz(inc)
                 do iz = 1, 20
                    kg = kg + 1
                    kc = inc + kg
                    if (Tt <= Temp(2, kc)) then
                       if (kg /= 0) then
                          jkg = j + kg
                          if (abs(kg) > 1 .or. Prod(j) == Prod(jkg)) &
                               go to 740
                          if (jkg == jsw) go to 720
                          if (Tt < Temp(1, inc) - gap .or. Tt > Temp(2, inc) + gap) go to 740
                          go to 720
                       end if
                       go to 710
                    else if (Ifz(kc+1) <= Ifz(kc)) then
                       go to 710
                    end if
                 end do
                 if (Tt > Temp(2, kc) * 1.2d0) go to 1000
              end if
710        continue
        end if
        sizeg = 0.
        szgj = 0.
        do inc = 1, Nc
           j = inc + Ng
           if (Debug(Npt)) write(IOOUT, '(/1X, A15, 2F10.3, 3X, E15.7)') Prod(j), Temp(1, inc), Temp(2, inc), En(j, Npt)
           if (En(j, Npt) <= 0.) then
              if (Tt > Temp(1, inc) .or. Temp(1, inc) == Tg(1)) then
                 if (Tt <= Temp(2, inc)) then
                    sum = 0.
                    do i = 1, Nlm
                       sum = sum + A(i, j) * X(i)
                    end do
                    delg = (H0(j) - S(j) - sum) / Mw(j)
                    if (delg < sizeg .and. delg < 0.) then
                       if (j /= jcons) then
                          sizeg = delg
                          jdelg = j
                       else
                          szgj = delg
                       end if
                       ipr = ipr - 1
                    end if
                    if (Debug(Npt)) write(IOOUT, '(" [G0j-SUM(Aij*PIi)]/Mj =", E15.7, 9X, "MAX NEG DELTA G =", E15.7)') delg, sizeg
                 end if
              end if
           end if
        end do
        if (sizeg == 0. .and. szgj == 0.) go to 750
        if (sizeg /= 0.) then
           j = jdelg
           go to 800
        else
           write(IOOUT, '(/" REINSERTION OF ", A16, " LIKELY TO CAUSE SINGULARITY, ", "(EQLBRM)")') Prod(jcons)
           go to 1500
        end if
720     kk = max(0, kg)
        tmelt = Temp(kk+1, inc)
        Tt = tmelt
        Tln = log(Tt)
        Jsol = min(j, jkg)
        Jliq = Jsol + 1
        En(jkg, Npt) = .5d0 * En(j, Npt)
        En(j, Npt) = En(jkg, Npt)
        j = jkg
        go to 800
! WRONG PHASE INCLUDED FOR T INTERVAL, SWITCH EN
740     En(jkg, Npt) = En(j, Npt)
        Jcond(ipr) = jkg
        En(j, Npt) = 0.
        jsw = j
        if (Prod(j) /= Prod(jkg) .and. .not. Short) &
             write(IOOUT, '(" PHASE CHANGE, REPLACE ", A16, "WITH ", A16)') Prod(j), Prod(jkg)
        j = jkg
        go to 900
     end if
! CONVERGED WITH NO CONDENSED CHANGES.  IF BOTH SOLID & LIQ PRESENT, 
! TEMPORARILY REMOVE LIQ TO PREVENT SINGULAR DERIVATIVE MATRIX.
750  Sumn = Enn
     if (Jsol /= 0) then
        ensol = En(Jsol, Npt)
        En(Jsol, Npt) = En(Jsol, Npt) + En(Jliq, Npt)
        Dlvtp(Npt) = 0.
        Cpr(Npt) = 0.
        Gammas(Npt) = 0.
        Pderiv = .true.
        do k = 1, Npr
           if (Jcond(k) == Jliq) go to 760
        end do
760     do i = k, Npr
           Jcond(i) = Jcond(i+1)
        end do
        Npr = Npr - 1
     end if
     go to 500
  end if
! ADD CONDENSED SPECIES
800 Npr = Npr + 1
  i = Npr
  do ix = 2, Npr
     Jcond(i) = Jcond(i-1)
     i = i - 1
  end do
  Jcond(1) = j
  if (.not. Short) write(IOOUT, '(" ADD ", A16)') Prod(j)
900 inc = j - Ng
  Convg = .false.
  if (Tp) cpcalc = .false.
  numb = -1
  go to 500
! REMOVE CONDENSED SPECIES
1000 En(j, Npt) = 0.d0
  Deln(j) = 0.d0
  Enln(j) = 0.d0
  do i = k, Npr
     Jcond(i) = Jcond(i+1)
  end do
  if (.not. Short) write(IOOUT, '(" REMOVE ", A16)') Prod(j)
  Npr = Npr - 1
  do i = 1, Nlm
     if (cmp(i) == Prod(j)) then
        numb = -1
        Convg = .false.
        if (Tp) cpcalc = .false.
        go to 1100
     end if
  end do
  go to 900
1100 newcom = .false.
  nn = Nlm
  if (Elmt(Nlm) == 'E') nn = Nlm - 1
! FIND ORDER OF SPECIES FOR COMPONENTS - BIGGEST TO SMALLEST
  njc = 0
  do lc = 1, nn
     lcs(lc) = 0
  end do
1200 bigen = -1.d-35
  do j = 1, Ng
     if (En(j, Npt) > bigen) then
        if (.not. Ions .or. A(ls, j) == 0.) then
           bigen = En(j, Npt)
           jbx = j
        end if
     end if
  end do
  if (bigen > 0.) then
     do 1250 lc = 1, nn
        if (jbx == 0) jbx = Jx(lc)
        if (A(lc, jbx) > smalno) then
           if (njc /= 0) then
              do 1205 i = 1, njc
                 l = lcs(i)
                 if (l == lc) go to 1250
                 if (l == 0) go to 1210
                 j = Jcm(l)
                 do l = 1, nn
                    if (A(l, jbx) /= A(l, j)) go to 1205
                 end do
                 go to 1250
1205          continue
           end if
1210       do i = 1, nn
              if (i /= lc) then
                 jex = Jx(i)
                 if (abs(A(lc, jbx) * A(i, jex) - A(lc, jex) * A(i, jbx)) <= smalno) go to 1250
              end if
           end do
           njc = njc + 1
           if (jbx /= Jcm(lc)) newcom = .true.
           Jcm(lc) = jbx
           lcs(njc) = lc
           go to 1300
        end if
1250 continue
1300 En(jbx, Npt) = -En(jbx, Npt)
     if (njc < nn) go to 1200
  end if
  do j = 1, Ng
     En(j, Npt) = abs(En(j, Npt))
  end do
  if (newcom) then
! SWITCH COMPONENTS
     do lc = 1, nn
        jb = Jcm(lc)
        if (A(lc, jb) == 0.) then
           jb = Jx(lc)
           Jcm(lc) = jb
        end if
        tem = A(lc, jb)
        if (tem /= 0.) then
           pisave(lc) = H0(jb) - S(jb)
           if (jb <= Ng) pisave(lc) = pisave(lc) + Enln(jb) + Tm
           cmp(lc) = trim(Prod(jb))
! CALCULATE NEW COEFFICIENTS
           if (tem /= 1.) then
              B0(lc) = B0(lc) / tem
              B0p(lc, 1) = B0p(lc, 1) / tem
              B0p(lc, 2) = B0p(lc, 2) / tem
              do j = 1, Nspx
                 A(lc, j) = A(lc, j) / tem
              end do
           end if
           do i = 1, nn
              if (A(i, jb) /= 0. .and. i /= lc) then
                 tem = A(i, jb)
                 do j = 1, Nspx
                    A(i, j) = A(i, j) - A(lc, j) * tem
                    if (abs(A(i, j)) < 1.E-5) A(i, j) = 0.
                 end do
                 B0(i) = B0(i) - B0(lc) * tem
                 B0p(i, 1) = B0p(i, 1) - B0p(lc, 1) * tem
                 B0p(i, 2) = B0p(i, 2) - B0p(lc, 2) * tem
              end if
           end do
        end if
     end do
     if (Debug(Npt)) then
        write(IOOUT, '(/" NEW COMPONENTS")')
        write(IOOUT, '(/2x, 6A12)') (cmp(k), k = 1, nn)
     end if
  end if
  if (Msing /= 0) then
! SWITCH ORDER OF MSING AND NLM COMPONENTS
     reduce = .true.
     lelim = Nlm
     lsing = Nlm
     if (Msing /= Nlm) then
        do j = 1, Nspx
           aa = A(Msing, j)
           A(Msing, j) = A(Nlm, j)
           A(Nlm, j) = aa
        end do
        ja = Jcm(Msing)
        Jcm(Msing) = Jcm(Nlm)
        Jcm(Nlm) = ja
        ae = cmp(Msing)
        cmp(Msing) = cmp(Nlm)
        cmp(Nlm) = ae
        ae = Elmt(Msing)
        Elmt(Msing) = Elmt(Nlm)
        Elmt(Nlm) = trim(ae)
        ja = Jx(Msing)
        Jx(Msing) = Jx(Nlm)
        Jx(Nlm) = ja
        aa = Atwt(Msing)
        Atwt(Msing) = Atwt(Nlm)
        Atwt(Nlm) = aa
        aa = B0(Msing)
        B0(Msing) = B0(Nlm)
        B0(Nlm) = aa
        aa = pisave(Msing)
        pisave(Msing) = pisave(Nlm)
        pisave(Nlm) = aa
        do i = 1, 2
           aa = B0p(Msing, i)
           B0p(Msing, i) = B0p(Nlm, i)
           B0p(Nlm, i) = aa
        end do
     end if
  else if (.not. newcom .and. Trace == 0.) then
     go to 600
  end if
  Msing = 0
  tsize = xsize
  go to 500
1400 Ttt(Npt) = Tt
  Ppp(Npt) = Pp
  Vlm(Npt) = Rr * Enn * Tt / Pp
  Hsum(Npt) = Hsum(Npt) * Tt
  Wm(Npt) = 1. / Enn
  gasfrc = Enn/Totn(Npt)
  if (gasfrc < .0001) write(IOOUT, '(/" WARNING!  RESULTS MAY BE WRONG FOR POINT", I3, " DUE TO", &
       & /" LOW MOLE FRACTION OF GASES (", E15.8, ") (EQLBRM)")') Npt, gasfrc
  if (Trace /= 0.) then
     do 1450 j = 1, Ng
        if (lelim /= 0) then
           do i = lelim, ls
              if (A(i, j) /= 0.) go to 1450
           end do
        end if
        if (Enln(j) > -87.) En(j, Npt) = exp(Enln(j))
1450 continue
  end if
  if (Debug(Npt)) write(IOOUT, '(/" POINT=", I3, 3X, "P=", E13.6, 3X, "T=", E13.6, /3X, "H/R=", &
       & E13.6, 3X, "S/R=", E13.6, /3X, "M=", E13.6, 3X, "CP/R=", E13.6, 3X, &
       & "DLVPT=", E13.6, /3X, "DLVTP=", E13.6, 3X, "GAMMA(S)=", E13.6, 3X, &
       & "V=", E13.6)') Npt, Pp, Tt, Hsum(Npt), &
       Ssum(Npt), Wm(Npt), Cpr(Npt), Dlvpt(Npt), &
       Dlvtp(Npt), Gammas(Npt), Vlm(Npt)
  if (Tt >= Tg(1) .and. Tt <= Tg(4)) go to 1600
  if (Shock) go to 1600
  write(IOOUT, '(" THE TEMPERATURE=", E12.4, " IS OUT OF RANGE FOR POINT", I5, "(EQLBRM)")') Tt, Npt
  if (Tt >= Tg(1) * .8d0 .and. Tt <= Tg(4) * 1.1d0) go to 1600
  Npt = Npt + 1
1500 Tt = 0.
  Npt = Npt - 1
  write(IOOUT, '(/" CALCULATIONS STOPPED AFTER POINT", I3, "(EQLBRM)")') Npt
1600 Lsave = Nlm
  Nlm = ls
  if (Npr > 0) Gonly = .false.
  return
end subroutine



subroutine FROZEN
!***********************************************************************
! CALCULATE PROPERTIES WITH FROZEN COMPOSITION AT ASSIGNED ENTROPY
! AND PRESSURE.  CALLED FROM ROCKET.
!***********************************************************************
  use cea
  implicit none
! LOCAL VARIABLES
  integer, save:: i, inc, iter, j, k, nnn
  real(8), save:: dlnt, dlpm

  Convg = .false.
  Tln = log(Tt)
  dlpm = log(Pp * Wm(Nfz))
  nnn = Npt
  Npt = Nfz
  do j = 1, Ng
     if (En(j, Nfz) /= 0.d0) Deln(j) = -(log(En(j, Nfz)) + dlpm)
  end do
  do iter = 1, 8
     Ssum(nnn) = 0.d0
     Cpsum = 0.d0
     call CPHS
     do j = 1, Ng
        Cpsum = Cpsum + En(j, Nfz) * Cp(j)
        Ssum(nnn) = Ssum(nnn) + En(j, Nfz) * (S(j) + Deln(j))
     end do
     if (Npr /= 0) then
        do k = 1, Npr
           j = Jcond(k)
           Cpsum = Cpsum + En(j, Nfz) * Cp(j)
           Ssum(nnn) = Ssum(nnn) + En(j, Nfz) * S(j)
        end do
     end if
     if (Convg) then
        Npt = nnn
        Hsum(Npt) = 0.d0
        do j = 1, Ngc
           Hsum(Npt) = Hsum(Npt) + En(j, Nfz) * H0(j)
        end do
        Hsum(Npt) = Hsum(Npt) * Tt
        Ttt(Npt) = Tt
        Gammas(Npt) = Cpsum/(Cpsum - 1. / Wm(Nfz))
        Vlm(Npt) = Rr * Tt / (Wm(Nfz) * Pp)
        Wm(Npt) = Wm(Nfz)
        Dlvpt(Npt) = -1.
        Dlvtp(Npt) = 1.
        Totn(Npt) = Totn(Nfz)
        Ppp(Npt) = Pp
        Cpr(Npt) = Cpsum
        if (Tt >= Tg(1) * .8d0) then
           do i = Ngp1, Ngc
              if (En(i, Nfz) /= 0.) then
                 inc = i - Ng
                 if (Tt < (Temp(1, inc)-50.) .or. Tt > (Temp(2, inc)+50.)) go to 100
              end if
           end do
           go to 200
        end if
        go to 100
     else
        dlnt = (Ssum(Nfz) - Ssum(nnn)) / Cpsum
        Tln = Tln + dlnt
        if (abs(dlnt) < 0.5d-4) Convg = .true.
        Tt = exp(Tln)
     end if
  end do
  write(IOOUT, '(/" FROZEN DID NOT CONVERGE IN 8 ITERATIONS (FROZEN)")')
100 Tt = 0.
  Npt = Npt - 1
200 return
end subroutine FROZEN



subroutine GAUSS
!***********************************************************************
! SOLVE ANY LINEAR SET OF UP TO maxMat EQUATIONS
! NUMBER OF EQUATIONS = IMAT
!***********************************************************************
  use cea
  implicit none
! LOCAL VARIABLES
  integer, save:: i, imatp1, j, k, nn, nnp1
  real(8):: bigno = 1e25
  real(8), save:: coefx(50), tmp

! BEGIN ELIMINATION OF NNTH VARIABLE
  imatp1 = Imat + 1
  do nn = 1, Imat
     if (nn /= Imat) then
! SEARCH FOR MAXIMUM COEFFICIENT IN EACH ROW
        nnp1 = nn + 1
        do i = nn, Imat
           coefx(i) = bigno
           if (G(i, nn) /= 0.) then
              coefx(i) = 0.
              do j = nnp1, imatp1
                 coefx(i) = max(coefx(i), abs(G(i, j)))
              end do
              tmp = abs(G(i, nn))
              if (bigno * tmp > coefx(i)) then
                 coefx(i) = coefx(i) / tmp
              else
                 coefx(i) = bigno
              end if
           end if
        end do
! LOCATE ROW WITH SMALLEST MAXIMUM COEFFICIENT
        tmp = bigno
        i = 0
        do j = nn, Imat
           if (coefx(j) < tmp) then
              tmp = coefx(j)
              i = j
           end if
        end do
        if (i == 0) then
           Msing = nn
           go to 99999
! INDEX I LOCATES EQUATION TO BE USED FOR ELIMINATING THE NTH
! VARIABLE FROM THE REMAINING EQUATIONS
! INTERCHANGE EQUATIONS I AND NN
        else if (nn /= i) then
           do j = nn, imatp1
              tmp = G(i, j)
              G(i, j) = G(nn, j)
              G(nn, j) = tmp
           end do
        end if
     else if (G(nn, nn) == 0) then
        Msing = nn
        go to 99999
     end if
! DIVIDE NTH ROW BY NTH DIAGONAL ELEMENT AND ELIMINATE THE NTH
! VARIABLE FROM THE REMAINING EQUATIONS
     k = nn + 1
     tmp = G(nn, nn)
     if (tmp == 0.) then
        Msing = nn
        go to 99999
     else
        do j = k, imatp1
           G(nn, j) = G(nn, j) / tmp
        end do
        if (k /= imatp1) then
           do i = k, Imat
!DIR$ IVDEP
              do j = k, imatp1
                 G(i, j) = G(i, j) - G(i, nn) * G(nn, j)
              end do
           end do
        end if
     end if
  end do
! BACKSOLVE FOR THE VARIABLES
  k = Imat
100 j = k + 1
  X(k) = 0.0d0
  tmp = 0.0
  if (Imat >= j) then
     do i = j, Imat
        tmp = tmp + G(k, i) * X(i)
     end do
  end if
  X(k) = G(k, imatp1) - tmp
  k = k - 1
  if (k /= 0) go to 100
99999 return
end subroutine GAUSS



subroutine HCALC
!***********************************************************************
! CALCULATE PROPERTIES FOR TOTAL REACTANT USING THERMO DATA FOR
! ONE OR MORE REACTANTS. USED ONLY FOR SHOCK AND DETON PROBLEMS.
!***********************************************************************
  use cea
  implicit none
! LOCAL VARIABLES
  character(6), save:: date(maxNgc)
  character(2), save:: el(5)
  character(15), save:: sub
  integer, save:: i, icf, ifaz, itot, j, k, l, m, n, nall, nint, ntgas, ntot
  real(8), save:: bb(5), enj, er, sj, t1, t2, tem, thermo(9, 3), tsave

  tsave = Tt
  Tm = 0.
  if (Pp > 0.) Tm = log(Pp * Wmix)
  Ssum(Npt) = 0.
  Hpp(1) = 0.
  Hpp(2) = 0.
  Hsub0 = 0.
  Cpmix = 0.
  tem = (1. + Oxfl)
! LOOP ON REACTANTS.
! if oxidant, k = 1
! if fuel,    k = 2
  Nspr = Nspx
  do n = 1, Nreac
     k = 2
     if (Fox(n)(:1) == 'O' .or. Fox(n)(:1) == 'o') k = 1
     if (Tt == 0.) Tt = Rtemp(n)
     j = Jray(n)
     if (j == 0) then
! SEARCH FOR REACTANT IN STORED THERMO SPECIES. STORE INDEX IN JRAY(N).
        ifaz = 0
        do j = 1, Ngc
           if (Rname(n) == Prod(j) .or. '*' // Rname(n) == Prod(j)) then
              Jray(n) = j
              if (j > Ng) then
                 write(IOOUT, '(/" REACTANTS MUST BE GASEOUS FOR THIS PROBLEM (HCALC)")')
                 go to 20
              end if
              go to 50
           end if
        end do
! SEARCH THERMO.LIB FOR SPECIES.
        rewind IOTHM
        read(IOTHM) Tg, ntgas, ntot, nall
        Nspr = Nspr + 1
        do itot = 1, nall
           if (itot <= ntot) then
              icf = 3
              if (itot > ntgas) icf = 1
              read(IOTHM) sub, nint, date(Nspr), (el(j), bb(j), j = 1, 5), ifaz, &
                   t1, t2, Mw(Nspr), ((thermo(l, m), l = 1, 9), m = 1, icf)
           else
              read(IOTHM) sub, nint, date(Nspr), (el(j), bb(j), j = 1, 5), ifaz, &
                   t1, t2, Mw(Nspr), er
              if (nint /= 0) then
                 read(IOTHM) ((thermo(i, j), i = 1, 9), j = 1, nint)
                 icf = nint
              end if
           end if
           if (sub == Rname(n) .or. sub == '*' // Rname(n)) then
              if (ifaz <= 0 .and. nint > 0) then
                 do j = 1, 5
                    if (bb(j) == 0.) go to 2
                    Nfla(n) = j
                    Ratom(n, j) = el(j)
                    Rnum(n, j) = bb(j)
                 end do
2                Jray(n) = Nspr
                 j = Nspr
                 do l = 1, icf
                    do m = 1, 9
                       Coef(j, m, l) = thermo(m, l)
                    end do
                 end do
                 go to 50
              else
                 if (ifaz > 0) write(IOOUT, '(/" REACTANTS MUST BE GASEOUS FOR THIS PROBLEM (HCALC)")')
                 if (nint == 0) write(IOOUT, '(/" COEFFICIENTS FOR ", A15, " ARE NOT AVAILABLE (HCALC)")') Rname(n)
                 go to 20
              end if
           end if
        end do
        Nspr = Nspr - 1
        write(IOOUT, '(/" ERROR IN DATA FOR ", A15, " CHECK NAME AND TEMPERATURE", &
             & " RANGE IN", /, " thermo.inp (HCALC)")') Rname(n)
        Energy(n) = ' '
20      Tt = 0.
        Cpmix = 0.
        go to 100
     end if
! CALCULATE EN FOR REACTANT AND CALCULATE PROPERTIES.
50   if (Moles) enj = Pecwt(n) / Wp(k)
     if (.not. Moles) enj = Pecwt(n) / Rmw(n)
     enj = enj / tem
     if (k == 1) enj = enj * Oxfl
     Tln = log(Tt)
     En(j, Npt) = enj
     l = 1
     if (ifaz <= 0) then
        if (Tt > Tg(2)) l = 2
        if (Tt > Tg(3) .and. ifaz < 0) l = 3
     end if
     S(j) = ((((Coef(j, 7, l) / 4.) * Tt + Coef(j, 6, l) / 3.) * Tt + Coef(j, 5, l) / 2.) * Tt &
          + Coef(j, 4, l)) * Tt - (Coef(j, 1, l) * .5d0 / Tt + Coef(j, 2, l)) &
          / Tt + Coef(j, 3, l) * Tln + Coef(j, 9, l)
     H0(j) = ((((Coef(j, 7, l) / 5.) * Tt + Coef(j, 6, l) / 4.) * Tt + Coef(j, 5, l) / 3.) * Tt &
          + Coef(j, 4, l) / 2.) * Tt &
          - (Coef(j, 1, l) / Tt - Coef(j, 2, l) * Tln - Coef(j, 8, l)) / Tt + Coef(j, 3, l)
     Cp(j) = (((Coef(j, 7, l) * Tt + Coef(j, 6, l)) * Tt + Coef(j, 5, l)) * Tt &
          + Coef(j, 4, l)) * Tt + (Coef(j, 1, l) / Tt + Coef(j, 2, l)) / Tt + Coef(j, 3, l)
     if (H0(j) > -.01 .and. H0(j) < .01) H0(j) = 0.
! ADD CONTRIBUTION TO CP, H, AND S OF TOTAL REACTANT.
     Cpmix = Cpmix + Cp(j) * enj
! FOR CONDENSED SPECIES:  SJ = S(J)
     sj = S(j) - log(enj) - Tm
     Ssum(Npt) = Ssum(Npt) + enj * sj
     er = H0(j) * enj * Tt
     Hsub0 = Hsub0 + er
     Hpp(k) = Hpp(k) + er
  end do
  if (tsave /= 0.) Tt = tsave
100 return
end subroutine HCALC



subroutine INFREE(readOK, Cin, Ncin, Lcin, Dpin)
!***********************************************************************
! FREE-FORM READ FOR CEA.  READS AND DECIPHERS DATA FOR ONE DATASET.
!
! DEFINITIONS:
!   Ch1  - INDIVIDUAL CHARACTERS IN RECORD, MAXIMUM 132.
!   Nch1 - COLUMN NUMBER FOR THE LAST NON-BLANK CHARACTER IN RECORD.
!   Ncin - NUMBER OF VARIABLES IN DATASET.
!   Cin  - CHARACTER STRINGS IN DATASET. MAXIMUM 15 CHARACTERS.
!   LCIN - NEG. LENGTH OF LITERALS.  FOR NUMERICS, INDEX OF PREVIOUS
!          LITERAL.  ZERO FOR UNACCEPTIBLE VARIABLES.  VARIABLE
!          FOLLOWING "CASE" IS ALWAYS ASSUMED TO BE LITERAL.
!   Nb   - NUMBER OF DELIMITERS IN STRING.
!   Nx   - NUMBER OF CHARACTERS IN STRING.
!   Dpin - NUMERICS AS DOUBLE PRECISION VARIABLE.
!   Cnum - CHARACTER STRING REPRESENTING DATASET NUMBERS. MAXIMUM 24
!          CHARACTERS.
!***********************************************************************
  use cea
  implicit none
! DUMMY ARGUMENTS
  character(15):: Cin(maxNgc)
  integer:: Ncin
  integer:: Lcin(maxNgc)
  logical:: readOK
  real(8):: Dpin(maxNgc)
! LOCAL VARIABLES
  character(1), parameter:: nums(13) = ['+', '-', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '.']
  character(2), parameter:: numg(24) = &
       [character(2):: '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', &
                       '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24']
  character(1):: ch1(132), cx
  character(24):: cnum
  character(3):: fmtl(3)
  character(4):: w1
  integer:: i, ich1, j, kcin, nb, nch1, nx

  fmtl = [character(3):: '(g', '16', '.0)']

  Ncin = 1
  Lcin(1) = 0
  kcin = 0
  Dpin(1) = 0
100 nb = 1
  nx = 0
  cnum = ' '
  Cin(Ncin) = ' '
  ch1(1) = ' '
  nch1 = 1
! READ CHARACTERS, ONE AT A TIME
  read(IOINP, '(132a1)', END=500, ERR=500) ch1
! FIND FIRST AND LAST NON-BLANK CHARACTER
  do i = 132, 1, - 1
     nch1 = i
     if (ch1(i) /= ' ' .and. ch1(i) /= '	') go to 200
  end do
200 do i = 1, nch1
     ich1 = i
     if (ch1(i) /= ' ' .and. ch1(i) /= '	') go to 300
  end do
300 if (nch1 == 1 .or. ch1(ich1) == '#' .or. ch1(ich1) == '!') then
     write(IOOUT, '(1x, 80a1)') (ch1(i), i = 1, nch1)
     go to 100
  end if
  w1 = ch1(ich1) // ch1(ich1+1) // ch1(ich1+2) // ch1(ich1+3)
! IS STRING A KEYWORD SIGNALLING START OR END OF DATASET?
  if (w1 == 'ther' .or. w1 == 'tran' .or. w1 == 'prob' .or.  &
       w1 == 'reac' .or. w1 == 'outp' .or. w1 == 'omit' .or.  &
       w1 == 'only' .or. w1 == 'inse' .or. w1(1:3) == 'end') then
     if (Ncin == 1) then
        Cin(Ncin) = w1
        if (w1(1:3) == 'end' .or. w1 == 'ther' .or. w1 == 'tran') then
           write(IOOUT, '(1x, 80a1)') (ch1(i), i = 1, nch1)
           return
        end if
        ich1 = ich1 + 4
        nx = 4
        Lcin(1) = -4
     else
! KEYWORD READ FOR NEXT DATASET. END PROCESSING
        backspace IOINP
        if (nx == 0) Ncin = Ncin - 1
        return
     end if
  else if (Ncin == 1) then
     write(IOOUT, '(/" FATAL ERROR IN INPUT format (INFREE)")')
     go to 500
  end if
  write(IOOUT, '(1x, 80a1)') (ch1(i), i = 1, nch1)
  do 400 i = ich1, nch1
     cx = ch1(i)
! LOOK FOR DELIMITER STRINGS
     if (cx == ',' .and. (Lcin(Ncin) > 0 .or. nx == 0)) cx = ' '
     if (cx == '=' .and. (Lcin(Ncin) < 0 .or. nx == 0)) cx = ' '
     if (cx /= ' ' .and. cx /= '	') then
! LOOK FOR CHARACTER STRINGS
        nx = nx + 1
        if (Ncin > 1) then
           cnum(nx:nx) = cx
           if (nx <= 15) Cin(Ncin) = trim(cnum)
           if (nx == 1) then
! IS THIS A NUMERIC?
              do j = 1, 13
                 if (ch1(i) == nums(j)) then
                    Lcin(Ncin) = kcin
                    go to 310
                 end if
              end do
              Lcin(Ncin) = -1
              kcin = Ncin
           else if (Lcin(Ncin) < 0) then
              Lcin(Ncin) = -nx
           end if
310        nb = 1
        end if
        if (i < nch1 .or. Lcin(Ncin) < 0) go to 400
     end if
     if (nb == 1. .and. nx > 0) then
        if (Ncin > 0 .and. Lcin(Ncin) > 0) then
! CONVERT NUMERIC CHARACTER STRINGS TO real(8) VARIABLES (DPIN)
           fmtl(2) = numg(min(24, nx))
! INTERNAL READ TO CONVERT TO NUMERIC
           read(cnum, fmtl, ERR=320) Dpin(Ncin)
        end if
        go to 340
320     if (Cin(Ncin-1)(:4) /= 'case') write(IOOUT, '(/" WARNING!!  UNACCEPTABLE NUMBER ", a15, " (INFREE)")') Cin(i)
        Lcin(Ncin) = 0
340     Ncin = Ncin + 1
        Cin(Ncin) = ' '
        Lcin(Ncin) = 0
        Dpin(Ncin) = 0
        nx = 0
        cnum = ' '
     end if
     nb = nb + 1
400 continue
  if (nx > 0) then
     Ncin = Ncin + 1
     Lcin(Ncin) = 0
     Dpin(Ncin) = 0
  end if
  go to 100
500 readOK = .false.
  return
end subroutine INFREE



subroutine INPUT(readOK, caseOK, Ensert)
!***********************************************************************
! DECIPHER KEYWORDS, LITERAL VARIABLES, & NUMERICAL VARIABLES IN INPUT.
!***********************************************************************
  use cea
  implicit none
! DUMMY ARGUMENTS
  logical:: caseOK, readOK
  character(15):: Ensert(20)
! LOCAL VARIABLES
  character(26), parameter:: lc = 'abcdefghijklmnopqrstuvwxyz', uc = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  character(15), save:: cin(maxNgc), cx15
  character(4), save:: code, cx4
  character(1), save:: cx1
  character(2), save:: cx2
  character(3), save:: cx3
  logical, save:: eqrats, incd, phi, pltdat, reacts, refl
  integer, save:: i, ifrmla, ii, in, iv, ix, j, jj, k, lcin(maxNgc), ncin, nmix
  real(8), save:: denmtr, dpin(maxNgc), eratio, hr, mix(maxNgc), ur, xyz

  write(IOOUT, '(/, /)')
  caseOK = .true.
  Nonly = 0
  Nomit = 0
  Nsert = 0
  reacts = .false.
  Trace = 0
  Short = .false.
  Massf = .false.
  do i = 1, Ncol
     Debug(i) = .false.
  end do
  Nplt = 0
  SIunit = .true.
  pltdat = .false.
! CALL INFREE TO READ DATASET
100 call INFREE(readOK, cin, ncin, lcin, dpin)
  if (.not. readOK) go to 400
  code = trim(cin(1))
  if (code /= '    ') then
! STORE PRODUCT NAMES FROM 'ONLY' DATASET
     if (code == 'only') then
        Nonly = min(maxNgc, ncin-1)
        do i = 1, Nonly
           Prod(i) = cin(i+1)
        end do
! STORE CONDENSED PRODUCT NAMES FROM 'INSERT' DATASET
     else if (code == 'inse') then
        Nsert = min(20, ncin-1)
        do i = 1, Nsert
           Ensert(i) = cin(i+1)
        end do
! STORE PRODUCT NAMES FROM 'OMIT' DATASET
     else if (code == 'omit') then
! CHECK OMIT DATASET
        Nomit = min(maxNgc, ncin-1)
        do i = 1, Nomit
           Omit(i) = cin(i+1)
        end do
! KEYWORD 'THER' READ
! CALL UTHERM TO CONVERT formatTED THERMODYNAMIC DATA
     else if (code == 'ther') then
        Newr = .true.
        rewind IOTHM
        call UTHERM(readOK)
        if (.not. readOK) then
           write(IOOUT, '(/" FATAL ERROR IN DATASET (INPUT)")')
           go to 400
        end if
! KEYWORD 'TRAN' READ
! CALL UTRAN TO CONVERT formatTED TRANSPORT PROPERTIES
     else if (code == 'tran') then
        call UTRAN(readOK)
        if (.not. readOK) then
           write(IOOUT, '(/" FATAL ERROR IN DATASET (INPUT)")')
           go to 400
        end if
! PROCESS 'OUTP' DATASET.
     else if (code == 'outp') then
        do i = 2, ncin
           if (lcin(i) < 0) then
              cx2 = cin(i)(1:2)
              cx3 = cin(i)(1:3)
              cx4 = cin(i)(1:4)
              if (cx3 == 'cal') then
                 SIunit = .false.
              else if (cx4 == 'tran' .or. cx3 == 'trn') then
                 Trnspt = .true.
              else if (cx4 == 'trac') then
                 Trace = dpin(i+1)
              else if (cin(i)(:5) == 'short') then
                 Short = .true.
              else if (cin(i)(:5) == 'massf') then
                 Massf = .true.
              else if (cx3 == 'deb' .or. cx3 == 'dbg') then
                 do j = i + 1, ncin
                    if (lcin(j) /= i) cycle
                    k = int(dpin(j))
                    if (k <= Ncol) Debug(k) = .true.
                    lcin(j) = 0
                 end do
              else if (cx2 == 'si') then
                 SIunit = .true.
              else if (pltdat .and. Nplt < 20) then
                 Nplt = Nplt + 1
                 Pltvar(Nplt) = cin(i)
              else if (cx2 == 'pl') then
                 pltdat = .true.
              else
                 write(IOOUT, '("  WARNING!!  DID NOT RECOGNIZE ", a15, " (INPUT)"/)') cin(i)
              end if
           end if
        end do
! SORT AND STORE DATA FROM 'REAC' DATASET.
     else if (code == 'reac') then
        reacts = .true.
        Moles = .false.
        Nreac = 0
        do i = 1, maxR
           Pecwt(i) = -1.
        end do
        i = 1
140     i = i + 1
        if (i <= ncin) then
           if ( lcin(i) /= 0) then
              if (lcin(i) > 0) then
                 write(IOOUT, '(/" WARNING!!  LITERAL EXPECTED FOR ", a15, "(INPUT)")') cin(i)
                 go to 140
              end if
              cx15 = cin(i)
              cx1 = cx15(:1)
              cx2 = cx15(:2)
              cx3 = cx15(:3)
              cx4 = cx15(:4)
! NEW REACTANT
              if (cx2 /= 'na' .and. cx2 /= 'ox' .and. cx2 /= 'fu') then
! LOOK FOR PERCENTS
                 if (cx1 == 'm' .or. cx1 == 'w') then
                    if (lcin(i+1) > 0) then
                       i = i + 1
                       Pecwt(Nreac) = dpin(i)
                    else
                       caseOK = .false.
                       write(IOOUT, '(/" REACTANT AMOUNT MISSING (INPUT)")')
                    end if
                    if (cx1 == 'm' .and. Nreac == 1) Moles = .true.
                    if (cx1 == 'm' .and. .not. Moles .or. cx1 == 'w' .and. Moles) then
                       caseOK = .false.
                       write(IOOUT, '(/" MOLES AND WEIGHT PERCENTS SHOULD NOT BE MIXED (INPUT)")')
                    end if
                    go to 140
                 end if
! LOOK FOR TEMPERATURES
                 if (cx1 == 't') then
                    if (lcin(i+1) > 0) then
                       i = i + 1
                       Rtemp(Nreac) = dpin(i)
                       if (lcin(i-1) < 1) then
                          if (index(cx15, 'r') > 0) Rtemp(Nreac) = Rtemp(Nreac)/1.8d0
                          if (index(cx15, 'c') > 0) Rtemp(Nreac) = Rtemp(Nreac) + 273.15d0
                          if (index(cx15, 'f') > 0) Rtemp(Nreac) = (Rtemp(Nreac)-32)/1.8d0 + 273.15d0
                       end if
                    else
                       write(IOOUT, '(/" REACTANT TEMPERATURE MISSING (INPUT) ")')
                       caseOK = .false.
                    end if
                    go to 140
                 end if
! LOOK FOR ENTHALPY
                 if (cx1 == 'h' .or. cx1 == 'u') then
                    Energy(Nreac) = cx15
                    if (lcin(i+1) > 0) then
                       i = i + 1
                       Enth(Nreac) = dpin(i)*1000/Rr
                       if (index(cin(i-1), 'c') > 0) Enth(Nreac) = Enth(Nreac)*4.184d0
                       if (index(cin(i-1), 'k') > 0) Enth(Nreac) = Enth(Nreac)*1000
                    end if
                    go to 140
                 end if
! LOOK FOR DENSITY
                 if (cx3 == 'rho' .or. cx3 == 'den') then
                    if (lcin(i+1) > 0) then
                       i = i + 1
                       Dens(Nreac) = dpin(i)
                       if (index(cx15, 'kg') > 0) Dens(Nreac) = Dens(Nreac)/1000
                    end if
                    go to 140
                 end if
! CHECK FOR CHEMICAL SYMBOLS IN EXPLODED FORMULA
                 if ((lcin(i) == -1 .or. lcin(i) == -2) .and. index(uc, cx1) > 0) then
                    Energy(Nreac) = ' '
                    ifrmla = ifrmla + 1
                    Nfla(Nreac) = ifrmla
                    if (lcin(i) == -2) then
                       ix = index(lc, cx2(2:2))
                       if (ix > 0) cx2(2:2) = uc(ix:ix)
                    end if
                    Ratom(Nreac, ifrmla) = cx2
                    if (lcin(i+1) == i) then
                       Rnum(Nreac, ifrmla) = dpin(i+1)
                    else
                       Rnum(Nreac, ifrmla) = 1.
                    end if
                    i = i + 1
                    go to 140
                 end if
                 write(IOOUT, '(/" WARNING!! ", a15, " NOT RECOGNIZED (INPUT)")') cin(i)
              else
                 Nreac = min(Nreac+1, maxR)
                 Fox(Nreac) = trim(cx15)
                 i = i + 1
                 if (lcin(i) < 0) Rname(Nreac) = cin(i)
                 ifrmla = 0
                 Nfla(Nreac) = 0
                 Energy(Nreac) = 'lib'
                 Enth(Nreac) = 0
                 Jray(Nreac) = 0
                 Pecwt(Nreac) = -1
                 Rnum(Nreac, 1) = 0
                 Rmw(Nreac) = 0
                 Rtemp(Nreac) = 0
              end if
           end if
           go to 140
        end if
! SORT AND STORE INPUT FROM 'PROB' DATASET
     else if (code == 'prob') then
        Case = ' '
        do i = 1, maxPv
           P(i) = 0
           V(i) = 0
        end do
        do i = 1, maxT
           T(i) = 0
        end do
        P(1) = 1
        Trace = 0
        Lsave = 0
        R = Rr/4184.d0
        S0 = 0
        hr = 1.d30
        ur = 1.d30
        Tp = .false.
        Hp = .false.
        Sp = .false.
        Rkt = .false.
        Shock = .false.
        Detn = .false.
        Vol = .false.
        Ions = .false.
        Eql = .false.
        Froz = .false.
        Fac = .false.
        Debugf = .false.
        Acat = 0
        Ma = 0
        Nfz = 1
        Nsub = 0
        Nsup = 0
        Npp = 0
        Tcest = 3800
        do i = 1, Ncol
           Pcp(i) = 0
           Pcp(i+Ncol) = 0
           Supar(i) = 0
           Subar(i) = 0
           Mach1(i) = 0
           U1(i) = 0
        end do
        Gamma1 = 0
        phi = .false.
        eqrats = .false.
        incd = .false.
        refl = .false.
        Shkdbg = .false.
        Incdeq = .false.
        Incdfz = .false.
        Refleq = .false.
        Reflfz = .false.
        Np = 0
        Nt = 1
        Trnspt = .false.
! PROCESS LITERAL VARIABLES IN 'PROB' DATASET THAT DO NOT HAVE
! ASSOCIATED NUMERICAL DATA.
        outerLoop: do i = 2, ncin
           if (lcin(i) < 0) then
              do j = i + 1, ncin
                 if (lcin(j) == i) cycle outerLoop
              end do
              cx15 = cin(i)
              cx2 = cx15(:2)
              cx3 = cx15(:3)
              cx4 = cx15(:4)
              if (cx4 == 'case') then
                 Case = cin(i+1)
                 lcin(i+1) = 0
              else if (cx2 == 'tp' .or. cx2 == 'pt') then
                 Tp = .true.
              else if (cx2 == 'hp' .or. cx2 == 'ph') then
                 Hp = .true.
              else if (cx2 == 'sp' .or. cx2 == 'ps') then
                 Sp = .true.
              else if (cx2 == 'sv' .or. cx2 == 'vs') then
                 Sp = .true.
                 Vol = .true.
              else if (cx2 == 'uv' .or. cx2 == 'vu') then
                 Hp = .true.
                 Vol = .true.
              else if (cx2 == 'tv' .or. cx2 == 'vt') then
                 Tp = .true.
                 Vol = .true.
              else if (cx2 == 'ro' .or. cx3 == 'rkt') then
                 Rkt = .true.
              else if (cx3 == 'dbg' .or. cx3 == 'deb') then
                 Debugf = .true.
                 Shkdbg = .true.
                 Detdbg = .true.
              else if (cx3 == 'fac') then
                 Rkt = .true.
                 Eql = .true.
                 Fac = .true.
                 Froz = .false.
              else if (cx2 == 'eq') then
                 Eql = .true.
              else if (cx2 == 'fr' .or. cx2 == 'fz') then
                 Froz = .true.
              else if (cx2 == 'sh') then
                 Shock = .true.
              else if (cx3 == 'inc') then
                 Shock = .true.
                 incd = .true.
                 if (index(cx15, 'eq') > 0) Eql = .true.
                 if (index(cx15, 'fr') > 0) Froz = .true.
                 if (index(cx15, 'fz') > 0) Froz = .true.
              else if (cx3 == 'ref') then
                 Shock = .true.
                 refl = .true.
                 if (index(cx15, 'eq') > 0) Eql = .true.
                 if (index(cx15, 'fr') > 0) Froz = .true.
                 if (index(cx15, 'fz') > 0) Froz = .true.
              else if (cx3 == 'det') then
                 Detn = .true.
              else if (cx4 == 'ions') then
                 Ions = .true.
              else
                 write(IOOUT, '("  WARNING!!  DID NOT RECOGNIZE ", a15, " (INPUT)"/)') cx15
              end if
              lcin(i) = 0
           end if
        end do outerLoop
        iv = 2
        Nof = 0
        go to 200
     else if (code(1:3) == 'end') then
        if (Shock) then
           if (incd .and. Froz) Incdfz = .true.
           if (incd .and. Eql) Incdeq = .true.
           if (refl .and. Froz) Reflfz = .true.
           if (refl .and. Eql) Refleq = .true.
        end if
        Hsub0 = min(hr, ur)
        Size = 0.
        if (hr > .9d30) hr = 0.d0
        if (ur > .9d30) ur = 0.d0
        if (Trnspt) Viscns = .3125*sqrt(1.E5*Boltz/(Pi*Avgdr))
        if (SIunit) R = Rr/1000
        if (Detn .or. Shock) Newr = .true.
        if (.not. Short) then
           write(IOOUT, '(/" OPTIONS: TP=", L1, "  HP=", L1, "  SP=", L1, "  TV=", L1, &
                & "  UV=", L1, "  SV=", L1, "  DETN=", L1, "  SHOCK=", L1, &
                & "  REFL=", L1, "  INCD=", L1, /" RKT=", L1, "  FROZ=", L1, &
                & "  EQL=", L1, "  IONS=", L1, "  SIUNIT=", L1, "  DEBUGF=", L1, &
                & "  SHKDBG=", L1, "  DETDBG=", L1, "  TRNSPT=", L1)') Tp, (Hp .and. .not. Vol), Sp, (Tp .and. Vol), &
                (Hp .and. Vol), (Sp .and. Vol), Detn, Shock, refl, &
                incd, Rkt, Froz, Eql, Ions, SIunit, Debugf, Shkdbg, &
                Detdbg, Trnspt
           if (T(1) > 0) write(IOOUT, '(/" T,K =", 7F11.4)') (T(jj), jj = 1, Nt)
           write(IOOUT, '(/1p, " TRACE=", E9.2, "  S/R=", E13.6, "  H/R=", E13.6, "  U/R=",  E13.6)') Trace, S0, hr, ur
           if (Np > 0 .and. Vol) write(IOOUT, '(/" SPECIFIC VOLUME,M**3/KG =", 1p, (4E14.7))') (V(jj)*1.d-05, jj = 1, Np)
        end if
        if (Rkt) then
           if (Nt == 0) Hp = .true.
           if (.not. Short) then
              write(IOOUT, '(/" Pc,BAR =", 7F13.6)') (P(jj), jj = 1, Np)
              write(IOOUT, '(/" Pc/P =", 9F11.4)') (Pcp(jj), jj = 1, Npp)
              write(IOOUT, '(/" SUBSONIC AREA RATIOS =", (5F11.4))') (Subar(i), i = 1, Nsub)
              write(IOOUT, '(/" SUPERSONIC AREA RATIOS =", (5F11.4))') (Supar(i), i = 1, Nsup)
              write(IOOUT, '(/" NFZ=", i3, 1p, "  Mdot/Ac=", e13.6, "  Ac/At=", e13.6)') Nfz, Ma, Acat
           end if
        else
           if (.not. Vol .and. .not. Short) write(IOOUT, '(/" P,BAR =", 7F13.6)') (P(jj), jj = 1, Np)
        end if
        if (reacts) call REACT
        if (Nreac == 0 .or. Nlm <= 0) then
           write(IOOUT, '(/" ERROR IN REACTANTS DATASET (INPUT)")')
           caseOK = .false.
           write(IOOUT, '(/" FATAL ERROR IN DATASET (INPUT)")')
           go to 400
        end if
        if (Nof == 0) then
           Nof = 1
           Oxf(1) = 0
           if (Wp(2) > 0) then
              Oxf(1) = Wp(1) / Wp(2)
           else
              caseOK = .false.
              write(IOOUT, '(/" REACTANT AMOUNT MISSING (INPUT)")')
              write(IOOUT, '(/" FATAL ERROR IN DATASET (INPUT)")')
              go to 400
           end if
        else if (phi .or. eqrats) then
           do i = 1, Nof
              eratio = Oxf(i)
              if (eqrats) then
                 xyz = -eratio*Vmin(2) - Vpls(2)
                 denmtr = eratio*Vmin(1) + Vpls(1)
              else
                 xyz = -Vmin(2) - Vpls(2)
                 denmtr = eratio * (Vmin(1) + Vpls(1))
              end if
              if (abs(denmtr) < 1.d-30) then
                 caseOK = .false.
                 write(IOOUT, '(/" UNABLE TO PROCESS EQUIVALENCE RATIO =", E11.4, "(INPUT)")') eratio
                 write(IOOUT, '(/" FATAL ERROR IN DATASET (INPUT)")')
                 go to 400
              end if
              Oxf(i) = xyz / denmtr
           end do
        end if
        if (.not. Sp .and. .not. Tp .and. .not. Hp .and. .not. Rkt .and. .not. Detn .and. .not. Shock) then
           caseOK = .false.
           write(IOOUT, '(/" TYPE OF PROBLEM NOT SPECIFIED (INPUT)")')
        else if (Tp .and. T(1) <= 0) then
           caseOK = .false.
           write(IOOUT, '(/" ASSIGNED VALUES OF TEMPERATURE ARE MISSING IN prob", " DATASET (INPUT)")')
        else if (Np <= 0) then
           caseOK = .false.
           write(IOOUT, '(/" ASSIGNED PRESSURE (OR DENSITY) MISSING IN prob", " DATASET (INPUT)")')
        end if
        if (.not. (caseOK .and. Nlm > 0)) write(IOOUT, '(/" FATAL ERROR IN DATASET (INPUT)")')
        go to 400
     else
        write(IOOUT, '(/" WARNING!!  A KEYWORD IS MISSING (INPUT)")')
     end if
  end if
  go to 100
! PROCESS NUMERICAL DATA FOLLOWING 'PROB' LITERALS
200 in = 0
  nmix = 0
  ii = iv
  do i = ii, ncin
     iv = i
     if (lcin(i) /= 0) then
        if (lcin(i) < 0) then
           if (in > 0) exit
           in = i
        else
           if (lcin(i) /= in) exit
           nmix = nmix + 1
           mix(nmix) = dpin(i)
           lcin(i) = 0
        end if
     end if
  end do
  if (nmix <= 0) then
     if (iv < ncin) go to 200
     go to 100
  end if
  cx15 = cin(in)
  cx1 = cx15(:1)
  cx2 = cx15(:2)
  cx3 = cx15(:3)
  cx4 = cx15(:4)
  if (cx1 == 't') then
     Nt = nmix
     if (nmix > maxMix) then
        Nt = maxMix
        write(IOOUT, '(/" NOTE!! MAXIMUM NUMBER OF ASSIGNED ", a5, " VALUES IS", I3, " (INPUT)", /)') 't', Nt
     end if
     do i = 1, Nt
        if (cx4 /= 'tces') then
           T(i) = mix(i)
           if (lcin(in) < -1) then
              if (index(cx15, 'r') > 0) T(i) = T(i)/1.8d0
              if (index(cx15, 'c') > 0) T(i) = T(i) + 273.15d0
              if (index(cx15, 'f') > 0) T(i) = (T(i)-32.d0) / 1.8d0 + 273.15d0
           end if
        end if
     end do
  else if ((cx2 == 'pc' .or. cx2 == 'pi') .and. index(cx15(3:15), 'p') > 0 .and. index(cx15, 'psi') == 0) then
     Npp = nmix
     if (nmix > 2*Ncol) then
        Npp = 2*Ncol
        write(IOOUT, '(/" NOTE!! MAXIMUM NUMBER OF ASSIGNED ", a5, " VALUES IS", I3, " (INPUT)", /)') 'pcp', Npp
     end if
     do i = 1, Npp
        Pcp(i) = mix(i)
     end do
  else if (cx1 == 'p' .and. cx3 /= 'phi') then
     Np = nmix
     if (nmix > maxPv) then
        Np = maxPv
        write(IOOUT, '(/" NOTE!! MAXIMUM NUMBER OF ASSIGNED ", a5, " VALUES IS", I3, " (INPUT)", /)') 'p', Np
     end if
     do i = 1, Np
        P(i) = mix(i)
        if (index(cx15, 'psi') /= 0) then
           P(i) = P(i)/14.696006d0
        else if (index(cx15, 'mmh') /= 0) then
           P(i) = P(i)/760.d0
        else if (index(cx15, 'atm') == 0) then
           cycle
        end if
        P(i) = P(i)*1.01325d0
     end do
  else if (cx3 == 'rho') then
     xyz = 1.d02
     if (index(cx15, 'kg') /= 0) xyz = 1.d05
     Np = nmix
     if (nmix > maxPv) then
        Np = maxPv
        write(IOOUT, '(/" NOTE!! MAXIMUM NUMBER OF ASSIGNED ", a5, " VALUES IS", I3, " (INPUT)", /)') 'rho', Np
     end if
     do i = 1, Np
        V(i) = xyz/mix(i)
     end do
  else if (cx1 == 'v') then
     xyz = 1.d02
     if (index(cx15, 'kg') /= 0) xyz = 1.d05
     Np = nmix
     if (nmix > maxPv) then
        Np = maxPv
        write(IOOUT, '(/" NOTE!! MAXIMUM NUMBER OF ASSIGNED ", a5, " VALUES IS", I3, " (INPUT)", /)') 'v', Np
     end if
     do i = 1, Np
        V(i) = mix(i)*xyz
     end do
  else if (cx3 == 'nfz' .or. cx3 == 'nfr') then
     Nfz = int(mix(1))
     Froz = .true.
  else if (cx4 == 'tces') then
     Tcest = mix(1)
  else if (cx4 == 'trac') then
     Trace = mix(1)
  else if (cx3 == 's/r') then
     S0 = mix(1)
  else if (cx3 == 'u/r' .or. cx2 == 'ur') then
     ur = mix(1)
  else if (cx3 == 'h/r' .or. cx2 == 'hr') then
     hr = mix(1)
  else if (cx2 == 'u1') then
     Nsk = nmix
     if (nmix > Ncol) then
        Nsk = Ncol
        write(IOOUT, '(/" NOTE!! MAXIMUM NUMBER OF ASSIGNED ", a5, " VALUES IS", I3, " (INPUT)", /)') 'u1', Nsk
     end if
     do i = 1, Nsk
        U1(i) = mix(i)
     end do
  else if (cx4 == 'mach') then
     Nsk = nmix
     if (nmix > Ncol) then
        Nsk = Ncol
        write(IOOUT, '(/" NOTE!! MAXIMUM NUMBER OF ASSIGNED ", a5, " VALUES IS", I3, " (INPUT)", /)') 'mach1', Nsk
     end if
     do i = 1, Nsk
        Mach1(i) = mix(i)
     end do
  else if (cx3 == 'sub') then
     Nsub = nmix
     if (nmix > 13) then
        Nsub = 13
        write(IOOUT, '(/" NOTE!! MAXIMUM NUMBER OF ASSIGNED ", a5, " VALUES IS", I3, " (INPUT)", /)') 'subar', Nsub
     end if
     do i = 1, Nsub
        Subar(i) = mix(i)
     end do
  else if (cx3 == 'sup') then
     Nsup = nmix
     if (nmix > 13) then
        Nsup = 13
        write(IOOUT, '(/" NOTE!! MAXIMUM NUMBER OF ASSIGNED ", a5, " VALUES IS", I3, " (INPUT)", /)') 'supar', Nsup
     end if
     do i = 1, Nsup
        Supar(i) = mix(i)
     end do
  else if (cx2 == 'ac') then
     Acat = mix(1)
  else if (cx4 == 'mdot' .or. cx2 == 'ma') then
     Ma = mix(1)
  else if (cx4 == 'case') then
     Case = cin(in+1)
     lcin(in+1) = 0
  else if (Nof == 0 .and. (cx3 == 'phi' .or. cx3 == 'o/f' .or. cx3 == 'f/a' .or. cx2 == '%f' .or. cx1 == 'r')) then
     Nof = nmix
     if (nmix > maxMix) then
        Nof = maxMix
        write(IOOUT, '(/" NOTE!! MAXIMUM NUMBER OF ASSIGNED ", a5, " VALUES IS", I3, " (INPUT)", /)') 'o/f', Nof
     end if
     do k = 1, Nof
        Oxf(k) = mix(k)
     end do
     if (cx3 == 'phi') then
        phi = .true.
     else if (cx1 == 'r') then
        eqrats = .true.
     else if (cx3 == 'f/a') then
        do k = 1, Nof
           if (Oxf(k) > 0) Oxf(k) = 1/Oxf(k)
        end do
     else if (cx4 == '%fue') then
        do k = 1, Nof
           if (Oxf(k) > 0) Oxf(k) = (100 - Oxf(k)) / Oxf(k)
        end do
     end if
  else
     write(IOOUT, '("  WARNING!!  DID NOT RECOGNIZE ", a15, " (INPUT)"/)') cx15
  end if
  if (iv >= ncin) go to 100
  go to 200
400 return
end subroutine INPUT



subroutine MATRIX
!***********************************************************************
! SET UP ITERATION OR DERIVATIVE MATRIX.
!***********************************************************************
  use cea
  implicit none
! LOCAL VARIABLES
  integer, save:: i, iq, iq2, iq3, isym, j, k, kk, kmat
  real(8), save:: energyl, f, h, ss, sss, term, term1

  iq = Nlm + Npr
  Iq1 = iq + 1
  iq2 = Iq1 + 1
  iq3 = iq2 + 1
  kmat = iq3
  if (.not. Convg .and. Tp) kmat = iq2
  Imat = kmat - 1
! CLEAR MATRIX STORAGES TO ZERO
  do i = 1, Imat
     do k = 1, kmat
        G(i, k) = 0
     end do
  end do
  G(iq2, Iq1) = 0
  sss = 0
  Hsum(Npt) = 0
! BEGIN SET-UP OF ITERATION OR DERIVATIVE MATRIX
  do j = 1, Ng
     Mu(j) = H0(j) - S(j) + Enln(j) + Tm
     if (En(j, Npt) /= 0) then
        h = H0(j) * En(j, Npt)
        f = Mu(j) * En(j, Npt)
        ss = h - f
        term1 = h
        if (kmat == iq2) term1 = f
        do i = 1, Nlm
           if (A(i, j) /= 0) then
              term = A(i, j) * En(j, Npt)
              do k = i, Nlm
                 G(i, k) = G(i, k) + A(k, j) * term
              end do
              G(i, Iq1) = G(i, Iq1) + term
              G(i, iq2) = G(i, iq2) + A(i, j) * term1
              if (.not. (Convg .or. Tp)) then
                 G(i, iq3) = G(i, iq3) + A(i, j) * f
                 if (Sp) G(iq2, i) = G(iq2, i) + A(i, j) * ss
              end if
           end if
        end do
        if (kmat /= iq2) then
           if (Convg .or. Hp) then
              G(iq2, iq2) = G(iq2, iq2) + H0(j) * h
              if (.not. Convg) then
                 G(iq2, iq3) = G(iq2, iq3) + H0(j) * f
                 G(Iq1, iq3) = G(Iq1, iq3) + f
              end if
           else
              G(iq2, Iq1) = G(iq2, Iq1) + ss
              G(iq2, iq2) = G(iq2, iq2) + H0(j) * ss
              G(iq2, iq3) = G(iq2, iq3) + Mu(j) * ss
              G(Iq1, iq3) = G(Iq1, iq3) + f
           end if
        end if
        G(Iq1, iq2) = G(Iq1, iq2) + term1
     end if
  end do
! CONDENSED SPECIES
  if (Npr /= 0) then
     do k = 1, Npr
        j = Jcond(k)
        kk = Nlm + k
        Mu(j) = H0(j) - S(j)
        do i = 1, Nlm
           G(i, kk) = A(i, j)
           G(i, kmat) = G(i, kmat) - A(i, j) * En(j, Npt)
        end do
        G(kk, iq2) = H0(j)
        G(kk, kmat) = Mu(j)
        Hsum(Npt) = Hsum(Npt) + H0(j) * En(j, Npt)
        if (Sp) then
           sss = sss + S(j) * En(j, Npt)
           G(iq2, kk) = S(j)
        end if
     end do
  end if
  sss = sss + G(iq2, Iq1)
  Hsum(Npt) = Hsum(Npt) + G(Iq1, iq2)
  G(Iq1, Iq1) = Sumn - Enn
! REFLECT SYMMETRIC PORTIONS OF THE MATRIX
  isym = Iq1
  if (Hp .or. Convg) isym = iq2
  do i = 1, isym
!DIR$ IVDEP
     do j = i, isym
        G(j, i) = G(i, j)
     end do
  end do
! COMPLETE THE RIGHT HAND SIDE
  if (.not. Convg) then
     do i = 1, Nlm
        G(i, kmat) = G(i, kmat) + B0(i) - G(i, Iq1)
     end do
     G(Iq1, kmat) = G(Iq1, kmat) + Enn - Sumn
! COMPLETE ENERGY ROW AND TEMPERATURE COLUMN
     if (kmat /= iq2) then
        if (Sp) energyl = S0 + Enn - Sumn - sss
        if (Hp) energyl = Hsub0/Tt - Hsum(Npt)
        G(iq2, iq3) = G(iq2, iq3) + energyl
        G(iq2, iq2) = G(iq2, iq2) + Cpsum
     end if
  else
     if (Pderiv) then
! PDERIV = .true.-- SET UP MATRIX TO SOLVE FOR DLVPT
        G(Iq1, iq2) = Enn
        do i = 1, iq
           G(i, iq2) = G(i, Iq1)
        end do
     end if
     G(iq2, iq2) = G(iq2, iq2) + Cpsum
  end if
  if (Vol .and. .not. Convg) then
! CONSTANT VOLUME MATRIX
     if (kmat == iq2) then
        do i = 1, iq
           G(i, Iq1) = G(i, iq2)
        end do
     else
!DIR$ IVDEP
        do i = 1, iq
           G(Iq1, i) = G(iq2, i) - G(Iq1, i)
           G(i, Iq1) = G(i, iq2) - G(i, Iq1)
           G(i, iq2) = G(i, iq3)
        end do
        G(Iq1, Iq1) = G(iq2, iq2) - G(Iq1, iq2) - G(iq2, Iq1)
        G(Iq1, iq2) = G(iq2, iq3) - G(Iq1, iq3)
        if (Hp) G(Iq1, iq2) = G(Iq1, iq2) + Enn
     end if
     kmat = Imat
     Imat = Imat - 1
  end if
end subroutine MATRIX



subroutine NEWOF
!***********************************************************************
! CALCULATE NEW VALUES OF B0 AND HSUB0 FOR NEW OF RATIO
!***********************************************************************
  use cea
  implicit none
! LOCAL VARIABLES
  integer, save:: i, j
  real(8), save:: assval, bigb, bratio, dbi, smalb, tem, v1, v2

  if (.not. Short) write(IOOUT, '(/" O/F = ", f10.6)') Oxfl
  Eqrat = 0
  tem = Oxfl + 1
  v2 = (Oxfl * Vmin(1) + Vmin(2)) / tem
  v1 = (Oxfl * Vpls(1) + Vpls(2)) / tem
  if (v2 /= 0.) Eqrat = abs(v1/v2)
  do i = 1, Nlm
     B0(i) = (Oxfl * B0p(i, 1) + B0p(i, 2)) / tem
     dbi = abs(B0(i))
     if (i == 1) then
        bigb = dbi
        smalb = dbi
     else if (dbi /= 0.) then
        if (dbi < smalb) smalb = dbi
        if (dbi > bigb) bigb = dbi
     end if
  end do
  Bcheck = bigb * .000001d0
! CALCUALTE MOLECULAR WEIGHT OF TOTAL REACTANT, WMIX.
  if (Am(1) /= 0 .and. Am(2) /= 0) then
     Wmix = (Oxfl + 1) * Am(1) * Am(2) / (Am(1) + Oxfl * Am(2))
  else
     Wmix = Am(2)
     if (Am(2) == 0.0) Wmix = Am(1)
  end if
  Npt = 1
! IF ASSIGNED U OR H NOT GIVEN IN PROB DATA, INITIAL HSUB0 = 1.d30
  if (Size == 0) assval = Hsub0
  if (assval >= 1.d30) Hsub0 = (Oxfl * Hpp(1) + Hpp(2)) / tem
! NOTE THAT "BRATIO" IS "BRATIO" IN SEC 3.2 IN RP-1311.
  bratio = smalb / bigb
  Size = 18.420681d0
  if (bratio < 1.d-5) Size = log(1000/bratio)
  Jsol = 0
  Jliq = 0
  if (.not. Short) then
     write(IOOUT, '(/, 23x, "EFFECTIVE FUEL", 5x, "EFFECTIVE OXIDANT", 8x, "MIXTURE")')
     if (Vol) write(IOOUT, '(" INTERNAL ENERGY", 11x, "u(2)/R", 14x, "u(1)/R", 14x, "u0/R")')
     if (.not. Vol) write(IOOUT, '(" ENTHALPY", 18x, "h(2)/R", 14x, "h(1)/R", 15x, "h0/R")')
     write(IOOUT, '(" (KG-MOL)(K)/KG", 4x, e18.8, 2e20.8)') Hpp(2), Hpp(1), Hsub0
     write(IOOUT, '(/" KG-FORM.WT./KG", 13x, "bi(2)", 15x, "bi(1)", 15x, "b0i")')
  end if
  do i = 1, Nlm
     j = Jcm(i)
     if (.not. Short) write(IOOUT, '(1x, a16, 3e20.8)') Prod(j), B0p(i, 2), B0p(i, 1), B0(i)
  end do
  return
end subroutine NEWOF



subroutine OUT1
!***********************************************************************
! OUT1 WRITES REACTANT AND FUEL-OXIDANT RATIO INformatION.
! ENTRY OUT2 WRITES THERMODYNAMIC PROPERTIES.
! ENTRY OUT3 WRITES MOLE FRACTIONS.
! ENTRY OUT4 WRITES TRANSPORT PROPERTIES.
!
! NOTE - ROCKET, SHOCK, AND DETON PROBLEMS HAVE ADDITIONAL OUTPUT.
!***********************************************************************
  use cea
  implicit none
! LOCAL VARIABLES
  character(15), save:: fc, fgi, fh, fp, frh, fs, fu
  character(4), save:: mamo
  integer, save:: i, im, ione, j, k, kin, m, mcond, mcondf, mcp, mdvp, mdvt, meq, mfa, &
       mg, mgam, mh, mie, mm, mmw, mof, mp, mpf, mph, mpn, mpnf, mrho, ms, &
       mson, mt, mvis, n, notuse
  logical, save:: kOK
  real(8), save:: pfactor, pfuel, phi, rho, tem, tra, vnum


  write(IOOUT, '(" CASE = ", a15)') Case
  if (Moles) then
     write(IOOUT, '(/13X, "REACTANT", 20x, a11, "      ENERGY", 6x, "TEMP")') '   MOLES   '
     if (.not. SIunit) write(IOOUT, '(57X, " CAL/MOL ", 6x, "K")')
     if (SIunit) write(IOOUT, '(57X, "KJ/KG-MOL", 6x, "K")')
  else
     write(IOOUT, '(/13X, "REACTANT", 20x, a11, "      ENERGY", 6x, "TEMP")') 'WT FRACTION'
     if (.not. SIunit) write(IOOUT, '(42X, "(SEE NOTE)      CAL/MOL       K  ")')
     if (SIunit) write(IOOUT, '(42X, "(SEE NOTE)     KJ/KG-MOL      K  ")')
  end if
  do n = 1, Nreac
     write(IOOUT, '(1x, a8, 4x, a15, 11x, f12.7, f14.3, f11.3)') Fox(n), Rname(n), Pecwt(n), Enth(n)*R, Rtemp(n)
  end do
  phi = 0.
  tem = (Vpls(1)+Vmin(1))*Oxfl
  if (ABS(tem) >= 1.d-3) phi = -(Vmin(2)+Vpls(2))/tem
  if (Fox(1) == 'NAME') then
     pfuel = 0.
  else
     pfuel = 100.d0/(1.d0+Oxfl)
  end if
  if (Rh(1) /= 0. .or. Rh(2) /= 0.) then
     if (Rh(1) == 0. .or. Rh(2) == 0.) then
        rho = max(Rh(1), Rh(2))
     else
        rho = (Oxfl+1.)*Rh(1)*Rh(2)/(Rh(1)+Oxfl*Rh(2))
     end if
     if (SIunit) then
        rho = rho*1000.d0
        write(IOOUT, '(/" REACTANT DENSITY=", F8.2, " KG/CU M")') rho
     else
        write(IOOUT, '(/" REACTANT DENSITY=", F8.4, " G/CC")') rho
     end if
  end if
  write(IOOUT, '(/" O/F=", F11.5, 2X, "%FUEL=", F10.6, 2X, "R,EQ.RATIO=", F9.6, 2X, &
       & "PHI,EQ.RATIO=", F9.6)') Oxfl, pfuel, Eqrat, phi
  return
!***********************************************************************
  entry OUT2
  ione = 0
  if (Rkt .and. .not. Page1) then
     ione = 2
     if (Iopt /= 0) ione = 3
  end if
! SET MXX ARRAY FOR PLOTTING PARAMETERS
  mp     = 0
  mt     = 0
  mrho   = 0
  mh     = 0
  mie    = 0
  mg     = 0
  ms     = 0
  mm     = 0
  mcp    = 0
  mgam   = 0
  mson   = 0
  mcond  = 0
  mvis   = 0
  mpn    = 0
  mpf    = 0
  mof    = 0
  mph    = 0
  meq    = 0
  mfa    = 0
  mmw    = 0
  mdvt   = 0
  mdvp   = 0
  mcondf = 0
  mpnf   = 0
  do 100 i = 1, Nplt
     if (index(Pltvar(i)(2:), '1') == 0) then
        if (index(Pltvar(i)(1:), 'dlnt') /= 0) then
           mdvt = i
        else if (index(Pltvar(i)(1:), 'dlnp') /= 0) then
           mdvp = i
        else if (Pltvar(i)(:4) == 'pran') then
           if (index(Pltvar(i)(3:), 'fz') /= 0  .or.  &
                index(Pltvar(i)(3:), 'fr') /= 0) then
              mpnf = i
           else
              mpn = i
           end if
        else if (Pltvar(i)(:4) == 'cond') then
           if (index(Pltvar(i)(3:), 'fz') /= 0  .or.  &
                index(Pltvar(i)(3:), 'fr') /= 0) then
              mcondf = i
           else
              mcond = i
           end if
        else if (Pltvar(i)(:3) == 'phi') then
           mph = i
        else if (Pltvar(i)(:2) == 'p ') then
           mp = i
        else if (Pltvar(i)(:1) == 't') then
           mt = i
        else if (Pltvar(i)(:3) == 'rho') then
           mrho = i
        else if (Pltvar(i)(:1) == 'h') then
           mh = i
        else if (Pltvar(i)(:1) == 'u') then
           mie = i
        else if (Pltvar(i)(:3) == 'gam') then
           mgam = i
        else if (Pltvar(i)(:3) == 'son') then
           mson = i
        else if (Pltvar(i)(:2) == 'g ') then
           mg = i
        else if (Pltvar(i)(:2) == 's ') then
           ms = i
        else if (Pltvar(i)(:1) == 'm' .and. Pltvar(i)(:2) /= 'ma') then
           if (.not. Gonly .and. Pltvar(i)(:2) == 'mw') then
              mmw = i
           else
              mm = i
           end if
        else if (Pltvar(i)(:2) == 'cp') then
           mcp = i
        else if (Pltvar(i)(:3) == 'vis') then
           mvis = i
        else if (Pltvar(i)(:3) == 'o/f') then
           mof = i
        else if (Pltvar(i)(:2) == '%f') then
           mpf = i
        else if (Pltvar(i)(:3) == 'f/a') then
           mfa = i
        else if (Pltvar(i)(:1) == 'r') then
           meq = i
        end if
     end if
100 continue
  do i = Iplt + 1, Iplt + Npt
     if (mof > 0) Pltout(i, mof) = Oxfl
     if (mpf > 0) Pltout(i, mpf) = pfuel
     if (mph > 0) Pltout(i, mph) = phi
     if (mfa > 0) Pltout(i, mfa) = 1.d0/Oxfl
     if (meq > 0) Pltout(i, meq) = Eqrat
  end do
  if (SIunit) then
     pfactor = 1.d0
     fp = 'P, BAR'
     vnum = 1.d05
     frh = 'RHO, KG/CU M'
     fh = 'H, KJ/KG'
     fu = 'U, KJ/KG'
     fgi = 'G, KJ/KG'
     fs = 'S, KJ/(KG)(K)'
     fc = 'Cp, KJ/(KG)(K)'
  else
     pfactor = 1.d0/1.01325d0
     fp = 'P, ATM'
     vnum = 100.d0
     frh = 'RHO, G/CC'
     fh = 'H, CAL/G'
     fu = 'U, CAL/G'
     fgi = 'G, CAL/G'
     fs = 'S, CAL/(G)(K)'
     fc = 'Cp, CAL/(G)(K)'
  end if
  Fmt(4) = Fmt(6)
! PRESSURE
  call VARFMT(Ppp)
  do i = 1, Npt
     X(i) = Ppp(i)*pfactor
     if (Nplt /= 0 .and. i > ione) then
        if (mp > 0) Pltout(i+Iplt-ione, mp) = X(i)
        if (mt > 0) Pltout(i+Iplt-ione, mt) = Ttt(i)
     end if
  end do
  write(IOOUT, Fmt) fp, (X(j), j=1, Npt)
! TEMPERATURE
  Fmt(4) = '13'
  Fmt(5) = ' '
  Fmt(7) = '2,'
  write(IOOUT, Fmt) 'T, K            ', (Ttt(j), j=1, Npt)
! DENSITY
  do i = 1, Npt
     if (Vlm(i) /= 0.) X(i) = vnum/Vlm(i)
     if (Nplt /= 0 .and. i > ione .and. mrho > 0) &
          Pltout(i+Iplt-ione, mrho) = X(i)
  end do
  call EFMT(Fmt(4), frh, X)
! ENTHALPY
  do i = 1, Npt
     X(i) = Hsum(i)*R
     if (Nplt /= 0 .and. i > ione .and. mh > 0) &
          Pltout(i+Iplt-ione, mh) = X(i)
  end do
  Fmt(4) = Fmt(6)
  call VARFMT(X)
  write(IOOUT, Fmt) fh, (X(j), j=1, Npt)
! INTERNAL ENERGY
  do i = 1, Npt
     X(i) = (Hsum(i)-Ppp(i)*Vlm(i)/Rr)*R
     if (Nplt /= 0 .and. i > ione .and. mie > 0) &
          Pltout(i+Iplt-ione, mie) = X(i)
  end do
  call VARFMT(X)
  write(IOOUT, Fmt) fu, (X(j), j=1, Npt)
! GIBBS ENERGY
  do i = 1, Npt
     X(i) = (Hsum(i)-Ttt(i)*Ssum(i))*R
     if (Nplt /= 0 .and. i > ione) then
        if (mg > 0) Pltout(i+Iplt-ione, mg) = X(i)
        if (mm > 0) Pltout(i+Iplt-ione, mm) = Wm(i)
        if (mmw > 0) Pltout(i+Iplt-ione, mmw) = 1.d0/Totn(i)
        if (ms > 0) Pltout(i+Iplt-ione, ms) = Ssum(i)*R
        if (mcp > 0) Pltout(i+Iplt-ione, mcp) = Cpr(i)*R
        if (mgam > 0) Pltout(i+Iplt-ione, mgam) = Gammas(i)
        if (mdvt > 0) Pltout(i+Iplt-ione, mdvt) = Dlvtp(i)
        if (mdvp > 0) Pltout(i+Iplt-ione, mdvp) = Dlvpt(i)
     end if
  end do
  call VARFMT(X)
  write(IOOUT, Fmt) fgi, (X(j), j=1, Npt)
! ENTROPY
  Fmt(4) = '13'
  Fmt(5) = ' '
  Fmt(7) = '4,'
  write(IOOUT, Fmt) fs, (Ssum(j)*R, j=1, Npt)
  write(IOOUT, *)
! MOLECULAR WEIGHT
  Fmt(7) = '3,'
  write(IOOUT, Fmt) 'M, (1/n)        ', (Wm(j), j=1, Npt)
  if (.not. Gonly) write(IOOUT, Fmt) 'MW, MOL WT      ', &
       (1.d0/Totn(j), j=1, Npt)
! (DLV/DLP)T
  Fmt(7) = '5,'
  if (Eql) write(IOOUT, Fmt) '(dLV/dLP)t      ', (Dlvpt(j), j=1, Npt)
! (DLV/DLT)P
  Fmt(7) = '4,'
  if (Eql) write(IOOUT, Fmt) '(dLV/dLT)p      ', (Dlvtp(j), j=1, Npt)
! HEAT CAPACITY
  write(IOOUT, Fmt) fc, (Cpr(j)*R, j=1, Npt)
! GAMMA(S)
  Fmt(7) = '4,'
  write(IOOUT, Fmt) 'GAMMAs          ', (Gammas(j), j=1, Npt)
! SONIC VELOCITY
  Fmt(7) = '1,'
  do i = 1, Npt
     Sonvel(i) = (Rr*Gammas(i)*Ttt(i)/Wm(i))**.5
     if (Nplt /= 0 .and. i > ione .and. mson > 0) &
          Pltout(i+Iplt-ione, mson) = Sonvel(i)
  end do
  write(IOOUT, Fmt) 'SON VEL,M/SEC   ', (Sonvel(j), j=1, Npt)
  return
!***********************************************************************
  entry OUT3
  tra = 5.d-6
  if (Trace /= 0.) tra = Trace
! MASS OR MOLE FRACTIONS 
  if (Massf) then
     mamo = 'MASS'
  else
     mamo = 'MOLE'
  end if
  if (Eql) then
     write(IOOUT, '(/1x, A4, " FRACTIONS"/)') mamo
     notuse = 0
     do k = 1, Ngc
        kOK = .true.
        if (k > Ng .and. k < Ngc .and. Prod(k) == Prod(k+1)) then
           kOK = .false.
           im = 0
           go to 120
        end if
        do m = 1, Nplt
           im = 0
           if (Pltvar(m) == Prod(k) .or. '*'//Pltvar(m) == Prod(k)) &
                then
              im = m
              go to 120
           end if
        end do
120     kin = 0
        do i = 1, Npt
           if (Massf) then
              tem = Mw(k)
           else
              tem = 1.d0/Totn(i)
           end if
           if (k <= Ng) then
              X(i) = En(k, i)*tem
           else
              if (Prod(k) /= Prod(k-1)) X(i) = 0.d0
              if (En(k, i) > 0.d0) X(i) = En(k, i)*tem
           end if
           if (Nplt /= 0 .and. i > ione .and. im > 0) &
                Pltout(i+Iplt-ione, im) = X(i)
           if (kOK .and. X(i) >= tra) kin = 1
        end do
        if (kin == 1) then
           if (Trace == 0.) then
              write(IOOUT, '(1x, A15, F9.5, 12F9.5)') Prod(k), (X(i), i=1, Npt)
           else
              call EFMT(Fmt(4), Prod(k), X)
           end if
           if (Prod(k) == Omit(notuse)) notuse = notuse - 1
        else if (Prod(k) /= Prod(k-1)) then
           notuse = notuse + 1
           Omit(notuse) = Prod(k)
        end if
     end do
  end if
  write(IOOUT, '(/"  * THERMODYNAMIC PROPERTIES FITTED TO", F7.0, "K")') Tg(4)
  if (.not. Short) then
     write(IOOUT, '(/"    PRODUCTS WHICH WERE CONSIDERED BUT WHOSE ", A4, &
          & " FRACTIONS", /"    WERE LESS THAN", 1PE13.6, &
          & " FOR ALL ASSIGNED CONDITIONS"/)') mamo, tra
     write(IOOUT, '(5(1x, A15))') (Omit(i), i=1, notuse)
  end if
  if (.not. Moles) write(IOOUT, '(/" NOTE. WEIGHT FRACTION OF FUEL IN TOTAL FUELS AND OF", &
       & " OXIDANT IN TOTAL OXIDANTS")')
  go to 200
!***********************************************************************
  entry OUT4
  write(IOOUT, *)
  write(IOOUT, '(" TRANSPORT PROPERTIES (GASES ONLY)")')
  if (SIunit) then
     write(IOOUT, '("   CONDUCTIVITY IN UNITS OF MILLIWATTS/(CM)(K)"/)')
  else
     write(IOOUT, '("   CONDUCTIVITY IN UNITS OF MILLICALORIES/(CM)(K)(SEC)"/)')
  end if
! TRANSPORT PROPERTIES
  Fmt(4) = Fmt(6)
  if (Nplt > 0) then
     do i = 1, Npt
        if (i > ione) then
           if (mvis > 0) Pltout(i+Iplt-ione, mvis) = Vis(i)
           if (mcond > 0) Pltout(i+Iplt-ione, mcond) = Coneql(i)
           if (mpn > 0) Pltout(i+Iplt-ione, mpn) = Preql(i)
           if (mcondf > 0) Pltout(i+Iplt-ione, mcondf) = Confro(i)
           if (mpnf > 0) Pltout(i+Iplt-ione, mpnf) = Prfro(i)
        end if
     end do
  end if
  call VARFMT(Vis)
  write(IOOUT, Fmt) 'VISC,MILLIPOISE', (Vis(j), j=1, Npt)
  Fmt(4) = '13'
  Fmt(5) = ' '
  Fmt(7) = '4,'
  if (Eql) then
     write(IOOUT, '(/"  WITH EQUILIBRIUM REACTIONS"/)')
! SPECIFIC HEAT
     write(IOOUT, Fmt) fc, (Cpeql(j), j=1, Npt)
! CONDUCTIVITY
     write(IOOUT, Fmt) 'CONDUCTIVITY    ', (Coneql(j), j=1, Npt)
! PRANDTL NUMBER
     write(IOOUT, Fmt) 'PRANDTL NUMBER  ', (Preql(j), j=1, Npt)
  end if
  write(IOOUT, '(/"  WITH FROZEN REACTIONS"/)')
! SPECIFIC HEAT
  write(IOOUT, Fmt) fc, (Cpfro(j), j=1, Npt)
! CONDUCTIVITY
  write(IOOUT, Fmt) 'CONDUCTIVITY    ', (Confro(j), j=1, Npt)
! PRANDTL NUMBER
  write(IOOUT, Fmt) 'PRANDTL NUMBER  ', (Prfro(j), j=1, Npt)
200 return
end subroutine



subroutine REACT
!***********************************************************************
! READ AND PROCESS REACTANT RECORDS.  CALLED FROM subroutine INPUT.
!***********************************************************************
  use cea
  implicit none
! LOCAL VARIABLES
  character(6), save:: date
  character(2), save:: el(5)
  character(15), save:: sub
  integer, save:: i, icf, ifaz, ifrmla, itot, j, jj, k, kk, kr, l, n, nall, nint, nj, ntgas, ntot
  logical, save:: fuel, rcoefs, wdone(2)
  logical:: hOK
  real(8), save:: bb(5), dat(35), dift, eform, pcwt, rcf(9, 3), rm, t1, t2
  real(8):: t1save, t2save

  do k = 1, 2
     wdone(k) = .false.
     Wp(k) = 0.
     Hpp(k) = 0.
     Vpls(k) = 0.
     Vmin(k) = 0.
     Am(k) = 0.
     Rh(k) = 0.
     do j = 1, maxEl
        Elmt(j) = ' '
        B0p(j, k) = 0.
     end do
  end do
  do i = 1, maxEl
     dat(i) = 0.
  end do
! IF OXIDANT, KR = 1
! IF FUEL, KR = 2
  do n = 1, Nreac
     hOK = .false.
     t1save = 20000.d0
     t2save = 0.d0
     rcoefs = .true.
     if (Energy(n) == 'lib' .or. Rnum(n, 1) == 0.) then
        Tt = Rtemp(n)
        rewind IOTHM
        read(IOTHM) Tg, ntgas, ntot, nall
        do 20 itot = 1, nall
           if (itot <= ntot) then
              icf = 3
              if (itot > ntgas) icf = 1
              read(IOTHM) sub, nint, date, (el(j), bb(j), j=1, 5), ifaz, t1, t2, &
                   rm, ((rcf(i, j), i=1, 9), j=1, icf)
           else
              read(IOTHM) sub, nint, date, (el(j), bb(j), j=1, 5), ifaz, t1, t2, &
                   rm, eform
              if (nint > 0) read(IOTHM) ((rcf(i, j), i=1, 9), j=1, nint)
           end if
           if (sub == Rname(n) .or. sub == '*'//Rname(n)) then
              if (nint == 0) then
                 rcoefs = .false.
                 hOK = .true.
                 Enth(n) = eform*1000.d0/Rr
                 if (Tt == 0) then
                    Tt = t1
                    Rtemp(n) = t1
                 else
                    dift = abs(Tt-t1)
                    if (dift > 01d0) then
                       if (dift > 10.d0) then
                          write(IOOUT, '(/" REACTANT ", A15, "HAS BEEN DEFINED FOR THE TEMPERATURE", &
                               & F8.2, "K ONLY."/" YOUR TEMPERATURE ASSIGNMENT", F8.2, &
                               & " IS MORE THAN 10 K FROM THIS VALUE. (REACT)")') Rname(n), t1, Tt
                          Nlm = 0
                          hOK = .false.
                          go to 200
                       else
                          write(IOOUT, '(/" NOTE! REACTANT ", A15, "HAS BEEN DEFINED FOR ", &
                               & "TEMPERATURE", F8.2, "K ONLY."/" YOUR TEMPERATURE ASSIGNMENT", &
                               & F8.2, " IS NOT = BUT <10 K FROM THIS VALUE. (REACT)")') Rname(n), t1, Tt
                          Tt = t1
                          Rtemp(n) = t1
                       end if
                    end if
                 end if
              else
                 if (ifaz <= 0) then
                    t1save = min(t1save, .8d0*tg(1))
                    t2save = max(t2save, 1.2d0*t2)
                 else
                    t1save = min(t1save, t1-.001d0)
                    t2save = max(t2save, t2+.001d0)
                 endif
                 if (t1save < Tt .and. Tt < t2save) hOK = .true.
              end if
              do j = 1, 5
                 if (bb(j) == 0.) go to 5
                 Nfla(n) = j
                 Ratom(n, j) = el(j)
                 Rnum(n, j) = bb(j)
              end do
5             if (Tt == 0.) then
                 if (.not. Hp) go to 50
                 write(IOOUT, '(/" TEMPERATURE MISSING FOR REACTANT NO.", I2, "(REACT)")') n
                 Nlm = 0
                 go to 200
              end if
              if (rcoefs .and. hOK) then
                 Tln = log(Tt)
                 l = 1
                 if (ifaz <= 0) then
                    if (Tt > Tg(2)) l = 2
                    if (Tt > Tg(3)) l = 3
                 end if
                 Enth(n) = (((((rcf(7, l)/5.d0)*Tt+rcf(6, l)/4.d0)*Tt+rcf(5 &
                      , l)/3.d0)*Tt+rcf(4, l)/2.d0)*Tt+rcf(3, l)) &
                      *Tt - rcf(1, l)/Tt + rcf(2, l)*Tln + rcf(8, l)
                 if (Vol .and. ifaz <= 0) Enth(n) = Enth(n) - Tt
              end if
              if (hOK) go to 50
           end if
20      continue
        if (.not. hOK) then
           write(IOOUT, '(/" YOUR ASSIGNED TEMPERATURE", F8.2, "K FOR ", A15, /, &
                & "IS OUTSIDE ITS TEMPERATURE RANGE", F8.2, " TO", F9.2, "K (REACT)")') Tt, Rname(n), t1save, t2save
           Energy(n) = ' '
           Nlm = 0
           goto 200
        endif
     end if
50   ifrmla = Nfla(n)
     if (Fox(n)(:1) == 'f') then
        fuel = .true.
        kr = 2
        Fox(n) = 'FUEL'
     else if (Fox(n)(:4) == 'name') then
        fuel = .true.
        kr = 2
        Fox(n) = 'NAME'
     else
        kr = 1
        Fox(n) = 'OXIDANT'
     end if
     do j = 1, maxEl
        dat(j) = 0.
     end do
! STORE ATOMIC SYMBOLS IN ELMT ARRAY.
! CALCULATE MOLECULAR WEIGHT.
! TEMPORARILY STORE ATOMIC VALENCE IN X.
     rm = 0.d0
     do 100 jj = 1, ifrmla
        do j = 1, maxEl
           nj = j
           if (Elmt(j) == ' ') go to 60
           if (Ratom(n, jj) == Elmt(j)) go to 80
        end do
60      Nlm = nj
        Elmt(j) = Ratom(n, jj)
80      do kk = 1, 100
           if (Symbol(kk) == Ratom(n, jj)) then
              rm = rm + Rnum(n, jj)*Atmwt(kk)
              Atwt(j) = Atmwt(kk)
              X(j) = Valnce(kk)
              dat(j) = dat(j) + Rnum(n, jj)
              go to 100
           end if
        end do
        write(IOOUT, '(/1x, a2, " NOT FOUND IN BLOCKDATA (REACT)")') Ratom(n, jj)
        Nlm = 0
        go to 200
100  continue
     if (Pecwt(n) < 0.) then
        Pecwt(n) = 0.
        if (.not. Moles .and. .not. wdone(kr)) then
           wdone(kr) = .true.
           Pecwt(n) = 100.
           write(IOOUT, '(/" WARNING!!  AMOUNT MISSING FOR REACTANT", I3, ".", &
                & /" PROGRAM SETS WEIGHT PERCENT = 100. (REACT)")') n
        else
           write(IOOUT, '(/" AMOUNT MISSING FOR REACTANT NO.", I2, "(REACT)")') n
           Nlm = 0
           go to 200
        end if
     end if
! ADD CONTRIBUTIONS TO WP(K), HPP(K), AM(K), AND B0P(I, K)
     if (Pecwt(n) > 0.) wdone(kr) = .true.
     pcwt = Pecwt(n)
     if (Moles) pcwt = pcwt*rm
     Wp(kr) = Wp(kr) + pcwt
     if (rm <= 0.d0) then
        Nlm = 0
        go to 200
     else
        Hpp(kr) = Hpp(kr) + Enth(n)*pcwt/rm
        Am(kr) = Am(kr) + pcwt/rm
        if (Dens(n) /= 0.) then
           Rh(kr) = Rh(kr) + pcwt/Dens(n)
        else
           Rh(1) = 0.
           Rh(2) = 0.
        end if
        do j = 1, Nlm
           B0p(j, kr) = dat(j)*pcwt/rm + B0p(j, kr)
        end do
        Rmw(n) = rm
     end if
  end do
  if (.not. fuel) then
! 100 PERCENT OXIDANT, SWITCH INDICES
     do n = 1, Nreac
        Fox(n) = ' '
     end do
     Wp(2) = Wp(1)
     Wp(1) = 0.
     Hpp(2) = Hpp(1)
     Am(2) = Am(1)
     Am(1) = 0.
     do j = 1, Nlm
        B0p(j, 2) = B0p(j, 1)
     end do
  end if
  if (Nlm /= 0) then
! NORMALIZE HPP(KKR), AM(KR), B0P(I, KR), AND PECWT(N).
! CALCULATE V+(KR), AND V-(KR)
     do kr = 1, 2
        if (Wp(kr) /= 0.) then
           Hpp(kr) = Hpp(kr)/Wp(kr)
           Am(kr) = Wp(kr)/Am(kr)
           if (Rh(kr) /= 0.) Rh(kr) = Wp(kr)/Rh(kr)
           do j = 1, Nlm
              B0p(j, kr) = B0p(j, kr)/Wp(kr)
              if (X(j) < 0.) Vmin(kr) = Vmin(kr) + B0p(j, kr)*X(j)
              if (X(j) > 0.) Vpls(kr) = Vpls(kr) + B0p(j, kr)*X(j)
           end do
           if (.not. Moles) then
              do n = 1, Nreac
                 if (Fox(n)(:1) /= 'O' .or. kr /= 2) then
                    if (Fox(n)(:1) == 'O' .or. kr /= 1) Pecwt(n) &
                         = Pecwt(n)/Wp(kr)
                 end if
              end do
           end if
        end if
     end do
     if (.not. Short) then
        if (Moles) then
           write(IOOUT, '(/4x, "REACTANT", 10x, A7, 3X, "(ENERGY/R),K", 3X, &
                & "TEMP,K  DENSITY"/, 8x, "EXPLODED FORMULA")') ' MOLES '
        else
           write(IOOUT, '(/4x, "REACTANT", 10x, A7, 3X, "(ENERGY/R),K", 3X, &
                & "TEMP,K  DENSITY"/, 8x, "EXPLODED FORMULA")') 'WT.FRAC'
        end if
        do n = 1, Nreac
           write(IOOUT, '(1x, a1, ": ", a15, f10.6, e15.6, f9.2, f8.4, /8x, 5(2x, a2, f8.5))') &
                Fox(n), Rname(n), Pecwt(n), Enth(n), &
                Rtemp(n), Dens(n), (Ratom(n, i), Rnum(n, i), i=1, Nfla(n))
        end do
     end if
  end if
200 return
end subroutine



subroutine RKTOUT
!***********************************************************************
! SPECIAL OUTPUT FOR ROCKET PROBLEMS.
!***********************************************************************
  use cea
  implicit none
! LOCAL VARIABLES
  character(4):: exit(11) = 'EXIT'
  character(15), save:: fi, fiv, fr, z(4)
  integer, save:: i, i23, i46, i57, i68, i79, ione, ixfr, ixfz, j, k, line, ln, mae, mcf, &
       misp, mivac, mmach, mppf, mppj, nex
  real(8), save:: agv, aw, gc, tem, tra, vaci(Ncol), ww


  if (.not. Eql) then
     write(IOOUT, '(/////10x, " THEORETICAL ROCKET PERFORMANCE ASSUMING FROZEN COMPOSITION")')
     if (Nfz > 1) write(IOOUT, '(33X, "AFTER POINT", I2)') Nfz
  else
     write(IOOUT, '(/////13x, " THEORETICAL ROCKET PERFORMANCE ASSUMING EQUILIBRIUM")')
     if (Iopt /= 0) write(IOOUT, '(/11x, " COMPOSITION DURING EXPANSION FROM FINITE AREA COMBUSTOR")')
     if (Iopt == 0) write(IOOUT, '(/10x, " COMPOSITION DURING EXPANSION FROM INFINITE AREA COMBUSTOR")')
  end if
  if (Ttt(1) == T(It)) write(IOOUT, '(25X, "AT AN ASSIGNED TEMPERATURE  ")')
  tem = Ppp(1)*14.696006d0/1.01325d0
  write(IOOUT, '(/1x, A3, " =", F8.1, " PSIA")') 'Pin', tem
  i23 = 2
  if (Iopt > 0) then
     if (Iopt == 1) write(IOOUT, '(" Ac/At =", F8.4, 6x, "Pinj/Pinf =", F10.6)') Subar(1), App(2)
     if (Iopt == 2) write(IOOUT, '(" MDOT/Ac =", F10.3, " (KG/S)/M**2", 6x, "Pinj/Pinf =", F10.6)') Ma, App(2)
     i23 = 3
  end if
  call OUT1
  Fmt(4) = Fmt(6)
  nex = Npt - 2
  if (Page1) then
     ione = 0
     i46 = 4
     i57 = 5
     i68 = 6
     i79 = 7
  else
     ione = i23
  end if
! PRESSURE RATIOS
  if (Iopt == 0) then
     write(IOOUT, '(/17X, "CHAMBER   THROAT", 11(5X, A4))') (exit(i), i=1, nex)
     call VARFMT(App)
     write(IOOUT, Fmt) 'Pinf/P         ', (App(j), j=1, Npt)
  else
     nex = nex - 1
     write(IOOUT, '(/, 17X, "INJECTOR  COMB END  THROAT", 10(5X, A4))') (exit(i), i=1, nex)
     X(1) = 1.d0
     do i = 2, Npt
        X(i) = Ppp(1)/Ppp(i)
     end do
     call VARFMT(X)
     write(IOOUT, Fmt) 'Pinj/P         ', (X(i), i=1, Npt)
  end if
  call OUT2
  mppf  = 0
  mppj  = 0
  mmach = 0
  mae   = 0
  mcf   = 0
  mivac = 0
  misp  = 0
  do 100 i = 1, Nplt
     ixfz = index(Pltvar(i)(2:), 'fz')
     ixfr = index(Pltvar(i)(2:), 'fr')
     if (ixfz /= 0 .or. ixfr /= 0) then
        if (Eql) go to 100
     else if (.not. Eql) then
        go to 100
     end if
     if (Pltvar(i)(:4) == 'pi/p' .or. Pltvar(i)(:3) == 'pip') then
        if (Iopt == 0) mppf = i
        if (Iopt /= 0) mppj = i
     else if (Pltvar(i)(:4) == 'mach') then
        mmach = i
     else if (Pltvar(i)(:2) == 'ae') then
        mae = i
     else if (Pltvar(i)(:2) == 'cf') then
        mcf = i
     else if (Pltvar(i)(:4) == 'ivac') then
        mivac = i
     else if (Pltvar(i)(:3) == 'isp') then
        misp = i
     end if
100 continue
  if (SIunit) then
     agv = 1.
     gc = 1.
     fr = 'CSTAR, M/SEC'
     fiv = 'Ivac, M/SEC'
     fi = 'Isp, M/SEC'
  else
     gc = 32.174
     agv = 9.80665
     fr = 'CSTAR, FT/SEC'
     fiv = 'Ivac,LB-SEC/LB'
     fi = 'Isp, LB-SEC/LB'
  end if
  do k = 2, Npt
     Spim(k) = (2.*Rr*(Hsum(1)-Hsum(k)))**.5/agv
! AW IS THE LEFT SIDE OF EQ.(6.12) IN RP-1311, PT I.
     aw = Rr*Ttt(k)/(Ppp(k)*Wm(k)*Spim(k)*agv**2)
     if (k == i23) then
        if (Iopt == 0) Cstr = gc*Ppp(1)*aw
        if (Iopt /= 0) Cstr = gc*Ppp(1)/App(2)*aw
     end if
     vaci(k) = Spim(k) + Ppp(k)*aw
     Vmoc(k) = 0.
     if (Sonvel(k) /= 0.) Vmoc(k) = Spim(k)*agv/Sonvel(k)
  end do
! MACH NUMBER
  Vmoc(1) = 0.
  if (Gammas(i23) == 0.) Vmoc(i23) = 0.
  Fmt(7) = '3,'
  write(IOOUT, Fmt) 'MACH NUMBER    ', (Vmoc(j), j=1, Npt)
  if (Trnspt) call OUT4
  write(IOOUT, '(/" PERFORMANCE PARAMETERS"/)')
! AREA RATIO
  Fmt(4) = '9x,'
  Fmt(i46) = '9x,'
  call VARFMT(Aeat)
  Fmt(5) = ' '
  Fmt(i57) = ' '
  write(IOOUT, Fmt) 'Ae/At          ', (Aeat(j), j=2, Npt)
! C*
  Fmt(i57) = '13'
  Fmt(i68) = Fmt(i68+2)
  Fmt(i79) = '1,'
  write(IOOUT, Fmt) fr, (Cstr, j=2, Npt)
! CF - THRUST COEFICIENT
  Fmt(i79) = '4,'
  do i = 2, Npt
     X(i) = gc*Spim(i)/Cstr
  end do
  write(IOOUT, Fmt) 'CF             ', (X(j), j=2, Npt)
! VACUUM IMPULSE
  Fmt(i57) = '13'
  Fmt(i79) = '1,'
  write(IOOUT, Fmt) fiv, (vaci(j), j=2, Npt)
! SPECIFIC IMPULSE
  write(IOOUT, Fmt) fi, (Spim(j), j=2, Npt)
  if (Nplt > 0) then
     Spim(1) = 0
     Aeat(1) = 0
     Vmoc(1) = 0
     vaci(1) = 0
     X(1) = 0
     Spim(1) = 0
     do i = ione + 1, Npt
        if (mppj > 0) Pltout(i+Iplt-ione, mppj) = Ppp(1)/Ppp(i)
        if (mppf > 0) Pltout(i+Iplt-ione, mppf) = App(i)
        if (mmach > 0) Pltout(i+Iplt-ione, mmach) = Vmoc(i)
        if (mae > 0) Pltout(i+Iplt-ione, mae) = Aeat(i)
        if (mcf > 0) Pltout(i+Iplt-ione, mcf) = X(i)
        if (mivac > 0) Pltout(i+Iplt-ione, mivac) = vaci(i)
        if (misp > 0) Pltout(i+Iplt-ione, misp) = Spim(i)
     end do
  end if
  write(IOOUT, *)
  Fmt(4) = ' '
  Fmt(5) = '13'
  Fmt(7) = '5,'
  if (Iopt /= 0) then
     Fmt(i46) = Fmt(8)
     Fmt(i57) = Fmt(9)
  end if
  if (.not. Eql) then
     if (Massf) then
        write(IOOUT, '(1x, A4, " FRACTIONS"/)') 'MASS'
     else
        write(IOOUT, '(1x, A4, " FRACTIONS"/)') 'MOLE'
        ww = 1.d0/Totn(Nfz)
     end if
! MOLE (OR MASS) FRACTIONS - FROZEN
     tra = 5.E-6
     if (Trace /= 0.) tra = Trace
     line = 0
     do k = 1, Ngc
        if (Massf) ww = Mw(k)
        X(line+1) = En(k, Nfz)*ww
        if (X(line+1) >= tra) then
           line = line + 1
           z(line) = Prod(k)
        end if
        if (line == 3 .or. k == Ngc) then
           if (line == 0) go to 200
           write(IOOUT, '(1X, 3(A15, F8.5, 3X))') (z(ln), X(ln), ln=1, line)
           line = 0
        end if
     end do
  end if
200 call OUT3
  return
end subroutine



subroutine ROCKET
!***********************************************************************
! EXECUTIVE ROUTINE FOR ROCKET PROBLEMS.
!***********************************************************************
  use cea
  implicit none
! LOCAL VARIABLES
  integer, save:: i, i01, i12, iof, iplt1, iplte, ipp, isub, isup1, isupsv, itnum, &
       itrot, nar, nipp, niter, nn, npr1, nptth
  logical, save:: done, seql, thi
  real(8):: a1l = -1.26505, b1 = 1.0257, c1 = -1.2318, pa = 1e5
  real(8), save:: acatsv, aeatl, appl, aratio, asq, check, cprf, dd, dh, &
       dlnp, dlnpe, dlt, dp, eln, mat, msq, p1, pcpa, pcplt, pinf, pinj, &
       pinjas, pjrat, ppa, pr, pracat, prat, pratsv, pvg, test, tmelt, usq

  iplte = Iplt
  isup1 = 1
  App(1) = 1.
  Iopt = 0
  Npp = Npp + 2
  nn = Npp
  i01 = 0
  i12 = 1
  nipp = 1
  nptth = 2
  if (Fac) then
     Eql = .true.
     Npp = Npp + 1
     if (Acat /= 0.) then
        Iopt = 1
     else if (Ma /= 0.) then
        Iopt = 2
     else
        write(IOOUT, '(/" FATAL ERROR!! EITHER mdot OR ac/at MISSING FOR fac PROBLEM (ROCKET)")')
        Tt = 0.
        go to 1400
     end if
     i01 = 1
     i12 = 2
     nipp = 2
     nptth = 3
     do i = Nsub, 1, - 1
        Subar(i+1) = Subar(i)
     end do
     Nsub = Nsub + 1
     if (Iopt /= 1) then
        if (Acat == 0.) Acat = 2.
     end if
     Subar(1) = Acat
  else if (.not. Eql .and. Nfz > 1 .and. Nsub > 0) then
     Nsub = 0
     write(IOOUT, '(/" WARNING!!  FOR FROZEN PERFORMANCE, SUBSONIC AREA ", /, &
          & " RATIOS WERE OMITTED SINCE nfz IS GREATER THAN 1 (ROCKET)")')
  end if
  nn = nn + Nsub + Nsup
  if (Nfz > 2 .and. nn > Ncol-2) then
     write(IOOUT, '(/" WARNING!!  nfz NOT ALLOWED TO BE > 2 IF THE TOTAL", /, &
          & " NUMBER OF POINTS IS >", i3, " (ROCKET)")') Ncol - 2
     Nfz = 1
     Froz = .false.
  end if
  seql = Eql
  iof = 0
  Tt = Tcest
  Pp = P(1)
  App(i12) = 1.
! LOOP FOR EACH O/F
100 It = 1
  iof = iof + 1
  Oxfl = Oxf(iof)
  if (T(1) /= 0.) then
     Tp = .true.
  else
     Hp = .true.
  end if
  Sp = .false.
  call NEWOF
  if (T(1) /= 0.) Tt = T(1)
! LOOP FOR CHAMBER PRESSURES
200 do Ip = 1, Np
     itnum = 0
     Area = .false.
     if (T(1) == 0.) Hp = .true.
     if (T(1) /= 0.) Tp = .true.
     Sp = .false.
     Eql = .true.
     isub = 1
     Isup = 1
     Pp = P(Ip)
     pinf = Pp
     ipp = 1
     itrot = 3
     isupsv = 1
     niter = 1
     Page1 = .true.
     iplt1 = iplte
     Iplt = iplte
     done = .false.
! LOOP FOR OUTPUT COLUMNS
250  nar = Npt
     if (Eql) then
        call EQLBRM
        if (Npt == Nfz) cprf = Cpsum
     else
        call FROZEN
     end if
! TT = 0 IF NO CONVERGENCE
     if (Tt /= 0.) then
! TEST FOR FINITE AREA COMBUSTOR
        if (.not. Fac) go to 400
        pinjas = P(Ip)*pa
        pinj = pinjas
        if (Npt <= 2) then
           if (Npt == 1 .and. Trnspt) call TRANP
           if (Npt == 2) pinf = Ppp(2)
        end if
        if (Npt /= 1) go to 400
! INITIAL ESTIMATE FOR PC (AND ACAT IF NOT ASSIGNED)
        do i = 1, 4
           prat = (b1+c1*Acat)/(1.+a1l*Acat)
           ppa = pinj*prat
           if (Iopt == 1) go to 260
           Acat = ppa/(Ma*2350.)
           if (Acat >= 1.) then
              pratsv = prat
              if (Debugf) then
                 if (i <= 1) write(IOOUT, '(/"  ITERATION", 9X, "PC", 7X, "CONTRACTION RATIO")')
                 write(IOOUT, '(5X, I2, 7X, F12.2, 3X, F12.6)') i, ppa, Acat
              end if
           else
              write(IOOUT, '(/" INPUT VALUE OF mdot/a =", F12.3, " IS TOO LARGE."/ &
                   & " GIVES CONTRACTION RATIO ESTIMATE LESS THAN 1 (ROCKET)")') Ma
              Tt = 0.
              go to 1400
           end if
        end do
        Subar(1) = Acat
260     Pp = ppa/pa
        App(1) = Pp/Ppp(1)
        go to 1100
     else
        if (Npt < 1) go to 1400
        if (.not. Area) go to 600
        Npt = nar - 1
        Isup = Nsup + 2
        Isv = 0
        itnum = 0
        go to 950
     end if
300  Hp = .true.
     Sp = .false.
     niter = niter + 1
     Isv = 0
     Npt = 2
     ipp = 2
     call SETEN
     go to 250
350  done = .true.
     App(1) = Ppp(2)/Ppp(1)
     Area = .false.
     if (Nsub > 1) isub = 2
     Isv = 4
     Npt = 2
     ipp = min(4, Npp)
     call SETEN
     Cpr(2) = Cpr(4)
     Dlvpt(2) = Dlvpt(4)
     Dlvtp(2) = Dlvtp(4)
     Gammas(2) = Gammas(4)
     Hsum(2) = Hsum(4)
     Ppp(2) = Ppp(4)
     App(2) = Ppp(1)/pinf
     Ssum(2) = Ssum(4)
     Totn(2) = Totn(4)
     Ttt(2) = Ttt(4)
     Vlm(2) = Vlm(4)
     Wm(2) = Wm(4)
     if (.not. Short) write(IOOUT, '(" END OF CHAMBER ITERATIONS")')
     go to 600
! INITIALIZE FOR THROAT
400  if (ipp > nipp) then
        usq = 2.*(Hsum(1)-Hsum(Npt))*Rr
        if (ipp > nptth) go to 600
! THROAT
        if (.not. thi) then
           Vv = Vlm(nptth)
           pvg = Pp*Vv*Gammas(nptth)
           if (pvg == 0.) then
              write(IOOUT, '(/" WARNING!!  DIFFICULTY IN LOCATING THROAT (ROCKET)")')
              go to 550
           else
              msq = usq/pvg
              if (Debug(1) .or. Debug(2)) write(IOOUT, '(/" USQ=", E15.8, 5X, "PVG=", E15.8)') usq, pvg
              dh = abs(msq-1.d0)
              if (dh <= 0.4d-4) go to 550
              if (itrot > 0) then
                 p1 = Pp
                 if (Jsol /= 0) then
                    tmelt = Tt
                    Pp = Pp*(1.d0+msq*Gammas(nptth))/(Gammas(nptth)+1.d0)
                 else if (tmelt == 0.) then
                    Pp = Pp*(1.d0+msq*Gammas(nptth))/(Gammas(nptth)+1.d0)
                 else
                    write(IOOUT, '(/" WARNING!!  DISCONTINUITY AT THE THROAT (ROCKET)")')
                    dlt = log(tmelt/Tt)
                    dd = dlt*Cpr(nptth)/(Enn*Dlvtp(nptth))
                    Pp = Pp*EXP(dd)
                    App(nptth) = P(Ip)/Pp
                    if (Fac) App(nptth) = pinf/Pp
                    if (Eql .and. .not. Short) write(IOOUT, '(" Pinf/Pt =", F9.6)') App(nptth)
                    thi = .true.
                    go to 250
                 end if
                 go to 500
              else if (itrot < 0) then
                 if (itrot < -19) then
                    write(IOOUT, '(/" WARNING!!  DIFFICULTY IN LOCATING THROAT (ROCKET)")')
                    go to 550
                 else
                    if (Npr /= npr1) go to 550
                    Pp = Pp - dp
                    go to 500
                 end if
              else if (Npr == npr1) then
                 write(IOOUT, '(/" WARNING!!  DIFFICULTY IN LOCATING THROAT (ROCKET)")')
                 go to 550
              else
                 dp = abs(Pp-p1)/20.
                 Pp = max(Pp, p1)
                 write(IOOUT, '(/" WARNING!!  DISCONTINUITY AT THE THROAT (ROCKET)")')
                 Pp = Pp - dp
                 go to 500
              end if
           end if
        else
           Gammas(nptth) = 0.
           go to 550
        end if
     else
        if (.not. Fac .and. Trnspt) call TRANP
        if (Npt == Nfz) Eql = seql
        Tp = .false.
        Hp = .false.
        Sp = .true.
        S0 = Ssum(i12)
     end if
450  tmelt = 0.
     itrot = 3
     thi = .false.
     App(nptth) = ((Gammas(i12)+1.)/2.) &
          **(Gammas(i12)/(Gammas(i12)-1.))
     if (Eql .and. .not. Short) write(IOOUT, '(" Pinf/Pt =", F9.6)') App(nptth)
     Pp = pinf/App(nptth)
     Isv = -i12
     go to 1200
500  npr1 = Npr
     App(nptth) = P(Ip)/Pp
     if (Fac) App(nptth) = pinf/Pp
     if (Eql .and. .not. Short) write(IOOUT, '(" Pinf/Pt =", F9.6)') App(nptth)
     itrot = itrot - 1
     go to 250
550  Awt = Enn*Tt/(Pp*usq**.5)
     pcplt = log(App(nptth))
600  Isv = 0
     Aeat(Npt) = Enn*Ttt(Npt)/(Pp*usq**.5*Awt)
     if (Tt == 0.) go to 1150
     if (Area) go to 750
     if (Trnspt .and. (.not. Fac .or. done .or. Npt > 2)) call TRANP
     if (Npt == Nfz) Eql = seql
     if (Fac) then
        if (Npt == nptth) then
           Area = .true.
           go to 750
        else if (Npt == 2 .and. done) then
           Npt = 3
!  The following statement was corrected 1/30/2004.  Only fac parameters 
!    after combustion were affected--generally extra or missing points.
!  (remove) if (ipp <= Npp) ipp = ipp - 1
           if (ipp < Npp .or. npp == 4) ipp = ipp - 1
        end if
     end if
650  if (ipp < Npp) go to 1100
700  if (Nsub == i01 .and. Nsup == 0) go to 1150
     Area = .true.
! PCP ESTIMATES FOR AREA RATIOS
750  if (itnum == 0) then
        dlnp = 1.
        itnum = 1
        aratio = Subar(isub)
        if ((.not. Fac .or. done) .and. Nsub <= i01) aratio = Supar(Isup)
        if (.not. Eql .and. Nfz >= 3) then
           if (aratio <= Aeat(Nfz)) then
              write(IOOUT, '(/, " WARNING!! FOR FROZEN PERFORMANCE, POINTS WERE OMITTED", &
                   & " WHERE THE ASSIGNED", /, " SUPERSONIC AREA RATIOS WERE ", &
                   & "LESS THAN THE VALUE AT POINT nfz =", I3, " (ROCKET)")') Nfz
              go to 1050
           end if
        end if
        if (aratio  <  1.d0) then
           write(IOOUT, '(/" AN ASSIGNED AREA RATIO IS < 1 (ROCKET)")')
           go to 1050
        end if
        eln = log(aratio)
        if (Fac) then
           if (.not. done) go to 800
        end if
        if (Nsub <= i01) then
           if (Nfz == ipp) isupsv = Isup
           if (Supar(Isup) < 2.) then
              appl = sqrt(eln*(1.535d0+3.294d0*eln)) + pcplt
              go to 1100
           else
              if (Isup > isup1 .and. Supar(Isup-1) >= 2.) go to 850
              appl = Gammas(nptth) + eln*1.4
              go to 1100
           end if
        end if
! TEST FOR CONVERGENCE ON AREA RATIO.
     else if (Gammas(Npt) > 0.) then
        check = .00004
        if (Debug(Npt)) write(IOOUT, '(/" ITER=", I2, 2X, "ASSIGNED AE/AT=", F14.7, 3X, "AE/AT=", F14.7, &
             & /, 2X, "PC/P=", F14.7, 2X, "DELTA LN PCP=", F14.7)') itnum, aratio, Aeat(Npt), &
             App(Npt), dlnp
        if (abs(Aeat(Npt)-aratio)/aratio <= check) go to 900
        if (ABS(dlnp) < .00004) go to 900
        aeatl = log(Aeat(Npt))
        itnum = itnum + 1
        if (itnum > 10) then
           write(IOOUT, '(/" WARNING!!  DID NOT CONVERGE FOR AREA RATIO =", F10.5, &
                & " (ROCKET)")') aratio
           go to 900
        else
! IMPROVED PCP ESTIMATES.
           asq = Gammas(Npt)*Enn*Rr*Tt
           dlnpe = Gammas(Npt)*usq/(usq-asq)
           go to 850
        end if
     else
        write(IOOUT, '(/" WARNING!!  AREA RATIO CALCULATION CANNOT BE DONE ", &
             & "BECAUSE GAMMAs", /, " CALCULATION IMPOSSIBLE. (ROCKET)")')
        Npt = Npt - 1
        if (Nsub <= 0) isup1 = 100
        if (Nsub < 0.) Nsup = Isup - 1
        if (Nsub > 0) Nsub = isub - 1
        go to 1000
     end if
800  appl = pcplt/(Subar(isub)+(10.587*eln**2+9.454)*eln)
     if (aratio < 1.09) appl = .9*appl
     if (aratio > 10.) appl = appl/aratio
     if (isub > 1 .or. Npt == Ncol) go to 1100
     go to 1200
850  dlnp = dlnpe*eln - dlnpe*aeatl
     appl = appl + dlnp
     if (itnum == 1) go to 1100
     if (appl < 0.) appl = .000001
     App(Npt) = EXP(appl)
     Pp = pinf/App(Npt)
     go to 250
! CONVERGENCE HAS BEEN REACHED FOR ASSIGNED AREA RATIO
900  Aeat(Npt) = aratio
     if (Fac) then
        if (.not. done) then
           if (Iopt == 1) then
! OPTION 1 FOR FINITE AREA COMBUSTOR. INPUT IS ASSIGNED INJECTOR
! PRESSURE AND CONTRACTION RATIO. IMPROVED ESTIMATE FOR PC
              Area = .false.
              itnum = 0
              ppa = Ppp(Npt)*pa
              pinj = ppa + 1.d05*usq/Vlm(Npt)
              test = (pinj-pinjas)/pinjas
              pcpa = pinf*pa
              if (Debugf) then
                 write(IOOUT, '(" ITER", 3X, "TEST", 3X, "ASSIGNED PINJ", 1x, "CALC PINJ", 5X, &
                      & "PC", 7X, "P AT ACAT", 3X, "PREV ACAT", 2X, "ACAT")')
                 write(IOOUT, '(I3, F10.6, 1x, 4F12.2, 2F9.5)') niter, test, pinjas, pinj, pcpa, ppa, &
                      acatsv, Acat
              end if
              if (ABS(test) < 0.00002) go to 350
              prat = pinjas/pinj
              Pp = pinf*prat
              go to 300
           else if (Iopt == 2) then
! OPTION 2 FOR FINITE AREA COMBUSTOR. INPUT IS ASSIGNED INJECTOR
! PRESSURE AND MASS FLOW PER UNIT AREA. IMPROVED ESTIMATE FOR PC
! AND ACAT
              acatsv = Acat
              pratsv = prat
              Area = .false.
              itnum = 0
              ppa = Ppp(4)*pa
              pinj = ppa + 1.d05*usq/Vlm(4)
              mat = pa/(Awt*Rr)
              Acat = mat/Ma
              prat = (b1+c1*Acat)/(1.+a1l*Acat)
              test = (pinj-pinjas)/pinjas
              pcpa = pinf*pa
              if (Debugf) then
                 write(IOOUT, '(" ITER", 3X, "TEST", 3X, "ASSIGNED PINJ", 1x, "CALC PINJ", 5X, &
                      & "PC", 7X, "P AT ACAT", 3X, "PREV ACAT", 2X, "ACAT")')
                 write(IOOUT, '(I3, F10.6, 1x, 4F12.2, 2F9.5)') niter, test, pinjas, pinj, pcpa, ppa, &
                      acatsv, Acat
              end if
              if (ABS(test) < 0.00002) go to 350
              pjrat = pinj/pinjas
              Pp = pinf
              do i = 1, 2
                 pracat = pratsv/prat
                 pr = pjrat*pracat
                 Pp = Pp/pr
                 pcpa = Pp*pa
                 Acat = Acat/pr
                 Subar(1) = Acat
                 pratsv = prat
                 pjrat = 1.
                 prat = (b1+c1*Acat)/(1.+a1l*Acat)
                 if (Debugf) write(IOOUT, '(" NEW PC = ", F10.2, 2X, "NEW ACAT = ", F9.6, 2X, "PJRAT =", &
                      & F10.7, " PRACAT =", F10.7)') pcpa, Acat, pjrat, pracat
              end do
              go to 300
           end if
        end if
     end if
950  if (Trnspt) call TRANP
     if (Npt == Nfz) Eql = seql
1000 itnum = 0
     if (Nsub > i01) then
        isub = isub + 1
        if (isub <= Nsub) go to 750
        isub = 1
        Nsub = -Nsub
        if (Isup <= Nsup) go to 750
        Area = .false.
        go to 1150
     end if
1050 Isup = Isup + 1
     itnum = 0
     if (Isup <= Nsup) go to 750
     Isup = isupsv
     Area = .false.
     go to 1150
! TEST FOR OUTPUT -- SCHEDULES COMPLETE OR NPT=Ncol
1100 Isv = Npt
     if (Npt /= Ncol) go to 1200
1150 if (.not. Eql) then
        if (Nfz <= 1) then
           Cpr(Nfz) = cprf
           Gammas(Nfz) = cprf/(cprf-1./Wm(Nfz))
        end if
     end if
     call RKTOUT
     Iplt = Iplt + Npt
     if (.not. Page1) then
        Iplt = Iplt - 2
        if (Iopt /= 0) Iplt = Iplt - 1
        Iplt = min(Iplt, 500)
     else
        Page1 = .false.
     end if
     iplte = max(iplte, Iplt)
     dlnp = 1.
     if (Tt == 0.) Area = .false.
     if (.not. Eql .and. Tt == 0.) write(IOOUT, '(/" WARNING!!  CALCULATIONS WERE STOPPED BECAUSE NEXT ", &
          & "POINT IS MORE", /, " THAN 50 K BELOW THE TEMPERATURE", &
          & " RANGE OF A CONDENSED SPECIES (ROCKET)")')
     if (Isv == 0) then
! PCP, SUBAR, AND SUPAR SCHEDULES COMPLETED
        if (Nsub < 0) Nsub = -Nsub
        if (.not. Froz .or. .not. Eql) go to 1300
! SET UP FOR FROZEN.
        if (Eql) Iplt = iplt1
        Eql = .false.
        Page1 = .true.
        call SETEN
        Tt = Ttt(Nfz)
        ipp = Nfz
        if (Nfz == Npt) go to 1150
        Npt = Nfz
        Enn = 1./Wm(Nfz)
        if (Nfz == 1) go to 450
        if (Nsub > 0) then
           Nsub = -Nsub
           write(IOOUT, '(/" WARNING!!  FOR FROZEN PERFORMANCE, SUBSONIC AREA ", /, &
                & " RATIOS WERE OMITTED SINCE nfz IS GREATER THAN 1 (ROCKET)")')
        end if
        if (App(Nfz) < App(nptth)) then
           write(IOOUT, '(/" WARNING!!  FREEZING IS NOT ALLOWED AT A SUBSONIC ", &
                & "PRESSURE RATIO FOR nfz GREATER"/" THAN 1. FROZEN ", &
                & "PERFORMANCE CALCULATIONS WERE OMITTED (ROCKET)")')
        else
           if (Nfz < Npp) go to 1200
           go to 700
        end if
        go to 1300
     else
        if (Eql) write(IOOUT, '(////)')
        Npt = nptth
     end if
! SET INDICES AND ESTIMATES FOR NEXT POINT.
1200 Npt = Npt + 1
     if (Eql .or. (Isv == -i12 .and. .not. seql)) then
! THE FOLLOWING STATEMENT WAS ADDED TO TAKE CARE OF A SITUATION
! WHERE EQLBRM WENT SINGULAR WHEN STARTING FROM ESTIMATES WHERE
! BOTH SOLID AND LIQUID WERE INCLUDED.  JULY 27, 1990.
        if (Jliq /= 0 .and. Isv > 0) Isv = 0
        call SETEN
     end if
1250 ipp = ipp + 1
     if (Npt > nptth) then
        if (Area) then
           App(Npt) = EXP(appl)
        else
           App(Npt) = Pcp(ipp-nptth)
           if (Fac) App(Npt) = App(Npt)*pinf/Ppp(1)
           if (.not. Eql .and. App(Npt) < App(Nfz)) then
              write(IOOUT, '(/, " WARNING!! FOR FROZEN PERFORMANCE, POINTS WERE OMITTED", &
                   & " WHERE THE ASSIGNED", /, &
                   & " PRESSURE RATIOS WERE LESS THAN ", &
                   & "THE VALUE AT POINT nfz =", I3, " (ROCKET)")') Nfz
              go to 1250
           end if
        end if
        Pp = pinf/App(Npt)
        if (Fac) then
           if (Area) then
              if (isub <= Nsub .and. isub > i01 .and. aratio >= Aeat(2)) &
                   then
                 write(IOOUT, '(/" WARNING!!  ASSIGNED subae/at =", f10.5, " IS NOT ", &
                      & "PERMITTED TO BE GREATER"/" THAN ac/at =", f9.5, &
                      & ".  POINT OMITTED (ROCKET)")') aratio, Aeat(2)
                 Npt = Npt - 1
                 go to 1000
              end if
           else if (Npt > nptth .and. Pcp(ipp-3) < Ppp(1)/Ppp(2)) then
              write(IOOUT, '(/" WARNING!!  ASSIGNED pip =", F10.5, &
                   & " IS NOT PERMITTED"/" TO BE LESS THAN  Pinj/Pc =", f9.5, &
                   & ". POINT OMITTED", " (ROCKET)")') Pcp(ipp-3), Ppp(1)/Ppp(2)
              Npt = Npt - 1
              go to 650
           end if
        end if
     end if
     go to 250
1300 Npt = 1
! CHECK FOR COMPLETED SCHEDULES -
! 1) CHAMBER PRESSURES(IP = NP)
! 2) CHAMBER TEMPERATURES(IT = NT)
! 3) O/F VALUES(IOF = NOF)
     if (Ip == Np .and. It == Nt .and. iof == Nof) go to 1400
     write(IOOUT, '(////)')
     call SETEN
     Tt = Ttt(i12)
  end do
  if (It < Nt) then
     It = It + 1
     Tt = T(It)
     go to 200
  else if (iof < Nof) then
     go to 100
  end if
1400 Iplt = max(Iplt, iplte)
  return
end subroutine ROCKET



subroutine SEARCH
!***********************************************************************
! SEARCH THERMO.LIB FOR THERMO DATA FOR SPECIES TO BE CONSIDERED.
!***********************************************************************
  use cea
  implicit none
! LOCAL VARIABLES
  character(16), save:: bin(2, 40), pure(6), spece(2)
  character(6), save:: date(maxNgc)
  character(2), save:: el(5)
  character(15), save:: sub
  integer, save:: i, i5, ifaz, ii, ir, itot, j, jj(2), jk, k, lineb, nall, ne, nint, &
       npure, nrec, ntgas, ntot
  real(8), save:: b(5), t1, t2, thermo(9, 3), trdata(36)

  Nc = 0
  ne = 0
  do i = 1, Nlm
     Jx(i) = 0
  end do
  do j = 1, maxNgc
     S(j) = 0.
     H0(j) = 0.
     Deln(j) = 0.
     do i = 1, Nlm
        A(i, j) = 0.
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
  rewind IOTHM
  read(IOTHM) Tg, ntgas, ntot, nall, Thdate
  Ngc = 1
  Nc = 1
! BEGIN LOOP FOR READING SPECIES DATA FROM THERMO.LIB.
  do 200 itot = 1, ntot
     if (itot > ntgas) then
        read(IOTHM) sub, nint, date(Ngc), (el(j), b(j), j=1, 5), Ifz(Nc), &
             Temp(1, Nc), Temp(2, Nc), Mw(Ngc), (Cft(Nc, k), k=1, 9)
     else
        read(IOTHM) sub, nint, date(Ngc), (el(j), b(j), j=1, 5), ifaz, t1, t2, &
             Mw(Ngc), thermo
     end if
     if (Nonly /= 0) then
        i = 1
20      if (Prod(i) /= sub .and. '*'//Prod(i) /= sub) then
           i = i + 1
           if (i <= Nonly) go to 20
           go to 200
        else
           if (sub == Prod(Ngc-1)) then
              Nonly = Nonly + 1
              do k = Nonly, i + 1, - 1
                 Prod(k) = Prod(k-1)
              end do
           else
              Prod(i) = Prod(Ngc)
           end if
           Prod(Ngc) = sub
        end if
     else if (Nomit /= 0) then
        do i = 1, Nomit
           if (Omit(i) == sub .or. '*'//Omit(i) == sub) go to 200
        end do
     end if
     do 50 k = 1, 5
        if (b(k) == 0.) go to 100
        do i = 1, Nlm
           if (Elmt(i) == el(k)) then
              A(i, Ngc) = b(k)
              go to 50
           end if
        end do
        do j = 1, Nlm
           A(j, Ngc) = 0.
        end do
        go to 200
50   continue
100  Prod(Ngc) = sub
     if (itot > ntgas) then
        Nc = Nc + 1
        if (Nc > maxNc) go to 400
     else
        Ng = Ngc
        if (Ng > maxNg) go to 400
        do i = 1, 3
           do j = 1, 9
              Coef(Ng, j, i) = thermo(j, i)
           end do
        end do
! IF SPECIES IS AN ATOMIC GAS, STORE INDEX IN JX
        if (b(2) == 0. .and. b(1) == 1.) then
           do i = 1, Nlm
              if (Elmt(i) == el(1)) then
                 ne = ne + 1
                 Jx(i) = Ngc
                 Jcm(i) = Ngc
                 go to 150
              end if
           end do
        end if
     end if
150  Ngc = Ngc + 1
     if (Ngc > maxNgc) go to 400
200 continue
! FINISHED READING THERMO DATA FROM I/O UNIT IOTHM.
  Ifz(Nc) = 0
  Nc = Nc - 1
  Ngc = Ngc - 1
  Ngp1 = Ng + 1
  if (Ngc < Nonly) then
     do k = Ngc + 1, Nonly
        write(IOOUT, '(/" WARNING!!  ", A15, " NOT A PRODUCT IN thermo.lib FILE (SEARCH)")') Prod(k)
     end do
  end if
! FIND MISSING ELEMENTS (IF ANY) FOR COMPONENTS
  Nspx = Ngc
  if (ne < Nlm) then
     do i = 1, Nlm
        if (Nspx > maxNgc) go to 400
        if (Jx(i) == 0) then
           Nspx = Nspx + 1
           do k = 1, Nlm
              A(k, Nspx) = 0.
           end do
           A(i, Nspx) = 1.
           Prod(Nspx) = Elmt(i)
           do k = 1, 100
              if (Elmt(i) == Symbol(k)) then
                 Mw(Nspx) = Atmwt(k)
                 Atwt(i) = Atmwt(k)
                 Cp(Nspx) = 2.5d0
                 go to 210
              end if
           end do
210        Jx(i) = Nspx
           Jcm(i) = Nspx
        end if
     end do
  end if
! ARE ALL ELEMENTS IN PRODUCT SPECIES?
  do 300 i = 1, Nlm
     do j = 1, Ngc
        if (A(i, j) /= 0.) go to 300
        ii = i
     end do
     write(IOOUT, '(/" PRODUCT SPECIES CONTAINING THE ELEMENT", A3, " MISSING", &
          & //, 13x, "FATAL ERROR (SEARCH)")') Elmt(ii)
     Ngc = 0
     go to 600
300 continue
! WRITE POSSIBLE PRODUCT LIST
  if (.not. Short) then
     write(IOOUT, '(/2x, "SPECIES BEING CONSIDERED IN THIS SYSTEM", &
          & /" (CONDENSED PHASE MAY HAVE NAME LISTED SEVERAL TIMES)", &
          & /"  LAST thermo.inp UPDATE: ", A10, /)') Thdate
     do i = 1, Ngc, 3
        i5 = i + 2
        if (Ngc < i5) i5 = Ngc
        write(IOOUT, '(3(2X, A6, 2X, A15))') (date(j), Prod(j), j=i, i5)
     end do
  end if
  go to 600
400 write(IOOUT, '(/" INSUFFICIENT STORAGE FOR PRODUCTS-SEE RP-1311,", &
       & /"   PART 2, PAGE 39. (SEARCH)")')
  Ngc = 0
  go to 600
! SEARCH FOR TRANSPORT PROPERTIES FOR THIS CHEMICAL SYSTEM
  entry READTR
  rewind IOTRN
  rewind IOSCH
  Ntape = 0
  npure = 0
  lineb = 1
  if (.not. Short) write(IOOUT, '(/" SPECIES WITH TRANSPORT PROPERTIES"//8X, "PURE SPECIES"/)')
  read(IOTRN) nrec
  do ir = 1, nrec
     read(IOTRN) spece, trdata
     k = 1
450  do j = 1, Ng
        if (spece(k) == Prod(j) .or. '*'//spece(k) == Prod(j)) then
           jj(k) = j
           if (k == 2) then
! STORE NAMES FOR BINARIES IN BIN ARRAY.
              do k = 1, 2
                 bin(k, lineb) = spece(k)
              end do
              lineb = lineb + 1
              go to 500
           else
              jj(2) = j
              if (spece(2) == ' ') then
! WRITE NAMES FOR PURE SPECIES.
                 npure = npure + 1
                 pure(npure) = spece(1)
                 go to 500
              else
                 k = 2
                 go to 450
              end if
           end if
        end if
     end do
     go to 550
500  write(IOSCH) jj, trdata
     Ntape = Ntape + 1
550  if (npure /= 0 .and. (npure >= 6 .or. ir >= nrec)) then
        if (.not. Short) write(IOOUT, '(4(2x, A16))') (pure(jk), jk=1, npure)
        npure = 0
     end if
  end do
  lineb = lineb - 1
  if (.not. Short) then
     write(IOOUT, '(/"     BINARY INTERACTIONS"/)')
     do j = 1, lineb
        write(IOOUT, '(5X, 2A16)') (bin(i, j), i=1, 2)
     end do
  end if
  write(IOOUT, *)
600 return
end subroutine



subroutine SETEN
!***********************************************************************
! USE COMPOSITIONS FROM PREVIOUS POINT AS INITIAL ESTIMATES FOR
! CURRENT POINT NPT.  IF -
!  ISV>0  USE COMPOSITIONS FROM POINT ISV.
!  ISV<0  SAVE COMPOSITIONS FROM POINT -ISV FOR POSSIBLE LATER USE.
!         ALSO USE COMPOSITIONS FROM POINT -ISV FOR NPT.
!  ISV=0  USE COMPOSITIONS SAVED WHEN ISV<0.
!***********************************************************************
  use cea
  implicit none
! LOCAL VARIABLES
  integer, save:: j, lsav
  real(8), save:: tsave

  if (Isv < 0) then
! FIRST T--SAVE COMPOSITIONS FOR FUTURE POINTS WITH THIS T
     Isv = -Isv
     tsave = Ttt(Isv)
     Ensave = Enn
     Enlsav = Ennl
     lsav = Lsave
     do j = 1, Ng
        Sln(j) = Enln(j)
     end do
     do j = 1, Ng
        En(j, Npt) = En(j, Isv)
     end do
     Npr = 0
     do j = Ngp1, Ngc
        Sln(j) = En(j, Isv)
        En(j, Npt) = Sln(j)
        if (Jliq == j) then
           En(Jsol, Npt) = En(Jsol, Isv) + En(Jliq, Isv)
           En(Jliq, Npt) = 0.
           Jsol = 0
           Jliq = 0
           tsave = tsave - 5.
           Tt = tsave
           Sln(j) = 0.
        else if (En(j, Npt) > 0.) then
           Npr = Npr + 1
           Jcond(Npr) = j
        end if
     end do
  else if (Isv == 0) then
! NEXT POINT FIRST T IN SCHEDULE, USE PREVIOUS COMPOSITIONS FOR THIS T
     Jsol = 0
     Jliq = 0
     Enn = Ensave
     Ennl = Enlsav
     Lsave = lsav
     Npr = 0
     do j = Ngp1, Ngc
        En(j, Npt) = Sln(j)
        if (En(j, Npt) > 0.d0) then
           Npr = Npr + 1
           Jcond(Npr) = j
        end if
     end do
     do j = 1, Ng
        En(j, Npt) = 0.
        Enln(j) = Sln(j)
        if (Sln(j) /= 0.) then
           if ((Enln(j)-Ennl+18.5) > 0.) En(j, Npt) = exp(Enln(j))
        end if
     end do
     if (.not. Tp) Tt = tsave
     Sumn = Enn
  else if (Isv > 0) then
! USE COMPOSITIONS FROM PREVIOUS POINT
     do j = 1, Ngc
        En(j, Npt) = En(j, Isv)
     end do
  end if
end subroutine SETEN



subroutine SHCK
!***********************************************************************
! PRIMARY ROUTINE FOR SHOCK PROBLEMS.
!***********************************************************************
  use cea
  implicit none
! LOCAL VARIABLES
  character(1), save:: cr12, cr52
  integer, save:: i, iof, it1, it2, itr, j, n
  logical, save:: refl, seql, srefl
  real(8), save:: ax, axx, b2, cormax, gg, hs, m2m1(Ncol), mis(13), mu12rt, p1, p21, &
       p21l, p2p1(Ncol), pmn, rho12, rho52, rrho(Ncol), sg(78), t1, t21, &
       t21l, t2t1(Ncol), ttmax, u1u2(Ncol), uis(13), utwo(Ncol), uu, wmx, ww

  if (Trace == 0.) Trace = 5.E-9
  Tp = .true.
  Cpmix = 0.
  srefl = .false.
  if (.not. Short) then
     write(IOOUT, '(/"   *** INPUT FOR SHOCK PROBLEMS ***")')
     write(IOOUT, '(/" INCDEQ =", L2, "   REFLEQ =", L2, "   INCDFZ =", L2, &
          & "    REFLFZ =", L2)') Incdeq, Refleq, Incdfz, Reflfz
  end if
  if (Refleq .or. Reflfz) srefl = .true.
  seql = Incdeq
  if (T(1) == 0.) T(1) = Rtemp(1)
  do i = 1, Nsk
     uis(i) = U1(i)
     mis(i) = Mach1(i)
     if (Mach1(i) == 0.0 .and. U1(i) == 0.0) go to 100
  end do
100 if (Nsk > Ncol) then
     write(IOOUT, '(/" WARNING!!  ONLY ", I2, " u1 OR mach1 VALUES ALLOWED (SHCK)")') Ncol
     Nsk = Ncol
  end if
  if (.not. Short) then
     write(IOOUT, '(/1p, " U1 =   ", 5E13.6, /(8X, 5E13.6))') (U1(i), i=1, Nsk)
     write(IOOUT, '(/1p, " MACH1 =", 5E13.6, /(8X, 5E13.6))') (Mach1(i), i=1, Nsk)
  end if
  iof = 0
200 iof = iof + 1
  Oxfl = Oxf(iof)
  call NEWOF
  Incdeq = seql
300 refl = .false.
  it2 = 2
  it1 = 1
  Pp = P(1)
  Tt = T(1)
  if (.not. Incdeq) then
! FROZEN
     do n = 1, Nsk
        Dlvtp(n) = 1.
        Dlvpt(n) = -1.
     end do
  end if
  do Npt = 1, Nsk
     Ppp(Npt) = P(Npt)
     Ttt(Npt) = T(Npt)
     if (Npt > 1) then
        if (Ppp(Npt) == 0.) Ppp(Npt) = Ppp(Npt-1)
        if (Ttt(Npt) == 0.) Ttt(Npt) = Ttt(Npt-1)
        Ssum(Npt) = Ssum(Npt-1)
        Hsum(Npt) = Hsum(Npt-1)
        if (Ttt(Npt) == Tt .and. Ppp(Npt) == Pp) go to 350
     end if
     Pp = Ppp(Npt)
     Tt = Ttt(Npt)
     if (Tt >= Tg(1)*.8d0) then
        call HCALC
        Hsum(Npt) = Hsub0
     else
        write(IOOUT, '(/" TEMPERATURE=", E12.4, " IS OUT OF EXTENDED RANGE ", &
             & "FOR POINT", I5, " (SHCK)")') Tt, Npt
        go to 1000
     end if
350  if (Cpmix /= 0.) Gamma1 = Cpmix/(Cpmix-1./Wmix)
     A1 = (Rr*Gamma1*Tt/Wmix)**.5
     if (U1(Npt) == 0.) U1(Npt) = A1*Mach1(Npt)
     if (Mach1(Npt) == 0.) Mach1(Npt) = U1(Npt)/A1
     Wm(Npt) = Wmix
     Cpr(Npt) = Cpmix
     Gammas(Npt) = Gamma1
     Vlm(Npt) = Rr*Tt/(Wmix*Pp)
  end do
  Npt = Nsk
! OUTPUT--1ST CONDITION
  write(IOOUT, '(////25X, "SHOCK WAVE PARAMETERS ASSUMING")')
  if (.not. Incdeq) then
     write(IOOUT, '(/, 17X, " FROZEN COMPOSITION FOR INCIDENT SHOCKED CONDITI1ONS"//)')
  else
     write(IOOUT, '(/, 16X, " EQUILIBRIUM COMPOSITION FOR INCIDENT SHOCKED CONDITIONS"//)')
  end if
  Eql = .false.
  call OUT1
  write(IOOUT, '(/" INITIAL GAS (1)")')
  Fmt(4) = '13'
  Fmt(5) = ' '
  Fmt(7) = '4,'
  write(IOOUT, Fmt) 'MACH NUMBER1   ', (Mach1(j), j=1, Npt)
  Fmt(7) = '2,'
  write(IOOUT, Fmt) 'U1, M/SEC      ', (U1(j), j=1, Npt)
  call OUT2
! BEGIN CALCULATIONS FOR 2ND CONDITION
  if (Incdeq) Eql = .true.
  Npt = 1
400 Gamma1 = Gammas(Npt)
  uu = U1(Npt)
  wmx = Wm(Npt)
  p1 = Ppp(Npt)
  t1 = Ttt(Npt)
  hs = Hsum(Npt)
  if (refl) uu = u1u2(Npt)
  mu12rt = wmx*uu**2/(Rr*t1)
  if (refl) then
! REFLECTED--SUBSCRIPTS 2=1, 5=2, P52=P21
     t21 = 2.
     b2 = (-1.-mu12rt-t21)/2.
     p21 = -b2 + SQRT(b2**2-t21)
  else
     p21 = (2.*Gamma1*Mach1(Npt)**2-Gamma1+1.)/(Gamma1+1.)
! THE FOLLOWING IMPROVED FORMULATION FOR THE INITIAL ESTIMATE FOR THE
! 2ND CONDITION WAS MADE AND TESTED BY S. GORDON 7/10/89.
     if (.not. Eql) then
        t21 = p21*(2./Mach1(Npt)**2+Gamma1-1.)/(Gamma1+1.)
     else
        Pp = p21*p1
        Tp = .false.
        Hp = .true.
        Hsub0 = hs + uu**2/(2.*Rr)
        call EQLBRM
        t21 = Ttt(Npt)/t1
        Hp = .false.
        Tp = .true.
     end if
  end if
  p21l = log(p21)
  ttmax = 1.05*Tg(4)/t1
  t21 = min(t21, ttmax)
  t21l = log(t21)
  itr = 1
500 if (Shkdbg) write(IOOUT, '(/" ITR NO.=", I3, 3X, "P", I1, "/P", I1, " =", F9.4, 3X, "T", I1, &
       & "/T", I1, " =", F9.4, "   RHO2/RHO1 =", F9.6)') itr, it2, it1, p21, it2, it1, t21, rho52
  Tt = t21*t1
  Pp = p21*p1
  if (.not. Eql) then
! FROZEN
     Tln = log(Tt)
     if (.not. Incdeq) then
        call HCALC
        if (Tt == 0.) go to 600
        Hsum(Npt) = Hsub0
        Cpr(Npt) = Cpmix
     else
        call CPHS
        Cpr(Npt) = Cpsum
        Hsum(Npt) = 0.
        do j = 1, Ng
           Hsum(Npt) = Hsum(Npt) + H0(j)*En(j, Npt)
        end do
        Hsum(Npt) = Hsum(Npt)*Tt
     end if
  else
     call EQLBRM
     if (Tt == 0.) go to 800
  end if
  rho12 = wmx*t21/(Wm(Npt)*p21)
  gg = rho12*mu12rt
  rho52 = 1./rho12
  if (refl) gg = -mu12rt*rho52/(rho52-1.)**2
  G(1, 1) = -gg*Dlvpt(Npt) - p21
  G(1, 2) = -gg*Dlvtp(Npt)
  G(1, 3) = p21 - 1. + gg - mu12rt
  if (refl) G(1, 3) = p21 - 1. + gg*(rho52-1.)
  gg = gg*t1/wmx
  if (.not. refl) gg = gg*rho12
  G(2, 1) = -gg*Dlvpt(Npt) + Tt*(Dlvtp(Npt)-1.)/Wm(Npt)
  G(2, 2) = -gg*Dlvtp(Npt) - Tt*Cpr(Npt)
  gg = 1. - rho12**2
  if (refl) gg = (rho52+1.)/(rho52-1.)
  G(2, 3) = Hsum(Npt) - hs - uu**2*gg/(2.*Rr)
  X(3) = G(1, 1)*G(2, 2) - G(1, 2)*G(2, 1)
  X(1) = (G(1, 3)*G(2, 2)-G(2, 3)*G(1, 2))/X(3)
  X(2) = (G(1, 1)*G(2, 3)-G(2, 1)*G(1, 3))/X(3)
  if (Shkdbg) then
     write(IOOUT, '(/" G(I,J)  ", 3E15.8)') G(1, 1), G(1, 2), G(1, 3)
     write(IOOUT, '(/" G(I,J)  ", 3E15.8)') G(2, 1), G(2, 2), G(2, 3)
     write(IOOUT, '(/" X       ", 2E15.8)') X(1), X(2)
     write(IOOUT, '(/" HSUM HS UU U2 ", 4E15.8)') Hsum(Npt), hs, uu, uu*rho12
  end if
  ax = abs(X(1))
  axx = abs(X(2))
  if (axx > ax) ax = axx
  if (ax >= .00005) then
     cormax = .40546511
     if (itr > 4) cormax = .22314355
     if (itr > 12) cormax = .09531018
     if (itr > 20) cormax = .04879016
     ax = ax/cormax
     if (ax > 1.) then
        X(1) = X(1)/ax
        X(2) = X(2)/ax
     end if
     p21l = p21l + X(1)
     t21l = t21l + X(2)
     p21 = exp(p21l)
     t21 = exp(t21l)
     if (Shkdbg) write(IOOUT, '(/" MAX.COR.=", e13.6, " X(1)=", e13.6, " X(2)=", e13.6)') cormax, X(1), X(2)
     if (itr /= 1 .or. t21 < ttmax) then
        itr = itr + 1
        if (itr < 61) go to 500
        write(IOOUT, '(/6x, " WARNING!!  NO CONVERGENCE FOR u1=", F8.1, &
             & /"  ANSWERS NOT RELIABLE, SOLUTION MAY NOT EXIST (SHCK)")') U1(Npt)
     else
        Tt = 0.
        Npt = Npt - 1
        go to 700
     end if
  end if
! CONVERGED OR TOOK 60 ITERATIONS WITHOUT CONVERGING.
! STORE RESULTS.
600 rrho(Npt) = rho52
  m2m1(Npt) = Wm(Npt)/wmx
  p2p1(Npt) = p21
  t2t1(Npt) = t21
  utwo(Npt) = uu*rho12
  u1u2(Npt) = uu - utwo(Npt)
  if (Tt >= Tg(1)*.8d0 .and. Tt <= Tg(4)*1.1d0) then
     if (.not. Eql) then
! FROZEN
        Ppp(Npt) = Pp
        Ttt(Npt) = Tt
        Gammas(Npt) = Cpr(Npt)/(Cpr(Npt)-1./wmx)
        Vlm(Npt) = Rr*Tt/(wmx*Pp)
        if (Incdeq) then
           Ssum(Npt) = 0.
           do j = 1, Ngc
              pmn = Pp*wmx*En(j, Npt)
              if (En(j, Npt) > 0.) Ssum(Npt) = Ssum(Npt) + En(j, Npt) &
                   *(S(j)-log(pmn))
           end do
        end if
     end if
     go to 900
  end if
700 write(IOOUT, '(/" TEMPERATURE=", E12.4, " IS OUT OF EXTENDED RANGE ", &
       & "FOR POINT", I5, " (SHCK)")') Tt, Npt
  Tt = 0.
800 if (Npt < 1) go to 1000
  Nsk = Npt
900 if (Trnspt) call TRANP
  Isv = 0
  if (Npt < Nsk) Isv = Npt
  if (Npt == 1) Isv = -1
  Npt = Npt + 1
  if (Eql) call SETEN
  if (Npt <= Nsk) go to 400
  Npt = Nsk
  if (refl) then
     if (.not. Eql) write(IOOUT, '(/" SHOCKED GAS (5)--REFLECTED--FROZEN")')
     if (Eql) write(IOOUT, '(/" SHOCKED GAS (5)--REFLECTED--EQUILIBRIUM")')
     cr12 = '2'
     cr52 = '5'
  else
     if (.not. Eql) write(IOOUT, '(/" SHOCKED GAS (2)--INCIDENT--FROZEN")')
     if (Eql) write(IOOUT, '(/" SHOCKED GAS (2)--INCIDENT--EQUILIBRIUM")')
     cr12 = '1'
     cr52 = '2'
  end if
  Fmt(7) = '2,'
  write(IOOUT, Fmt) 'U'//cr52//', M/SEC      ', (utwo(j), j=1, Npt)
  call OUT2
  if (Trnspt) call OUT4
  write(IOOUT, *)
  Fmt(7) = '3,'
  write(IOOUT, Fmt) 'P'//cr52//'/P'//cr12//'           ', &
       (p2p1(j), j=1, Npt)
  write(IOOUT, Fmt) 'T'//cr52//'/T'//cr12//'           ', &
       (t2t1(j), j=1, Npt)
  Fmt(7) = '4,'
  write(IOOUT, Fmt) 'M'//cr52//'/M'//cr12//'           ', &
       (m2m1(j), j=1, Npt)
  write(IOOUT, Fmt) 'RHO'//cr52//'/RHO'//cr12//'       ', &
       (rrho(j), j=1, Npt)
  Fmt(7) = '2,'
  if (.not. refl) write(IOOUT, Fmt) 'V2, M/SEC      ', (u1u2(j), &
       j=1, Npt)
  if (refl) write(IOOUT, Fmt) 'U5+V2,M/SEC    ', (u1u2(j), j=1, Npt)
  if (.not. Eql) then
! WRITE FROZEN MOLE (OR MASS) FRACTIONS
     Fmt(7) = '5,'
     if (.not. Incdeq) then
        if (Massf) then
           write(IOOUT, '(/1x, A4, " FRACTIONS"/)') 'MASS'
        else
           write(IOOUT, '(/1x, A4, " FRACTIONS"/)') 'MOLE'
           ww = wmx
        end if
        do n = 1, Nreac
           j = Jray(n)
           if (Massf) ww = Mw(j)
           write(IOOUT, '(" ", A16, F8.5, 12F9.5)') Prod(j), (En(j, i)*ww, i=1, Npt)
        end do
     else
        Eql = .true.
        call OUT3
        Eql = .false.
     end if
  else
     call OUT3
  end if
  Iplt = min(Iplt+Npt, 500)
  if (srefl) then
     if (.not. refl) then
        refl = .true.
        it2 = 5
        it1 = 2
        Eql = .true.
        if (Reflfz) then
           Eql = .false.
           if (Refleq) then
              j = 0
              do i = 1, Npt
                 j = j + 1
                 sg(j) = u1u2(i)
                 j = j + 1
                 sg(j) = Wm(i)
                 j = j + 1
                 sg(j) = Ppp(i)
                 j = j + 1
                 sg(j) = Ttt(i)
                 j = j + 1
                 sg(j) = Hsum(i)
                 j = j + 1
                 sg(j) = Gammas(i)
              end do
           end if
        end if
        Npt = 1
        go to 400
     else if (.not. Eql .and. Refleq) then
        j = 1
        do i = 1, Npt
           u1u2(i) = sg(j)
           Wm(i) = sg(j+1)
           Ppp(i) = sg(j+2)
           Ttt(i) = sg(j+3)
           Hsum(i) = sg(j+4)
           Gammas(i) = sg(j+5)
           j = j + 6
        end do
        Eql = .true.
        Npt = 1
        go to 400
     end if
  end if
  if (Incdeq .and. Incdfz) then
     Incdeq = .false.
     Eql = .false.
     go to 300
  else if (iof >= Nof) then
     Tp = .false.
     do n = 1, Nreac
        Rtemp(n) = T(1)
     end do
  else
     do i = 1, Nsk
        U1(i) = uis(i)
        Mach1(i) = mis(i)
     end do
     go to 200
  end if
1000 return
end subroutine SHCK



subroutine THERMP
!***********************************************************************
! ASSIGNED THERMODYNAMIC STATES.  HP, SP, TP, UV, SV, AND TV PROBLEMS.
!***********************************************************************
  use cea
  implicit none
! LOCAL VARIABLES
  integer, save:: iof
  logical:: Uv, Tv, Sv

  Uv = transfer(Hp, Uv)
  Tv = transfer(Tp, Tv)
  Sv = transfer(Sp, Sv)
  Eql = .true.
  do 100 iof = 1, Nof
     Oxfl = Oxf(iof)
     call NEWOF
! SET ASSIGNED P OR VOLUME
     do Ip = 1, Np
        Pp = P(Ip)
! SET ASSIGNED T
        do It = 1, Nt
           Vv = V(Ip)
           Tt = T(It)
           call EQLBRM
           if (Npt == 0) go to 200
           if (Trnspt .and. Tt /= 0.) call TRANP
           Isv = 0
           if (Ip /= Np .or. It /= Nt .and. Tt /= 0.) then
              Isv = Npt
              if (Npt /= Ncol) go to 10
           end if
           if (.not. Hp) write(IOOUT, '(////15X, "THERMODYNAMIC EQUILIBRIUM PROPERTIES AT ASSIGNED")')
           if (Hp) write(IOOUT, '(////9X, "THERMODYNAMIC EQUILIBRIUM COMBUSTION PROPERTIES AT ASSIGNED")')
           if (.not. Vol) then
              if (Hp) write(IOOUT, '(/34X, " PRESSURES"/)')
              if (Tp) write(IOOUT, '(/27X, "TEMPERATURE AND PRESSURE"/)')
              if (Sp) write(IOOUT, '(/29X, "ENTROPY AND PRESSURE"/)')
           else
              if (Uv) write(IOOUT, '(/36X, " VOLUME"/)')
              if (Tv) write(IOOUT, '(/28X, "TEMPERATURE AND VOLUME"/)')
              if (Sv) write(IOOUT, '(/30X, "ENTROPY AND VOLUME"/)')
           end if
           call OUT1
           write(IOOUT, '(/" THERMODYNAMIC PROPERTIES"/)')
           call OUT2
           if (Trnspt) call OUT4
           call OUT3
           Iplt = min(Iplt+Npt, 500)
           if (Isv == 0 .and. iof == Nof) go to 200
           write(IOOUT, '(////)')
           Npt = 0
10         Npt = Npt + 1
           if (.not. Tp .and. Tt /= 0.) T(1) = Tt
           if (Nt == 1 .and. Np == 1) go to 100
           if (Ip == 1 .and. It == 1) Isv = -Isv
           if (Nt /= 1) then
              if (It == Nt .or. Tt == 0.) Isv = 0
           end if
           call SETEN
        end do
     end do
100 continue
200 return
end subroutine



subroutine TRANIN
!***********************************************************************
! BRINGS IN AND SORTS OUT INPUT FOR TRANSPORT CALCULATIONS
!***********************************************************************
  use cea
  implicit none
! LOCAL VARIABLES
  integer, save:: i, ii, inds(maxTr), ir, j, jtape(2), k, k1, k2, kt, kvc, l, loop, m, nms
  logical, save:: change, elc1, elc2, ion1, ion2, setx
  real(8), save:: coeff, debye, ekt, enel, enmin, ionic, lamda, omega, prop, qc, ratio, &
       stcf(maxTr, maxTr), stcoef(maxTr), te, testen, testot, total, &
       trc(6, 3, 2), wmols(maxTr), wmred, xsel, xss(maxTr)

  if (.not. Eql) then
     if (.not. Shock) then
        if (.not. setx) then
           setx = .true.
           Nm = nms
           do i = 1, Nm
              Xs(i) = xss(i)
              Wmol(i) = wmols(i)
              Ind(i) = inds(i)
           end do
        end if
        go to 300
     else if (.not. Incdeq) then
        if (Npt <= 1) then
           Nm = Nreac
           do i = 1, Nm
              j = Jray(i)
              Ind(i) = j
              Wmol(i) = Mw(j)
              Xs(i) = En(j, 1)*Wm(1)
           end do
        end if
        go to 300
     end if
  end if
! PICK OUT IMPORTANT SPECIES
  Nm = 0
  total = 0.d0
  enmin = 1.0d-11/Wm(Npt)
  testot = 0.999999999d0/Wm(Npt)
  do i = 1, Lsave
     j = Jcm(i)
     if (En(j, Npt) <= 0.d0 .and. j <= Ngc) then
        if ((Enln(j)-Ennl+25.328436d0) > 0.d0) En(j, Npt) &
             = exp(Enln(j))
     end if
     Nm = Nm + 1
     Ind(Nm) = j
     total = total + En(j, Npt)
     if (Mw(j) < 1.0d0) enel = En(j, Npt)
     En(j, Npt) = -En(j, Npt)
  end do
  testen = 1.d0/(Ng*Wm(Npt))
  loop = 0
100 if (total <= testot .and. loop <= Ng) then
     loop = loop + 1
     testen = testen/10.
     do j = 1, Ng
        if (En(j, Npt) >= testen) then
           if (Nm >= maxTr) then
              write(IOOUT, '(/" WARNING!!  MAXIMUM ALLOWED NO. OF SPECIES", I3, " WAS USED IN ", &
                   & /" TRANSPORT PROPERTY CALCULATIONS FOR POINT", I3, "(TRANIN))")') Nm, Npt
              go to 200
           else
              total = total + En(j, Npt)
              Nm = Nm + 1
              Ind(Nm) = j
              En(j, Npt) = -En(j, Npt)
           end if
        end if
     end do
     if (testen > enmin) go to 100
  end if
! CALCULATE MOLE FRACTIONS FROM THE EN(J, NPT)
200 do j = 1, Ng
     En(j, Npt) = abs(En(j, Npt))
  end do
  do i = 1, Nm
     j = Ind(i)
     Wmol(i) = Mw(j)
     Xs(i) = En(j, Npt)/total
  end do
  if (Npt == Nfz) then
     nms = Nm
     do i = 1, Nm
        xss(i) = Xs(i)
        wmols(i) = Wmol(i)
        inds(i) = Ind(i)
     end do
     setx = .false.
  end if
! REWRITE REACTIONS TO ELIMINATE TRACE ELEMENTS
  Nr = Nm - Lsave
  if (Nr /= 0) then
     do k = 1, maxTr
        do m = 1, maxTr
           Stc(k, m) = 0.0d0
        end do
     end do
     k = 1
     do i = Lsave + 1, Nm
        Stc(k, i) = -1.0d0
        j = Ind(i)
        do m = 1, Lsave
           Stc(k, m) = A(m, j)
        end do
        k = k + 1
     end do
     do i = 1, Nm
        if (Xs(i) < 1.0d-10) then
           m = 1
           change = .false.
           do 210 j = 1, Nr
              coeff = Stc(j, i)
              if (ABS(coeff) > 1.0d-05) then
                 if (.not. change) then
                    change = .true.
                    do k = 1, Nm
                       stcoef(k) = Stc(j, k)/coeff
                    end do
                    go to 210
                 else
                    do k = 1, Nm
                       Stc(j, k) = (Stc(j, k)/coeff) - stcoef(k)
                    end do
                 end if
              end if
              do k = 1, Nm
                 stcf(m, k) = Stc(j, k)
              end do
              m = m + 1
210        continue
           do ii = 1, Nm
              do j = 1, Nr
                 Stc(j, ii) = stcf(j, ii)
              end do
           end do
           Nr = m - 1
        end if
     end do
  end if
! FIND TRANSPORT DATA FOR IMPORTANT INTERACTIONS
300 do i = 1, Nm
     Con(i) = 0.0
     do j = 1, Nm
        Eta(i, j) = 0.0
     end do
  end do
  rewind IOSCH
  do 400 ir = 1, Ntape
     read(IOSCH) jtape, trc
     do 350 k = 1, 2
        do i = 1, Nm
           j = Ind(i)
           if (j == jtape(k)) then
              l = i
              if (k == 2) then
                 kvc = 1
302              kt = 1
                 if (trc(2, 1, kvc) /= 0.E0) then
                    if (trc(2, 2, kvc) /= 0.E0) then
                       if (Tt > trc(2, 1, kvc)) kt = 2
                       if (trc(2, 3, kvc) /= 0.) then
                          if (Tt > trc(2, 2, kvc)) kt = 3
                       end if
                    end if
                    prop = EXP(trc(6, kt, kvc) &
                         +(trc(5, kt, kvc)/Tt+trc(4, kt, kvc)) &
                         /Tt+trc(3, kt, kvc)*Tln)
                    if (kvc == 2) then
                       Con(l) = prop
                       go to 400
                    else
                       Eta(l, m) = prop
                       if (l /= m) Eta(m, l) = Eta(l, m)
                    end if
                 else if (kvc == 2) then
                    go to 400
                 end if
                 kvc = 2
                 go to 302
              else
                 m = i
                 go to 350
              end if
           end if
        end do
        go to 400
350  continue
400 continue
! MAKE ESTIMATES FOR MISSING DATA
!
! INCLUDES ION CROSS SECTION ESTIMATES
! ESTIMATES FOR  E-ION, ION-ION, E-NEUTRAL, ION-NEUTRAL
! DEBYE SHIELDING WITH IONIC CUTOFF DISTANCE
  if (Ions) then
     te = Tt/1000.d0
     ekt = 4.8032d0**2/(Boltz*te)
     qc = 100.d0*(ekt**2)
     xsel = enel/total
     if (xsel < 1.0d-12) xsel = 1.0d-12
     debye = ((22.5d0/Pi)*(Rr/Avgdr*100.d0)*(te/xsel))/ekt**3
     ionic = ((810.d0/(4.0d0*Pi))*(Rr/Avgdr*100d0)*(te/xsel)) &
          **(2.0/3.0)/ekt**2
     lamda = sqrt(debye+ionic)
     lamda = max(lamda, 2.71828183d0)
  end if
  do i = 1, Nm
     k = Ind(i)
     Cprr(i) = Cp(k)
     if (.not. (Ions .and. (abs(A(Nlm, k)) == 1.d0) .and.  &
          (Eta(i, i) == 0.d0))) then
        if (Eta(i, i) == 0.d0) then
           omega = log(50.d0*Wmol(i)**4.6/Tt**1.4)
           omega = max(omega, 1.d0)
           Eta(i, i) = Viscns*sqrt(Wmol(i)*Tt)/omega
        end if
        if (Con(i) == 0.d0) Con(i) = Eta(i, i) &
             *Rr*(.00375d0+.00132d0*(Cprr(i)-2.5d0))/Wmol(i)
     end if
  end do
  do i = 1, Nm
     do 450 j = i, Nm
        ion1 = .false.
        ion2 = .false.
        elc1 = .false.
        elc2 = .false.
        omega = 0.0
        if (Eta(i, j) == 0.) Eta(i, j) = Eta(j, i)
        if (Eta(j, i) == 0.) Eta(j, i) = Eta(i, j)
        if (Eta(i, j) == 0.) then
           if (Ions) then
! ESTIMATE FOR IONS
              k1 = Ind(i)
              k2 = Ind(j)
              if (ABS(A(Nlm, k1)) == 1.0) ion1 = .true.
              if (ABS(A(Nlm, k2)) == 1.0) ion2 = .true.
              if (Wmol(i) < 1.0) elc1 = .true.
              if (Wmol(j) < 1.0) elc2 = .true.
              if (ion1 .and. ion2) omega = 1.36d0*qc*log(lamda)
              if ((ion1 .and. elc2) .or. (ion2 .and. elc1)) &
                   omega = 1.29d0*qc*log(lamda)
              if ((ion1 .and. .not. ion2) .or. (ion2 .and. .not. ion1)) &
                   omega = EXP(6.776-0.4*Tln)
              if (omega /= 0.) then
                 wmred = sqrt(2.0*Tt*Wmol(i)*Wmol(j)/(Wmol(i)+Wmol(j)))
                 Eta(i, j) = Viscns*wmred*Pi/omega
                 Eta(j, i) = Eta(i, j)
                 if (i == j) then
                    Cprr(i) = Cp(k1)
                    Con(i) = Eta(i, i) &
                         *Rr*(.00375d0+.00132d0*(Cprr(i)-2.5d0)) &
                         /Wmol(i)
                 end if
                 go to 450
              end if
           end if
! ESTIMATE FOR UNLIKE INTERACTIONS FROM RIGID SPHERE ANALOGY
           ratio = sqrt(Wmol(j)/Wmol(i))
           Eta(i, j) = 5.656854d0*Eta(i, i) &
                *SQRT(Wmol(j)/(Wmol(i)+Wmol(j)))
           Eta(i, j) = Eta(i, j)/(1.d0+sqrt(ratio*Eta(i, i)/Eta(j, j)))**2
           Eta(j, i) = Eta(i, j)
        end if
450  continue
  end do
  return
end subroutine



subroutine TRANP
!***********************************************************************
! CALCULATES GAS TRANSPORT PROPERTIES
!
!   NUMBER OF GASEOUS SPECIES = NM   (MAXIMUM maxTr)
!   NUMBER OF CHEMICAL REACTIONS = NR (NM - NLM)
!   ARRAY OF STOICHIOMETRIC COEFFICIENTS = STC
!***********************************************************************
  use cea
  implicit none
! LOCAL VARIABLES
  integer, save:: i, i1, j, jj, k, m, mm, nlmm, nmm
  real(8), save:: cpreac, delh(maxTr), gmat(maxMat, maxMat+1), phi(maxTr, maxTr), &
       psi(maxTr, maxTr), reacon, rtpd(maxTr, maxTr), stx(maxTr), &
       stxij(maxTr, maxTr), sumc, sumv, wtmol, xskm(maxTr, maxTr)

  call TRANIN
! CALCULATE VISCOSITY AND FROZEN THERMAL CONDUCTIVITY
  nmm = Nm - 1
  do i = 1, Nm
     rtpd(i, i) = 0.d0
     phi(i, i) = 1.d0
     psi(i, i) = 1.d0
  end do
  Confro(Npt) = 0.d0
  Vis(Npt) = 0.d0
  do i = 1, nmm
     i1 = i + 1
!DIR$ IVDEP
     do j = i1, Nm
        sumc = 2.d0/(Eta(i, j)*(Wmol(i)+Wmol(j)))
        phi(i, j) = sumc*Wmol(j)*Eta(i, i)
        phi(j, i) = sumc*Wmol(i)*Eta(j, j)
        sumc = (Wmol(i)+Wmol(j))**2
        psi(i, j) = phi(i, j) &
             *(1.d0+2.41d0*(Wmol(i)-Wmol(j))*(Wmol(i)-.142d0* &
             Wmol(j))/sumc)
        psi(j, i) = phi(j, i) &
             *(1.d0+2.41d0*(Wmol(j)-Wmol(i))*(Wmol(j)-.142d0* &
             Wmol(i))/sumc)
     end do
  end do
  do i = 1, Nm
     sumc = 0.d0
     sumv = 0.d0
     do j = 1, Nm
        sumc = sumc + psi(i, j)*Xs(j)
        sumv = sumv + phi(i, j)*Xs(j)
     end do
     Vis(Npt) = Vis(Npt) + Eta(i, i)*Xs(i)/sumv
     Confro(Npt) = Confro(Npt) + Con(i)*Xs(i)/sumc
  end do
  if (Eql .and. Nr > 0) then
! CALCULATE REACTION HEAT CAPACITY AND THERMAL CONDUCTIVITY
     m = Nr + 1
     do i = 1, Nr
        delh(i) = 0.0d0
        do k = 1, Lsave
           j = Jcm(k)
           delh(i) = Stc(i, k)*H0(j) + delh(i)
        end do
        nlmm = Lsave + 1
        do k = nlmm, Nm
           j = Ind(k)
           delh(i) = Stc(i, k)*H0(j) + delh(i)
        end do
        G(i, m) = delh(i)
     end do
     do i = 1, maxTr
        do j = 1, maxTr
           if (abs(Stc(i, j)) < 1.0d-6) Stc(i, j) = 0.0d0
        end do
     end do
     jj = Nm - 1
     do k = 1, jj
        mm = k + 1
        do m = mm, Nm
           rtpd(k, m) = Wmol(k)*Wmol(m)/(1.1*Eta(k, m)*(Wmol(k)+Wmol(m)))
           xskm(k, m) = Xs(k)*Xs(m)
           xskm(m, k) = xskm(k, m)
           rtpd(m, k) = rtpd(k, m)
        end do
     end do
     do i = 1, Nr
        do j = i, Nr
           G(i, j) = 0.0d0
           gmat(i, j) = 0.0d0
        end do
     end do
     do k = 1, jj
        mm = k + 1
        do m = mm, Nm
           if (Xs(k) >= 1.0d-10 .and. Xs(m) >= 1.0d-10) then
              do j = 1, Nr
                 if ((Stc(j, k) == 0.d0) .and. (Stc(j, m) == 0.d0)) stx(j) &
                      = 0.d0
                 if ((Stc(j, k) /= 0.d0) .or. (Stc(j, m) /= 0.d0)) stx(j) &
                      = Xs(m)*Stc(j, k) - Xs(k)*Stc(j, m)
              end do
              do i = 1, Nr
                 do j = i, Nr
                    stxij(i, j) = stx(i)*stx(j)/xskm(k, m)
                    G(i, j) = G(i, j) + stxij(i, j)
                    gmat(i, j) = gmat(i, j) + stxij(i, j)*rtpd(k, m)
                 end do
              end do
           end if
        end do
     end do
     m = 1 + Nr
     do i = 1, Nr
!DIR$ IVDEP
        do j = i, Nr
           G(j, i) = G(i, j)
        end do
        G(i, m) = delh(i)
     end do
     Imat = Nr
     call GAUSS
     cpreac = 0.d0
     do i = 1, Nr
        G(i, m) = delh(i)
        cpreac = cpreac + R*delh(i)*X(i)
!DIR$ IVDEP
        do j = i, Nr
           G(i, j) = gmat(i, j)
           G(j, i) = G(i, j)
        end do
     end do
     call GAUSS
     reacon = 0.d0
     do i = 1, Nr
        reacon = reacon + R*delh(i)*X(i)
     end do
     reacon = .6d0*reacon
  else
     cpreac = 0.d0
     reacon = 0.d0
  end if
! CALCULATE OTHER ANSWERS
  Cpfro(Npt) = 0.d0
  wtmol = 0.d0
  do i = 1, Nm
     Cpfro(Npt) = Cpfro(Npt) + Xs(i)*Cprr(i)
     wtmol = wtmol + Xs(i)*Wmol(i)
  end do
  Cpfro(Npt) = Cpfro(Npt)*R/wtmol
  Confro(Npt) = Confro(Npt)/1000.d0
  if (.not. SIunit) Confro(Npt) = Confro(Npt)/4.184d0
  Vis(Npt) = Vis(Npt)/1000.d0
  Prfro(Npt) = Vis(Npt)*Cpfro(Npt)/Confro(Npt)
  if (Eql) then
     cpreac = cpreac/wtmol
     reacon = reacon/1000.d0
     Cpeql(Npt) = cpreac + Cpfro(Npt)
     Coneql(Npt) = Confro(Npt) + reacon
     Preql(Npt) = Vis(Npt)*Cpeql(Npt)/Coneql(Npt)
  end if
end subroutine TRANP



subroutine UTHERM(readOK)
!***********************************************************************
! READ THERMO DATA FROM I/O UNIT 7 IN RECORD format AND WRITE
! UNformatTED ON I/O UNIT IOTHM.  DATA ARE REORDERED GASES FIRST.
!
! USES SCRATCH I/O UNIT IOSCH.
!
! UTHERM IS CALLED FROM subroutine INPUT.
!
! NOTE:  THIS ROUTINE MAY BE CALLED DIRECTLY AND USED BY ITSELF TO
! PROCESS THE THERMO DATA.
!
! GASEOUS SPECIES:
! THE STANDARD TEMPERATURE RANGES TGL ARE GIVEN ON THE FIRST
! RECORD, FOLLOWED BY THE DATE OF THE LAST DATA CHANGE THDATE.
!
! WHEN COEFFICIENTS ARE NOT GIVEN FOR THE THIRD TEMPERATURE
! INTERVAL, A STRAIGHT LINE FOR CP/R IS USED.  FOR HIGH TEMPS, 
! THE EXTRAPOLATION GOES BETWEEN THE LAST POINT GIVEN AND THE
! FOLLOWING VALUES AA AT TINF=1.d06 K:
!      MONATOMICS  2.5
!      DIATOMICS   4.5
!      POLYATOMICS 3*N-1.75  (AVERAGE 1.5 AND 2)
!
! THE FOLLOWING EXTRAPOLATION IS NOT CURRENTLY PROGRAMED (12/9/98):
!   FOR LOW TEMPS, THE EXTRAPOLATION GOES BETWEEN THE FIRST VALUE
!   DOWN TO THE FOLLOWING VALUES AA AT 0 K:
!      MONATOMICS  2.5
!      DIATOMICS   3.5
!      POLYATOMICS 3.75 (AVERAGE 3.5 AND 4.0)
!
! IF DATA ARE AVAILABLE FOR THE THIRD T INTERVAL, IFAZ (SEE
! DEFINITION) IS SET TO -1 AND THE NAME IS ALTERED TO START WITH *.
!
! CONDENSED SPECIES:
! NO EXTRAPOLATIONS ARE DONE.  TEMP INTERVALS VARY.
!
! SOME DEFINITIONS:
! TGL(I)  - TEMPERATURE INTERVALS FOR GASES (I.E. 200, 1000, 6000, 20000).
! FILL(I) - IF TRUE, DATA MISSING FOR INTERVAL.  CURRENTLY ONLY 3RD
!           INTERVAL CHECKED.
! NGL     - NUMBER OF GASEOUS PRODUCTS.
! NS      - NGL + NUMBER OF CONDENSED PRODUCT PHASES.
! NALL    - NS + NUMBER OF REACTANT SPECIES.
! IFAZ    - PHASE INDICATOR. GASES ARE 0, CONDENSED PHASES ARE NUMBERED
!           STARTING WITH 1 FOR THE LOWEST T RANGE, 2 FOR THE NEXT
!           CONTIGUOUS PHASE, ETC.
! NTL     - NUMBER OF T INTERVALS FOR A SPECIES SET.
!***********************************************************************
  use cea
  implicit none
! DUMMY ARGUMENTS
  logical:: readOK
! LOCAL VARIABLES
  character(15):: name
  character(16):: namee
  character(65):: notes
  character(2):: sym(5)
  character(6):: date
  integer:: i, ifaz, ifzm1, inew, int, j, k, kk, l, nall, ncoef, ngl, ns, ntl
  logical:: fill(3)
  real(8):: aa, atms, cpfix, dlt, expn(8), fno(5), hform, hh, mwt, templ(9), tex, &
       tgl(4), thermo(9, 3), tinf, tl(2), ttl, tx

  ngl = 0
  ns = 0
  nall = 0
  ifzm1 = 0
  inew = 0
  tinf = 1.d06
  rewind IOSCH
  read(IOINP, '(4F10.3, a10)') tgl, Thdate
100 do i = 1, 3
     fill(i) = .true.
     do j = 1, 9
        thermo(j, i) = 0.
     end do
  end do
  hform = 0.
  tl(1) = 0.
  tl(2) = 0.
  read(IOINP, '(a15, a65)', END=300, ERR=400) name, notes
  if (name(:3) == 'END' .or. name(:3) == 'end') then
     if (index(name, 'ROD') == 0 .and. index(name, 'rod') == 0) &
          go to 300
     ns = nall
     go to 100
  end if
  read(IOINP, '(i2, 1x, a6, 1x, 5(a2, f6.2), i2, f13.5, f15.3)', ERR=400) ntl, date, (sym(j), fno(j), j=1, 5), &
       ifaz, mwt, hform
  write(IOOUT, '(" ", a15, 2x, a6, e15.6, 2x, a65)') name, date, hform, notes
! IF NTL=0, REACTANT WITHOUT COEFFICIENTS
  if (ntl == 0) then
     if (ns == 0) go to 300
     nall = nall + 1
     read(IOINP, '(2F11.3, i1, 8F5.1, 2x, f15.3)', ERR=400) tl, ncoef, expn, hh
     thermo(1, 1) = hform
     write(IOSCH) name, ntl, date, (sym(j), fno(j), j=1, 5), ifaz, tl, mwt, &
          thermo
     go to 100
  else if (name == 'Air') then
     sym(1) = 'N'
     fno(1) = 1.56168d0
     sym(2) = 'O'
     fno(2) = .419590d0
     sym(3) = 'AR'
     fno(3) = .009365d0
     sym(4) = 'C'
     fno(4) = .000319d0
  else if (name == 'e-') then
     mwt = 5.48579903d-04
  end if
  do 200 i = 1, ntl
     read(IOINP, '(2F11.3, i1, 8F5.1, 2x, f15.3)', ERR=400) tl, ncoef, expn, hh
     read(IOINP, '(5d16.8/2d16.8, 16x, 2d16.8)', ERR=400) templ
     if (ifaz == 0 .and. i > 3) go to 400
     if (ifaz <= 0) then
        if (tl(2) > tgl(4)-.01d0) then
           ifaz = -1
           namee = '*'//name
           name = namee(:15)
        end if
        if (tl(1) >= tgl(i+1)) go to 200
        int = i
        fill(i) = .false.
     else
        int = 1
        if (i > 1) then
           do k = 1, 7
              thermo(k, 1) = 0.d0
           end do
        end if
     end if
     do 150 l = 1, ncoef
        do k = 1, 7
           if (expn(l) == real(k-3)) then
              thermo(k, int) = templ(l)
              go to 150
           end if
        end do
150  continue
     thermo(8, int) = templ(8)
     thermo(9, int) = templ(9)
     if (ifaz > 0) then
        nall = nall + 1
        if (ifaz > ifzm1) then
           inew = inew + 1
        else
           inew = i
        end if
        write(IOSCH) name, ntl, date, (sym(j), fno(j), j=1, 5), inew, tl, mwt, &
             thermo
     end if
200 continue
  ifzm1 = ifaz
  if (ifaz <= 0) then
     inew = 0
     nall = nall + 1
     if (ifaz <= 0 .and. ns == 0) then
        ngl = ngl + 1
        if (fill(3)) then
           atms = 0.
           do i = 1, 5
              if (sym(i) == ' ' .or. sym(i) == 'E') go to 210
              atms = atms + fno(i)
           end do
! FOR GASES WITH NO COEFFICIENTS FOR TGL(3)-TGL(4) INTERVAL, 
! CALCULATE ESTIMATED COEFFICIENTS. (STRAIGHT LINE FOR CP/R)
210        aa = 2.5d0
           if (atms > 1.9) aa = 4.5d0
           if (atms > 2.1) aa = 3.*atms - 1.75d0
           ttl = tl(2)
           tx = ttl - tinf
           cpfix = 0
           templ(8) = 0.
           templ(9) = 0.
           dlt = log(ttl)
           do k = 7, 1, - 1
              kk = k - 3
              if (kk == 0) then
                 cpfix = cpfix + thermo(k, 2)
                 templ(8) = templ(8) + thermo(k, 2)
                 templ(9) = templ(9) + thermo(k, 2)*dlt
              else
                 tex = ttl**kk
                 cpfix = cpfix + thermo(k, 2)*tex
                 templ(9) = templ(9) + thermo(k, 2)*tex/kk
                 if (kk == -1) then
                    templ(8) = templ(8) + thermo(k, 2)*dlt/ttl
                 else
                    templ(8) = templ(8) + thermo(k, 2)*tex/(kk+1)
                 end if
              end if
           end do
           templ(2) = (cpfix-aa)/tx
           thermo(4, 3) = templ(2)
           templ(1) = cpfix - ttl*templ(2)
           thermo(3, 3) = templ(1)
           thermo(8, 3) = thermo(8, 2) &
                + ttl*(templ(8)-templ(1)-.5*templ(2)*ttl)
           thermo(9, 3) = -templ(1)*dlt + thermo(9, 2) + templ(9) &
                - templ(2)*ttl
        end if
     end if
! WRITE COEFFICIENTS ON SCRATCH I/O UNIT IOSCH
     write(IOSCH) name, ntl, date, (sym(j), fno(j), j=1, 5), ifaz, tl, mwt, &
          thermo
  end if
  go to 100
! END OF DATA. COPY CONDENSED & REACTANT DATA FROM IOSCH & ADD TO IOTHM.
300 rewind IOSCH
  if (ns == 0) ns = nall
  write(IOTHM) tgl, ngl, ns, nall, Thdate
! WRITE GASEOUS PRODUCTS ON IOTHM
  if (ngl /= 0) then
     do i = 1, ns
        read(IOSCH) name, ntl, date, (sym(j), fno(j), j=1, 5), ifaz, tl, mwt, &
             thermo
        if (ifaz <= 0) write(IOTHM) name, ntl, date, &
             (sym(j), fno(j), j=1, 5), ifaz, tl, mwt, &
             thermo
     end do
  end if
  if (ngl /= nall) then
! WRITE CONDENSED PRODUCTS AND REACTANTS ON IOTHM
     rewind IOSCH
     do i = 1, nall
        read(IOSCH) name, ntl, date, (sym(j), fno(j), j=1, 5), ifaz, tl, mwt, &
             thermo
        if (i > ns) then
           write(IOTHM) name, ntl, date, (sym(j), fno(j), j=1, 5), ifaz, tl, &
                mwt, thermo(1, 1)
           if (ntl > 0) write(IOTHM) thermo
        else if (ifaz > 0) then
           write(IOTHM) name, ntl, date, (sym(j), fno(j), j=1, 5), ifaz, tl, &
                mwt, (thermo(k, 1), k=1, 9)
        end if
     end do
  end if
  return
400 write(IOOUT, '(/" ERROR IN PROCESSING thermo.inp AT OR NEAR ", A15, " (UTHERM)")') name
  readOK = .false.
  return
end subroutine



subroutine UTRAN(readOK)
!***********************************************************************
! READ TRANSPORT PROPERTIES FORM I/O UNIT 7 IN RECORD format AND WRITE
! UNformatTED ON I/O UNIT IOTRN.  USES SCRATCH I/O UNIT IOSCH.
!
! UTRAN IS CALLED FROM subroutine INPUT AFTER A RECORD WITH 'tran'
! IN COLUMNS 1-4 HAS BEEN READ.
!
! NOTE:  THIS ROUTINE MAY BE CALLED DIRECTLY  AND USED BY ITSELF TO
! PROCESS THE TRANSPORT PROPERTY DATA.
!***********************************************************************
  use cea
  implicit none
! DUMMY ARGUMENTS
  logical, intent(out):: readOK
! LOCAL VARIABLES
  character(16):: tname(2)
  character(1):: vorc
  integer:: i, ic, in, iv, j, k, ncc, nn, ns, nv
  real(8):: cc, tcin(6), trcoef(6, 3, 2), vvl


  ns = 0
  rewind IOSCH
  outerLoop: do
     trcoef(:, :, :) = 0.
     read(IOINP, '(2A16, 2X, A1, I1, A1, I1)') tname, vvl, nv, cc, ncc
     if (tname(1) == 'end' .or. tname(1) == 'LAST') then
        write(IOTRN) ns
        rewind IOSCH
        do i = 1, ns
           read(IOSCH, ERR=200) tname, trcoef
           write(IOTRN) tname, trcoef
        end do
        return
     else
        ic = 0
        iv = 0
        nn = nv + ncc
        if (nv <= 3 .and. ncc <= 3) then
           do in = 1, nn
              read(IOINP, '(1X, A1, 2F9.2, 4E15.8)') vorc, tcin
              if (vorc == 'C') then
                 k = 2
                 ic = ic + 1
                 j = ic
              else
                 k = 1
                 iv = iv + 1
                 j = iv
              end if
              if (j > 3) exit outerLoop
              do i = 1, 6
                 trcoef(i, j, k) = tcin(i)
              end do
           end do
           ns = ns + 1
           write(IOSCH) tname, trcoef
           cycle
        end if
     end if
     exit
  end do outerLoop
200 write(IOOUT, '(/" ERROR IN PROCESSING trans.inp AT OR NEAR (UTRAN)", /1X, 2A16)') tname
  readOK = .false.
  return
end subroutine UTRAN



subroutine VARFMT(Vx)
!***********************************************************************
! SET DECIMAL PLACES ACCORDING TO NUMBER SIZE FOR F-format IN
! VARIABLE format FMT.
!***********************************************************************
  use cea
  implicit none
! DUMMY ARGUMENTS
  real(8), intent(in):: Vx(Ncol)
! LOCAL VARIABLES
  integer, save:: i, k
  real(8), save:: vi

  do i = 1, Npt
     vi = abs(Vx(i))
     k = 2*i + 3
     Fmt(k) = '5,'
     if (vi >= 1.) Fmt(k) = '4,'
     if (vi >= 10.) Fmt(k) = '3,'
     if (vi >= 100.) Fmt(k) = '2,'
     if (vi >= 10000.) Fmt(k) = '1,'
     if (vi >= 1000000.) Fmt(k) = '0,'
  end do
  Fmt(29)(2:) = ' '
end subroutine VARFMT
