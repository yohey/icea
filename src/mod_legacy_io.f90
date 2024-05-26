module mod_legacy_io
  use mod_cea
  implicit none

  character(15), private:: fc
  integer, private:: ione, mcond, mcondf, mpn, mpnf, mvis
  real(8), private:: pfuel, phi, tem

contains

  subroutine INPUT(readOK, caseOK, Ensert)
    !***********************************************************************
    ! DECIPHER KEYWORDS, LITERAL VARIABLES, & NUMERICAL VARIABLES IN INPUT.
    !***********************************************************************
    use mod_legacy_cea
    implicit none

    ! DUMMY ARGUMENTS
    logical:: caseOK, readOK
    character(15):: Ensert(20)

    ! LOCAL VARIABLES
    character(26), parameter:: lc = 'abcdefghijklmnopqrstuvwxyz', uc = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(15):: cin(maxNgc), cx15
    character(4):: code, cx4
    character(1):: cx1
    character(2):: cx2
    character(3):: cx3
    logical:: eqrats, incd, phi, pltdat, reacts, refl
    integer:: i, ifrmla, ii, in, iv, ix, j, jj, k, lcin(maxNgc), ncin, nmix
    real(8):: denmtr, dpin(maxNgc), eratio, hr, mix(maxNgc), ur, xyz


    write(IOOUT, '(/, /)')

    caseOK = .true.
    Nonly = 0
    Nomit = 0
    Nsert = 0
    reacts = .false.
    Trace = 0
    Short = .false.
    Massf = .false.
    Debug(1:Ncol) = .false.
    Nplt = 0
    SIunit = .true.
    pltdat = .false.

    do
       do
          ! CALL INFREE TO READ DATASET
          call INFREE(readOK, cin, ncin, lcin, dpin)

          if (.not. readOK) return

          code = trim(cin(1))

          if (code == '    ') then
             cycle
          end if

          ! STORE PRODUCT NAMES FROM 'ONLY' DATASET
          if (code == 'only') then
             Nonly = min(maxNgc, ncin-1)
             forall(i = 1:Nonly) Prod(i) = cin(i+1)

             ! STORE CONDENSED PRODUCT NAMES FROM 'INSERT' DATASET
          else if (code == 'inse') then
             Nsert = min(20, ncin-1)
             forall(i = 1:Nsert) Ensert(i) = cin(i+1)

             ! STORE PRODUCT NAMES FROM 'OMIT' DATASET
          else if (code == 'omit') then
             ! CHECK OMIT DATASET
             Nomit = min(maxNgc, ncin-1)
             forall(i = 1:Nomit) Omit(i) = cin(i+1)

             ! KEYWORD 'THER' READ
             ! CALL UTHERM TO CONVERT FORMATTED THERMODYNAMIC DATA
          else if (code == 'ther') then
             Newr = .true.
             rewind IOTHM
             call UTHERM(readOK)
             if (.not. readOK) then
                write(IOOUT, '(/" FATAL ERROR IN DATASET (INPUT)")')
                return
             end if

             ! KEYWORD 'TRAN' READ
             ! CALL UTRAN TO CONVERT FORMATTED TRANSPORT PROPERTIES
          else if (code == 'tran') then
             call UTRAN(readOK)
             if (.not. readOK) then
                write(IOOUT, '(/" FATAL ERROR IN DATASET (INPUT)")')
                return
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
                      forall(j = i+1:ncin, lcin(j) == i .and. int(dpin(j)) <= Ncol) Debug(int(dpin(j))) = .true.
                      forall(j = i+1:ncin, lcin(j) == i) lcin(j) = 0

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
             Pecwt(1:maxR) = -1

             i = 1
             do while (i < ncin)
                i = i + 1

                if (lcin(i) == 0) then
                   cycle
                end if

                if (lcin(i) > 0) then
                   write(IOOUT, '(/" WARNING!!  LITERAL EXPECTED FOR ", a15, "(INPUT)")') cin(i)
                   cycle
                end if

                cx15 = cin(i)
                cx1 = cx15(:1)
                cx2 = cx15(:2)
                cx3 = cx15(:3)
                cx4 = cx15(:4)

                ! NEW REACTANT
                if (cx2 == 'na' .or. cx2 == 'ox' .or. cx2 == 'fu') then
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

                else
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

                      ! LOOK FOR TEMPERATURES
                   else if (cx1 == 't') then
                      if (lcin(i+1) > 0) then
                         i = i + 1
                         Rtemp(Nreac) = dpin(i)

                         if (lcin(i-1) < 1) then
                            if (index(cx15, 'r') > 0) Rtemp(Nreac) = Rtemp(Nreac) / 1.8d0
                            if (index(cx15, 'c') > 0) Rtemp(Nreac) = Rtemp(Nreac) + 273.15d0
                            if (index(cx15, 'f') > 0) Rtemp(Nreac) = (Rtemp(Nreac) - 32) / 1.8d0 + 273.15d0
                         end if

                      else
                         write(IOOUT, '(/" REACTANT TEMPERATURE MISSING (INPUT) ")')
                         caseOK = .false.
                      end if

                      ! LOOK FOR ENTHALPY
                   else if (cx1 == 'h' .or. cx1 == 'u') then
                      Energy(Nreac) = cx15

                      if (lcin(i+1) > 0) then
                         i = i + 1
                         Enth(Nreac) = dpin(i) * 1000 / Rr

                         if (index(cin(i-1), 'c') > 0) Enth(Nreac) = Enth(Nreac) * 4.184d0
                         if (index(cin(i-1), 'k') > 0) Enth(Nreac) = Enth(Nreac) * 1000
                      end if

                      ! LOOK FOR DENSITY
                   else if (cx3 == 'rho' .or. cx3 == 'den') then
                      if (lcin(i+1) > 0) then
                         i = i + 1
                         Dens(Nreac) = dpin(i)
                         if (index(cx15, 'kg') > 0) Dens(Nreac) = Dens(Nreac) / 1000
                      end if

                      ! CHECK FOR CHEMICAL SYMBOLS IN EXPLODED FORMULA
                   else if ((lcin(i) == -1 .or. lcin(i) == -2) .and. index(uc, cx1) > 0) then
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
                         Rnum(Nreac, ifrmla) = 1
                      end if

                      i = i + 1

                   else
                      write(IOOUT, '(/" WARNING!! ", a15, " NOT RECOGNIZED (INPUT)")') cin(i)
                   end if
                end if
             end do

             ! SORT AND STORE INPUT FROM 'PROB' DATASET
          else if (code == 'prob') then
             Case = ' '
             P(1:maxPv) = 0
             V(1:maxPv) = 0
             T(1:maxT) = 0
             P(1) = 1
             Trace = 0
             Lsave = 0
             R = Rr / 4184
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
             forall(i = 1:Ncol)
                Pcp(i) = 0
                Pcp(i+Ncol) = 0
                Supar(i) = 0
                Subar(i) = 0
                Mach1(i) = 0
                U1(i) = 0
             end forall
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
             do i = 2, ncin
                if (lcin(i) < 0) then
                   if (any(lcin(i+1:ncin) == i)) cycle

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
             end do

             iv = 2
             Nof = 0
             exit

          else if (code(1:3) == 'end') then
             if (Shock) then
                if (incd .and. Froz) Incdfz = .true.
                if (incd .and. Eql) Incdeq = .true.
                if (refl .and. Froz) Reflfz = .true.
                if (refl .and. Eql) Refleq = .true.
             end if

             Hsub0 = min(hr, ur)
             Size = 0

             if (hr > 0.9d30) hr = 0
             if (ur > 0.9d30) ur = 0
             if (Trnspt) Viscns = 0.3125 * sqrt(1.e5 * Boltz / (pi * Avgdr))
             if (SIunit) R = Rr / 1000
             if (Detn .or. Shock) Newr = .true.

             if (.not. Short) then
                write(IOOUT, '(/" OPTIONS: TP=", l1, "  HP=", l1, "  SP=", l1, "  TV=", l1, &
                     & "  UV=", l1, "  SV=", l1, "  DETN=", l1, "  SHOCK=", l1, &
                     & "  REFL=", l1, "  INCD=", l1, /" RKT=", l1, "  FROZ=", l1, &
                     & "  EQL=", l1, "  IONS=", l1, "  SIUNIT=", l1, "  DEBUGF=", l1, &
                     & "  SHKDBG=", l1, "  DETDBG=", l1, "  TRNSPT=", l1)') &
                     Tp, (Hp .and. .not. Vol), Sp, (Tp .and. Vol), (Hp .and. Vol), (Sp .and. Vol), Detn, Shock, refl, &
                     incd, Rkt, Froz, Eql, Ions, SIunit, Debugf, Shkdbg, Detdbg, Trnspt
                if (T(1) > 0) write(IOOUT, '(/" T,K =", 7f11.4)') (T(jj), jj = 1, Nt)
                write(IOOUT, '(/1p, " TRACE=", e9.2, "  S/R=", e13.6, "  H/R=", e13.6, "  U/R=",  e13.6)') Trace, S0, hr, ur
                if (Np > 0 .and. Vol) write(IOOUT, '(/" SPECIFIC VOLUME,M**3/KG =", 1p, (4e14.7))') (V(jj)*1.d-05, jj = 1, Np)
             end if

             if (Rkt) then
                if (Nt == 0) Hp = .true.
                if (.not. Short) then
                   write(IOOUT, '(/" Pc,BAR =", 7f13.6)') (P(jj), jj = 1, Np)
                   write(IOOUT, '(/" Pc/P =", 9f11.4)') (Pcp(jj), jj = 1, Npp)
                   write(IOOUT, '(/" SUBSONIC AREA RATIOS =", (5f11.4))') (Subar(i), i = 1, Nsub)
                   write(IOOUT, '(/" SUPERSONIC AREA RATIOS =", (5f11.4))') (Supar(i), i = 1, Nsup)
                   write(IOOUT, '(/" NFZ=", i3, 1p, "  Mdot/Ac=", e13.6, "  Ac/At=", e13.6)') Nfz, Ma, Acat
                end if
             else
                if (.not. Vol .and. .not. Short) write(IOOUT, '(/" P,BAR =", 7f13.6)') (P(jj), jj = 1, Np)
             end if

             if (reacts) call REACT
             if (Nreac == 0 .or. Nlm <= 0) then
                write(IOOUT, '(/" ERROR IN REACTANTS DATASET (INPUT)")')
                caseOK = .false.
                write(IOOUT, '(/" FATAL ERROR IN DATASET (INPUT)")')
                return
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
                   return
                end if

             else if (phi .or. eqrats) then
                do i = 1, Nof
                   eratio = Oxf(i)
                   if (eqrats) then
                      xyz = -eratio * Vmin(2) - Vpls(2)
                      denmtr = eratio * Vmin(1) + Vpls(1)
                   else
                      xyz = -Vmin(2) - Vpls(2)
                      denmtr = eratio * (Vmin(1) + Vpls(1))
                   end if

                   if (abs(denmtr) < 1.d-30) then
                      caseOK = .false.
                      write(IOOUT, '(/" UNABLE TO PROCESS EQUIVALENCE RATIO =", E11.4, "(INPUT)")') eratio
                      write(IOOUT, '(/" FATAL ERROR IN DATASET (INPUT)")')
                      return
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
             return

          else
             write(IOOUT, '(/" WARNING!!  A KEYWORD IS MISSING (INPUT)")')
          end if

       end do

       ! PROCESS NUMERICAL DATA FOLLOWING 'PROB' LITERALS
       do
          in = 0
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
             if (iv < ncin) cycle
             exit
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
                      if (index(cx15, 'r') > 0) T(i) = T(i) / 1.8d0
                      if (index(cx15, 'c') > 0) T(i) = T(i) + 273.15d0
                      if (index(cx15, 'f') > 0) T(i) = (T(i) - 32) / 1.8d0 + 273.15d0
                   end if
                end if
             end do

          else if ((cx2 == 'pc' .or. cx2 == 'pi') .and. index(cx15(3:15), 'p') > 0 .and. index(cx15, 'psi') == 0) then
             Npp = nmix
             if (nmix > 2*Ncol) then
                Npp = 2 * Ncol
                write(IOOUT, '(/" NOTE!! MAXIMUM NUMBER OF ASSIGNED ", a5, " VALUES IS", I3, " (INPUT)", /)') 'pcp', Npp
             end if
             Pcp(1:Npp) = mix(1:Npp)

          else if (cx1 == 'p' .and. cx3 /= 'phi') then
             Np = nmix
             if (nmix > maxPv) then
                Np = maxPv
                write(IOOUT, '(/" NOTE!! MAXIMUM NUMBER OF ASSIGNED ", a5, " VALUES IS", I3, " (INPUT)", /)') 'p', Np
             end if

             do i = 1, Np
                P(i) = mix(i)
                if (index(cx15, 'psi') /= 0) then
                   P(i) = P(i) / 14.696006d0

                else if (index(cx15, 'mmh') /= 0) then
                   P(i) = P(i) / 760

                else if (index(cx15, 'atm') == 0) then
                   cycle
                end if

                P(i) = P(i) * 1.01325d0
             end do

          else if (cx3 == 'rho') then
             xyz = 1.d02
             if (index(cx15, 'kg') /= 0) xyz = 1.d05
             Np = nmix
             if (nmix > maxPv) then
                Np = maxPv
                write(IOOUT, '(/" NOTE!! MAXIMUM NUMBER OF ASSIGNED ", a5, " VALUES IS", i3, " (INPUT)", /)') 'rho', Np
             end if
             V(1:Np) = xyz / mix(1:Np)

          else if (cx1 == 'v') then
             xyz = 1.d02
             if (index(cx15, 'kg') /= 0) xyz = 1.d05
             Np = nmix
             if (nmix > maxPv) then
                Np = maxPv
                write(IOOUT, '(/" NOTE!! MAXIMUM NUMBER OF ASSIGNED ", a5, " VALUES IS", i3, " (INPUT)", /)') 'v', Np
             end if
             V(1:Np) = mix(1:Np) * xyz

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
                write(IOOUT, '(/" NOTE!! MAXIMUM NUMBER OF ASSIGNED ", a5, " VALUES IS", i3, " (INPUT)", /)') 'u1', Nsk
             end if
             U1(1:Nsk) = mix(1:Nsk)

          else if (cx4 == 'mach') then
             Nsk = nmix
             if (nmix > Ncol) then
                Nsk = Ncol
                write(IOOUT, '(/" NOTE!! MAXIMUM NUMBER OF ASSIGNED ", a5, " VALUES IS", i3, " (INPUT)", /)') 'mach1', Nsk
             end if
             Mach1(1:Nsk) = mix(1:Nsk)

          else if (cx3 == 'sub') then
             Nsub = nmix
             if (nmix > 13) then
                Nsub = 13
                write(IOOUT, '(/" NOTE!! MAXIMUM NUMBER OF ASSIGNED ", a5, " VALUES IS", i3, " (INPUT)", /)') 'subar', Nsub
             end if
             Subar(1:Nsub) = mix(1:Nsub)

          else if (cx3 == 'sup') then
             Nsup = nmix
             if (nmix > 13) then
                Nsup = 13
                write(IOOUT, '(/" NOTE!! MAXIMUM NUMBER OF ASSIGNED ", a5, " VALUES IS", i3, " (INPUT)", /)') 'supar', Nsup
             end if
             Supar(1:Nsup) = mix(1:Nsup)

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
                write(IOOUT, '(/" NOTE!! MAXIMUM NUMBER OF ASSIGNED ", a5, " VALUES IS", i3, " (INPUT)", /)') 'o/f', Nof
             end if
             Oxf(1:Nof) = mix(1:Nof)

             if (cx3 == 'phi') then
                phi = .true.

             else if (cx1 == 'r') then
                eqrats = .true.

             else if (cx3 == 'f/a') then
                forall(k = 1:Nof, Oxf(k) > 0) Oxf(k) = 1/Oxf(k)

             else if (cx4 == '%fue') then
                forall(k = 1:Nof, Oxf(k) > 0) Oxf(k) = (100 - Oxf(k)) / Oxf(k)
             end if

          else
             write(IOOUT, '("  WARNING!!  DID NOT RECOGNIZE ", a15, " (INPUT)"/)') cx15
          end if

          if (iv >= ncin) exit
       end do
    end do

    return
  end subroutine INPUT


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
    use mod_legacy_cea
    implicit none

    ! DUMMY ARGUMENTS
    character(15), intent(out):: Cin(maxNgc)
    integer, intent(out):: Ncin
    integer, intent(out):: Lcin(maxNgc)
    logical, intent(out):: readOK
    real(8), intent(out):: Dpin(maxNgc)

    ! LOCAL VARIABLES
    character(1), parameter:: nums(13) = ['+', '-', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '.']
    character(2), parameter:: numg(24) = &
         [character(2):: '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', &
         '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24']
    character(1):: ch1(132), cx_
    character(24):: cnum
    character(3):: fmtl(3)
    character(4):: w1
    integer:: i, ich1, j, kcin, nb, nch1, nx

    fmtl = [character(3):: '(g', '16', '.0)']

    Ncin = 1
    Lcin(1) = 0
    kcin = 0
    Dpin(1) = 0

    do
       nb = 1
       nx = 0
       cnum = ' '
       Cin(Ncin) = ' '
       ch1(1) = ' '
       nch1 = 1

       ! READ CHARACTERS, ONE AT A TIME
       read(IOINP, '(132a1)', END = 500, ERR = 500) ch1

       ! FIND FIRST AND LAST NON-BLANK CHARACTER
       do i = 132, 1, - 1
          nch1 = i
          if (ch1(i) /= ' ' .and. ch1(i) /= '	') exit
       end do

       do i = 1, nch1
          ich1 = i
          if (ch1(i) /= ' ' .and. ch1(i) /= '	') exit
       end do

       if (nch1 == 1 .or. ch1(ich1) == '#' .or. ch1(ich1) == '!') then
          write(IOOUT, '(1x, 80a1)') (ch1(i), i = 1, nch1)
          cycle
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

       do i = ich1, nch1
          cx_ = ch1(i)

          ! LOOK FOR DELIMITER STRINGS
          if (cx_ == ',' .and. (Lcin(Ncin) > 0 .or. nx == 0)) cx_ = ' '
          if (cx_ == '=' .and. (Lcin(Ncin) < 0 .or. nx == 0)) cx_ = ' '
          if (cx_ /= ' ' .and. cx_ /= '	') then

             ! LOOK FOR CHARACTER STRINGS
             nx = nx + 1
             if (Ncin > 1) then
                cnum(nx:nx) = cx_
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

310             nb = 1
             end if

             if (i < nch1 .or. Lcin(Ncin) < 0) cycle
          end if

          if (nb == 1. .and. nx > 0) then
             if (Ncin > 0 .and. Lcin(Ncin) > 0) then

                ! CONVERT NUMERIC CHARACTER STRINGS TO real(8) VARIABLES (DPIN)
                fmtl(2) = numg(min(24, nx))

                ! INTERNAL READ TO CONVERT TO NUMERIC
                read(cnum, fmtl, ERR=320) Dpin(Ncin)
             end if

             go to 340

320          if (Cin(Ncin-1)(:4) /= 'case') write(IOOUT, '(/" WARNING!!  UNACCEPTABLE NUMBER ", a15, " (INFREE)")') Cin(i)
             Lcin(Ncin) = 0

340          Ncin = Ncin + 1
             Cin(Ncin) = ' '
             Lcin(Ncin) = 0
             Dpin(Ncin) = 0
             nx = 0
             cnum = ' '
          end if

          nb = nb + 1

       end do

       if (nx > 0) then
          Ncin = Ncin + 1
          Lcin(Ncin) = 0
          Dpin(Ncin) = 0
       end if

    end do

500 readOK = .false.

    return
  end subroutine INFREE


  subroutine OUT1
    !***********************************************************************
    ! OUT1 WRITES REACTANT AND FUEL-OXIDANT RATIO INformatION.
    !
    ! NOTE - ROCKET, SHOCK, AND DETON PROBLEMS HAVE ADDITIONAL OUTPUT.
    !***********************************************************************
    use mod_legacy_cea
    implicit none

    integer:: i, n
    real(8):: rho

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

    phi = 0
    tem = (Vpls(1) + Vmin(1)) * Oxfl
    if (ABS(tem) >= 1.d-3) phi = -(Vmin(2) + Vpls(2)) / tem

    if (Fox(1) == 'NAME') then
       pfuel = 0
    else
       pfuel = 100 / (1 + Oxfl)
    end if

    if (any(Rh /= 0)) then
       if (any(Rh == 0)) then
          rho = max(Rh(1), Rh(2))
       else
          rho = (Oxfl + 1) * Rh(1) * Rh(2) / (Rh(1) + Oxfl * Rh(2))
       end if
       if (SIunit) then
          rho = rho * 1000
          write(IOOUT, '(/" REACTANT DENSITY=", F8.2, " KG/CU M")') rho
       else
          write(IOOUT, '(/" REACTANT DENSITY=", F8.4, " G/CC")') rho
       end if
    end if

    write(IOOUT, '(/" O/F=", F11.5, 2X, "%FUEL=", F10.6, 2X, "R,EQ.RATIO=", F9.6, 2X, "PHI,EQ.RATIO=", F9.6)') &
         Oxfl, pfuel, Eqrat, phi
    return
  end subroutine OUT1


  subroutine OUT2
    !***********************************************************************
    ! OUT2 WRITES THERMODYNAMIC PROPERTIES.
    !***********************************************************************
    use mod_legacy_cea
    implicit none

    character(15):: fgi, fh, fp, frh, fs, fu
    integer:: i, j
    integer, save:: mcp, mdvp, mdvt, meq, mfa, mg, mgam, mh, mie, mm, mmw, mof, mp, mpf, mph, mrho, ms, mson, mt
    real(8):: pfactor, vnum

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
             if (index(Pltvar(i)(3:), 'fz') /= 0 .or. index(Pltvar(i)(3:), 'fr') /= 0) then
                mpnf = i
             else
                mpn = i
             end if
          else if (Pltvar(i)(:4) == 'cond') then
             if (index(Pltvar(i)(3:), 'fz') /= 0 .or. index(Pltvar(i)(3:), 'fr') /= 0) then
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
       pfactor = 1
       fp = 'P, BAR'
       vnum = 1.d05
       frh = 'RHO, KG/CU M'
       fh = 'H, KJ/KG'
       fu = 'U, KJ/KG'
       fgi = 'G, KJ/KG'
       fs = 'S, KJ/(KG)(K)'
       fc = 'Cp, KJ/(KG)(K)'
    else
       pfactor = 1 / 1.01325d0
       fp = 'P, ATM'
       vnum = 100
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
       X(i) = Ppp(i) * pfactor
       if (Nplt /= 0 .and. i > ione) then
          if (mp > 0) Pltout(i + Iplt - ione, mp) = X(i)
          if (mt > 0) Pltout(i + Iplt - ione, mt) = Ttt(i)
       end if
    end do

    write(IOOUT, Fmt) fp, (X(j), j = 1, Npt)

    ! TEMPERATURE
    Fmt(4) = '13'
    Fmt(5) = ' '
    Fmt(7) = '2,'
    write(IOOUT, Fmt) 'T, K            ', (Ttt(j), j = 1, Npt)

    ! DENSITY
    do i = 1, Npt
       if (Vlm(i) /= 0) X(i) = vnum / Vlm(i)
       if (Nplt /= 0 .and. i > ione .and. mrho > 0) Pltout(i+Iplt-ione, mrho) = X(i)
    end do
    call EFMT(Fmt(4), frh, X)

    ! ENTHALPY
    do i = 1, Npt
       X(i) = Hsum(i) * R
       if (Nplt /= 0 .and. i > ione .and. mh > 0) Pltout(i+Iplt-ione, mh) = X(i)
    end do
    Fmt(4) = Fmt(6)
    call VARFMT(X)
    write(IOOUT, Fmt) fh, (X(j), j = 1, Npt)

    ! INTERNAL ENERGY
    do i = 1, Npt
       X(i) = (Hsum(i) - Ppp(i) * Vlm(i) / Rr) * R
       if (Nplt /= 0 .and. i > ione .and. mie > 0) Pltout(i+Iplt-ione, mie) = X(i)
    end do
    call VARFMT(X)
    write(IOOUT, Fmt) fu, (X(j), j = 1, Npt)

    ! GIBBS ENERGY
    do i = 1, Npt
       X(i) = (Hsum(i) - Ttt(i) * Ssum(i)) * R
       if (Nplt /= 0 .and. i > ione) then
          if (mg > 0) Pltout(i+Iplt-ione, mg) = X(i)
          if (mm > 0) Pltout(i+Iplt-ione, mm) = Wm(i)
          if (mmw > 0) Pltout(i+Iplt-ione, mmw) = 1 / Totn(i)
          if (ms > 0) Pltout(i+Iplt-ione, ms) = Ssum(i) * R
          if (mcp > 0) Pltout(i+Iplt-ione, mcp) = Cpr(i) * R
          if (mgam > 0) Pltout(i+Iplt-ione, mgam) = Gammas(i)
          if (mdvt > 0) Pltout(i+Iplt-ione, mdvt) = Dlvtp(i)
          if (mdvp > 0) Pltout(i+Iplt-ione, mdvp) = Dlvpt(i)
       end if
    end do
    call VARFMT(X)
    write(IOOUT, Fmt) fgi, (X(j), j = 1, Npt)

    ! ENTROPY
    Fmt(4) = '13'
    Fmt(5) = ' '
    Fmt(7) = '4,'
    write(IOOUT, Fmt) fs, (Ssum(j) * R, j = 1, Npt)
    write(IOOUT, *)

    ! MOLECULAR WEIGHT
    Fmt(7) = '3,'
    write(IOOUT, Fmt) 'M, (1/n)        ', (Wm(j), j = 1, Npt)
    if (.not. Gonly) write(IOOUT, Fmt) 'MW, MOL WT      ', (1/Totn(j), j = 1, Npt)

    ! (DLV/DLP)T
    Fmt(7) = '5,'
    if (Eql) write(IOOUT, Fmt) '(dLV/dLP)t      ', (Dlvpt(j), j = 1, Npt)

    ! (DLV/DLT)P
    Fmt(7) = '4,'
    if (Eql) write(IOOUT, Fmt) '(dLV/dLT)p      ', (Dlvtp(j), j = 1, Npt)

    ! HEAT CAPACITY
    write(IOOUT, Fmt) fc, (Cpr(j) * R, j = 1, Npt)

    ! GAMMA(S)
    Fmt(7) = '4,'
    write(IOOUT, Fmt) 'GAMMAs          ', (Gammas(j), j = 1, Npt)

    ! SONIC VELOCITY
    Fmt(7) = '1,'
    do i = 1, Npt
       Sonvel(i) = sqrt(Rr * Gammas(i) * Ttt(i) / Wm(i))
       if (Nplt /= 0 .and. i > ione .and. mson > 0) Pltout(i+Iplt-ione, mson) = Sonvel(i)
    end do
    write(IOOUT, Fmt) 'SON VEL,M/SEC   ', (Sonvel(j), j = 1, Npt)
    return
  end subroutine OUT2


  subroutine OUT3
    !***********************************************************************
    ! OUT3 WRITES MOLE FRACTIONS.
    !***********************************************************************
    use mod_legacy_cea
    implicit none

    character(4):: mamo
    logical:: kOK
    integer:: i, k, m, im, kin
    integer, save:: notuse
    real(8):: tra

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
          else
             do m = 1, Nplt
                im = 0
                if (Pltvar(m) == Prod(k) .or. '*' // Pltvar(m) == Prod(k)) then
                   im = m
                   exit
                end if
             end do
          end if

          kin = 0
          do i = 1, Npt
             if (Massf) then
                tem = Mw(k)
             else
                tem = 1 / Totn(i)
             end if
             if (k <= Ng) then
                X(i) = En(k, i) * tem
             else
                if (Prod(k) /= Prod(k-1)) X(i) = 0
                if (En(k, i) > 0) X(i) = En(k, i) * tem
             end if
             if (Nplt /= 0 .and. i > ione .and. im > 0) Pltout(i+Iplt-ione, im) = X(i)
             if (kOK .and. X(i) >= tra) kin = 1
          end do

          if (kin == 1) then
             if (Trace == 0) then
                write(IOOUT, '(1x, A15, F9.5, 12F9.5)') Prod(k), (X(i), i = 1, Npt)
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

    write(IOOUT, '(/"  * THERMODYNAMIC PROPERTIES FITTED TO", f7.0, "K")') Tg(4)
    if (.not. Short) then
       write(IOOUT, '(/"    PRODUCTS WHICH WERE CONSIDERED BUT WHOSE ", a4, &
            & " FRACTIONS", /"    WERE LESS THAN", 1pe13.6, &
            & " FOR ALL ASSIGNED CONDITIONS"/)') mamo, tra
       write(IOOUT, '(5(1x, a15))') (Omit(i), i = 1, notuse)
    end if

    if (.not. Moles) write(IOOUT, '(/" NOTE. WEIGHT FRACTION OF FUEL IN TOTAL FUELS AND OF", &
         & " OXIDANT IN TOTAL OXIDANTS")')
    return
  end subroutine OUT3


  subroutine OUT4
    !***********************************************************************
    ! OUT4 WRITES TRANSPORT PROPERTIES.
    !***********************************************************************
    use mod_legacy_cea
    implicit none

    integer:: i, j

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

    write(IOOUT, Fmt) 'VISC,MILLIPOISE', (Vis(j), j = 1, Npt)

    Fmt(4) = '13'
    Fmt(5) = ' '
    Fmt(7) = '4,'

    if (Eql) then
       write(IOOUT, '(/"  WITH EQUILIBRIUM REACTIONS"/)')
       ! SPECIFIC HEAT
       write(IOOUT, Fmt) fc, (Cpeql(j), j = 1, Npt)
       ! CONDUCTIVITY
       write(IOOUT, Fmt) 'CONDUCTIVITY    ', (Coneql(j), j = 1, Npt)
       ! PRANDTL NUMBER
       write(IOOUT, Fmt) 'PRANDTL NUMBER  ', (Preql(j), j = 1, Npt)
    end if

    write(IOOUT, '(/"  WITH FROZEN REACTIONS"/)')
    ! SPECIFIC HEAT
    write(IOOUT, Fmt) fc, (Cpfro(j), j = 1, Npt)
    ! CONDUCTIVITY
    write(IOOUT, Fmt) 'CONDUCTIVITY    ', (Confro(j), j = 1, Npt)
    ! PRANDTL NUMBER
    write(IOOUT, Fmt) 'PRANDTL NUMBER  ', (Prfro(j), j = 1, Npt)

    return
  end subroutine OUT4

end module mod_legacy_io
!!$
!!$
!!$  subroutine read_legacy_input(filename, problems)
!!$    character(*), intent(in):: filename
!!$    type(CEA_Problem), dimension(:), intent(out):: problems
!!$
!!$  end subroutine read_legacy_input
!!$
!!$  subroutine write_legacy_output(filename, problems)
!!$    character(*), intent(in):: filename
!!$    type(CEA_Problem), dimension(:), intent(in):: problems
!!$
!!$  end subroutine write_legacy_output
