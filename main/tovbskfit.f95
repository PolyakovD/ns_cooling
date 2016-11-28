! Integration of the Tolman-Oppenheimer-Volkoff equation
!   with the fit to one of the equations of state BSk19, BSk20, or BSk21
! This code should be linked with bskfit.f
! References: http://www.ioffe.ru/astro/NSG/BSk/
!   A.Y. Potekhin, A.F. Fantina, N. Chamel, J.M. Pearson, & S. Goriely,
!     2013, Astron. Astrophys. ...
! Contact: Alexander Potekhin <palex@astro.ioffe.ru>
! Last change: 14.03.13
    !  implicit double precision (A-H), double precision (O-Z)
    !  character KEY
    !  data EPSMASS/1.d-7/
!      write(*,'('' EPSMASS: ''$)')
!      read*,EPSMASS
    !  write(*, '('' BSk: 19, 20, or 21? ''$)')
    !  read*, KEOS
    !  write(*,110) KEOS, EPSMASS
    !1 continue
    !  write(*,'('' rho_c (g/cc) (<0 to stop: ''$)')
    !  read*, RHOC
    !  if (RHOC .le. 0.) stop
    !  call TOVINT(RHOC,KEOS,EPSMASS,R6,SMASS)
    !  write(*, 112)
    !  write(*, 111) RHOC, (R6 * 10.d0), SMASS
    !  goto 1
  !110 format(' Output of tovbskfit v.11.03.13'/'KEOS   EPS'/I3,1PE11.2)
  !111 format(1PE12.4,0PF9.3,F11.6)
  !112 format('rho_c [g/cm^3]  R [km]  M [MSun]')
      !end
      
      subroutine TOVINT(RHOC, KEOS, isDeep, EPSMASS, R6, SMASS, length, outputFile)
!                                                       Version 14.03.13
! Integration of the Tolman-Oppenheimer-Volkoff equation
!   with the fit to one of the equations of state BSk19, BSk20, or BSk21
!   from the stellar center with mass density RHOC .
! Input: RHOC - central grav.mass density [g/cc]
!        KEOS=19,20,21 - EOS number
!        isDeep - "1", if the model estimated to 10^{10} g/cc, "2", if the model estimated to 10^{8} g/cc
!        EPSMASS - relative accuracy in total stellar mass SMASS
! Output: R6 - stellar radius [10^6 cm]
!         SMASS - stellar mass / MSun
      implicit double precision (A-H), double precision (O-Z)
      
      integer, parameter :: fout = 101 
      character(*) :: outputFile
      integer(4) :: isDeep
      integer(4) :: length
      
      parameter (EPS = .001d0, EPS1 = .05d0)
      real*8 buffer(3, 4) ! r, rho, m_{r}
      
      real(8) :: surfaceRho
      
      if (isDeep .eq. 1) then
          surfaceRho = 1.0d+10
      else
          surfaceRho = 1.0d+8
      endif
      open(fout, file = outputFile)
      
      length = 1
      write(fout, '(e22.15, 5x, e22.15, 5x, e22.15)') 0.0d0, RHOC, 0.0d0
      RLG = dlog10(RHOC)
      call BSKfit(KEOS, RLG, PLG, CHIrho, Gamma)
      R6 = 0.
      SM = 0.
      P12 = 10.d0**(PLG - 12.d0)
      RHO = 10.d0**RLG
      SMASS = 0.
      H = EPS
    1 continue
      RHO0 = RHO
      SMASS0 = SMASS
      R60 = R6
    2 continue
      do ITRY = 1, 2
        do IBIN = 1, ITRY
           HH = H*.5d0
! 1:
           RLG = dlog10(RHO)
           call BSKfit(KEOS, RLG, PLG, CHIrho, Gamma)
           P12fit = 10.d0**(PLG - 12.d0)
           P12out = 3.5d2 * RHO
           P12 = P12fit + P12out
           CHIrho = (P12fit * CHIrho + P12out) / P12
           call SUBRHS(P12, RHO, R6, SMASS, DP12, DSMASS)
           DRHO = DP12 / P12 * RHO / CHIrho
           C1R = H * DRHO
           C1M = H * DSMASS
           
           if (ITRY .eq. 2) then
               buffer(1,1 + (IBIN - 1) * 2) = (R6 + HH / 2.0d0) * 1.0d+6
               buffer(2,1 + (IBIN - 1) * 2) = RHO + C1R 
               buffer(3,1 + (IBIN - 1) * 2) = SMASS + C1M
           endif
! 2:
           R6A = R6+HH
           RHOA = RHO+C1R*.5d0
           SMASSA = SMASS+C1M*.5d0
           RLGA = dlog10(RHOA)
           call BSKfit(KEOS, RLGA, PLGA, CHIrho, Gamma)
           P12fit = 10.d0**(PLGA - 12.d0)
           P12out = 3.5d2 * RHOA
           P12A = P12fit + P12out
           CHIrho = (P12fit * CHIrho + P12out) / P12A
           call SUBRHS(P12A, RHOA, R6A, SMASSA, DP12, DSMASS)
           DRHO = DP12/P12A*RHOA/CHIrho
           C2R = H * DRHO
           C2M = H * DSMASS
           
           if (ITRY .eq. 2) then
               buffer(1,2 + (IBIN - 1) * 2) = (R6A + HH / 2.0d0) * 1.0d+6
               buffer(2,2 + (IBIN - 1) * 2) = RHOA + C2R  
               buffer(3,2 + (IBIN - 1) * 2) = SMASSA + C2M
           endif
! 3:
           RHOA = RHO + C2R * .5d0
           SMASSA = SMASS + C2M * .5d0
           RLGA = dlog10(RHOA)
           call BSKfit(KEOS, RLGA, PLGA, CHIrho, Gamma)
           P12fit = 10.d0**(PLGA - 12.d0)
           P12out = 3.5d2 * RHOA
           P12A = P12fit + P12out
           CHIrho = (P12fit * CHIrho + P12out) / P12A
           call SUBRHS(P12A, RHOA, R6A, SMASSA, DP12, DSMASS)
           DRHO = DP12 / P12A * RHOA / CHIrho
           C3R = H * DRHO
           C3M = H * DSMASS
           
           !buffer(1,3) = (R6A + H) * 1.0d+6
           !buffer(2,3) = RHOA + C3R
! 4:
           R6A = R6+H
           RHOA = RHO+C3R
           SMASSA = SMASS + C3M
           RLGA = dlog10(RHOA)
           call BSKfit(KEOS, RLGA, PLGA, CHIrho, Gamma)
           P12fit = 10.d0**(PLGA - 12.d0)
           P12out = 3.5d2 * RHOA
           P12A = P12fit + P12out
           CHIrho = (P12fit * CHIrho + P12out) / P12A
           call SUBRHS (P12A, RHOA, R6A, SMASSA, DP12, DSMASS)
           DRHO = DP12 / P12A * RHOA / CHIrho
           C4R = H * DRHO
           C4M = H * DSMASS
          
! Final:
           RHOA = RHO + (C1R + 2.d0 * (C2R + C3R) + C4R) / 6.d0
           SMASSA = SMASS + (C1M + 2.d0 * (C2M + C3M) + C4M) / 6.d0
          if (ITRY .eq. 1) then
             RHO1 = RHOA
             SMASS1 = SMASSA
             H = HH ! the 2nd try will be with the diminished step
          else
             RHO = RHOA
             SMASS = SMASSA
             R6 = R6A
          endif
        enddo ! next half, if ITRY=2
      enddo ! next ITRY
! Check accuracy:
      DELTA = dmax1(dabs(RHO - RHO1) / RHO,(SMASS - SMASS1) / SMASS)
      if (DELTA .lt. EPSMASS) then ! sufficient accuracy
         do i = 1, 4
            if (buffer(2, i) > surfaceRho) then
                length = length + 1 
                write(fout, '(e22.15, 5x, e22.15, 5x, e22.15)') buffer(1, i), buffer(2, i), buffer(3, i)
            end if
         end do
         H = 4.d0 * H ! increase H for the next interval
      else ! decrease H and retry
         R6 = R60
         RHO = RHO0
         SMASS = SMASS0
         H = 1.5d0 * H * (EPSMASS / DELTA)**.2
         goto 2
      endif
      if (DRHO .ne. 0.) H = dmin1(H, EPS1 * RHO / dabs(DRHO))
      if (RHO .gt. 7.86d0) goto 1 ! next step
      close(fout)
      !write(*, *) 'rho_c [g/cm^3]  R [km]  M [MSun]'
      !write(*, '(1PE12.4,0PF9.3,F11.6)') RHOC, (R6 * 10.d0), SMASS
      return
      end

      subroutine SUBRHS(P12,RHO,R6,SMASS,DP12,DSMASS)
!                                                       Version 24.07.12
! Stems from SUBRHS7 v.08.11.07.
! Input: P12 - pressure P in Mbar
!        RHO - density in g.cc
!        R6 - radial coordinate r in 10^6 cm [= 10 km]
!        SMASS - m(r)/Msun
! Output: DP12= d P12 / d R6
!         DSMASS= d m(r) / d R6 in Msun
!         G14 - local gravitational acceleration in 10^{14} cm/s^2
      implicit double precision (A-H), double precision (O-Z)
      save
      if (R6 .gt. 0.) then
         GR_R = dsqrt(1.d0 - .295423 * SMASS / R6) ! radius (inverse vol.) corr.
         G14 = 1.32757 * SMASS / R6**2 / GR_R ! g(r) in 10^{14} cm/s^2
      else
         G14 = 0.
         GR_R = 1.d0
      endif
      if (SMASS .gt. 0.) then
         GR_G = 1.d0 + 7.0293d-24 * R6**3 * P12 / SMASS ! gravit.-accel.correction
      else
         GR_G = 1.d0
      endif
      Prhoc2 = 1.11265d-9 * P12 / RHO
      GR_H = 1.d0 + Prhoc2 ! enthalpy correction
      DP12 = -RHO * G14 * GR_H * GR_G / GR_R * 1.d8
      DSMASS = 6.3182d-15 * R6**2 * RHO
      return
      end