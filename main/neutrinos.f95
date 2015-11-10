module neutrinos
   
    real(8), private, parameter :: PI = 3.14159265d0
    real(8), private, parameter :: DECLN = 2.302585093d0
    real(8), private, parameter :: eV = 1.602176487d-12 ! erg/eV
    real(8), private, parameter :: AUM = 1822.9 ! a.m.u./m_e
    real(8), private, parameter :: AUD = 15819.3832d0 ! relativistic unit of density in g/cc
    real(8), private, parameter :: SBCREL = 0.1644933 ! Stefan-Boltzmann constant in rel.un.
    real(8), private, parameter :: UAMU = 1.660538782d-24 ! unified atomic mass unit in grams
    real(8), private, parameter :: UNIPR = 1.4217752544d13 ! relativistic unit of pressure in Mbar
    real(8), private, parameter :: UNIB = 44.14004916 ! rel.un.of magnetic field [TG]
    real(8), private, parameter :: GAMIMELT = 175. ! Gamma_i value at which the melting occurs
    real(8), private, parameter :: RSIMELT = 140. ! ion density parameter of quantum melting
    ! [approximate value, acc.to Jones & Ceperley (1996) PRL 76, 4572]
    real(8), private, parameter :: Ryd = 13.6057d0 ! Rydberg energy / eV
    real(8), private, parameter :: UN_T6 = .3157746d0 ! = 2*(Rydberg energy)/10^6 K = a.u. of T [MK]
    real(8), private, parameter :: BOHR = 137.03599968d0 ! inverse fine-structure constant
    real(8), private, parameter :: AURHO = 11.2058713d0 ! (atomic unit of mass density)*AUM in g/cc
    real(8), private, parameter :: BOHRAD = .52917720859d-8 ! Bohr radius in cm
    real(8), private, parameter :: COMPTON = 3.8615926459d-11 ! Compton length over 2\pi in cm
    real(8), private, parameter :: AUTIME = 2.4188843265d-17 ! atomic unit of time in s
    
    real(8), private, parameter :: C13 = 1.d0 / 3.d0, C53 = 5.d0 / 3.d0, C73 = 7.d0 / 3.d0
    real(8), private, parameter :: TwoPI = 2.d0 * PI, FourPI = 4.d0 * PI
    real(8), public, parameter :: UNITEMP = UN_T6 * BOHR**2 ! rel.un. of T [MK] = 5930
    real(8), private, parameter :: UNIOP = 1.d0 / (AUD * COMPTON) ! rel.un.of opacity [cm^2/g]=1.637d6
    real(8), private, parameter :: AUNUMDEN = 1.d0 / BOHRAD**3 ! a.u. of number density
    real(8), private, parameter :: UN_B12 = BOHR**2 / UNIB ! inverse a.u. of magnetic field [TG]
    
    real(8), private, parameter :: M_P = 1.6726231d-24   ! the rest proton mass 
    real(8), private, parameter :: M_N = 1.6749286d-24   ! the rest neutron mass 
    real(8), private, parameter :: N0 = 0.16d0 ! [fm -3] is the nucleon number density in standard nuclear matter
    real(8), private, parameter :: c = 3.0d10 ! the velocity of light in vacuum
    real(8), private, parameter :: nNu = 3.0d0 ! number of neutrino flavors
    real(8), private, parameter :: aNB = 4.17d0 ! coefficien for Cooper pairing of neutrons in core
    real(8), private, parameter :: aNA = 1.0d0 ! coefficien for Cooper pairing of neutrons in crust

contains   
      
      subroutine SINQpair(X, T, Qpair)
!                                                       Version 09.05.11
! Neutrino emission due to e-e+ annihilation: YKGH(2001):Eq.(22)
! Input: X - relativity parameter
!        T - temperature in relativistic units
! Output: Qpair - neutrino emission rate in erg/(s cm^3)
      implicit double precision (A-H), double precision (O-Z)
      !include 'const.inc' ! i.o.PI,TwoPI 09.05.11
      parameter(Qc = 1.023d23, Cplus2 = 1.7037)
      parameter(TwoPI2 = 2.d0 * PI**2, Qc0 = Qc / 36.d0 / PI, Cmin2 = Cplus2 - 1.5)
      Y = dsqrt(1.d0 + X**2)
      EY = T / TwoPI2 * dexp(-(1.d0 + Y) / T)
      PHIm1 = dsqrt(T * (TwoPI + T * (7.662 + 1.92 * T)) / (1. + .48 * T)) * EY ! YKGH:27
      PHI0 = dsqrt(T * (TwoPI + T * (23.61 + T * (32.11 + T * 16.)))) * EY
      PHIp1 = dsqrt(T * (TwoPI + T * (42.44 + T * (140.8 + T * (265.2 + T * (287.9 + T * 144.)))))) * EY ! YKGH:27
      T2 = T**2
      PHIp2 = dsqrt(T * (TwoPI + 61.33 * T + 321.9 * T2 + 1153. * T * T2 + 2624. * T2**2 + &
                    4468. * T**5 + T2**3 * (4600. + 2304. * T))) * EY
      CLN = dlog(X + Y)
      Um1 = (Y * X - CLN) / TwoPI2
      U0 = X**3 / (3. * PI**2)
      Up1 = (Y * X * (X**2 + Y**2) - CLN) / TwoPI2 / 4.
      Up2 = X**3 * (5. + 3. * X**2) / (15. * PI**2)
      Q = 8. * (PHIp1 * Up2 + PHIp2 * Up1) -2. *(PHIm1 * Up2 + PHIp2 * Um1) + &
          7. * (PHI0 * Up1 + PHIp1 * U0) + 5. *(PHI0 * Um1 + PHIm1 * U0)
      Q = Cplus2 * Q + 9. * Cmin2 * (PHI0 * (Up1 + Um1) + (PHIm1 + PHIp1) * U0)
      Qpair = Qc0 * Q
      return
      end subroutine
        
      subroutine SINQpl(X, T, Qpl)
!                                                       Version 04.05.11
! Neutrino emission due to plasmon decay: YKGH(2001):Eq.(38)
! Input: X - relativity parameter
!        T - temperature in relativistic units
! Output: Qpl - neutrino emission rate in erg/(s cm^3)
      implicit double precision (A-H), double precision (O-Z)
      !include 'const.inc' ! i.o."PI=3.14159265d0,BOHR=137.036" 04.05.11
      parameter(Qc = 1.023d23, CV2 = 0.9537)
      parameter(Cfp = 4. / (3. * BOHR * PI), Qc1 = Qc * BOHR / (96.d0 * PI**4) * CV2) ! 4/3 alpha_f/pi  
      Y = dsqrt(1.d0 + X**2)
      fp = dsqrt(Cfp * X**3 / Y) / T
      Qpl = Qc1 * T**9 * fp**6 * (16.23 + 4.604 * dsqrt(fp) * fp) * dexp(-fp)
      return
      end subroutine

      subroutine SINQsyn(X, T, B12, Qsyn)
!                                                       Version 09.05.11
! Neutrino electron-synchrotron: YKGH(2001):Eq.(56)
! Input: X - relativity parameter
!        T - temperature in relativistic units
!        B12 - m.f. / 10^12 G
! Output: Qsyn - neutrino emission rate in erg/(s cm^3)
      implicit double precision (A-H), double precision (O-Z)
      !include 'const.inc' ! PI,UNIB 09.05.11
      parameter(TwoThird = 2.d0 / 3.d0, zeta5 = 1.036928d0, Cplus2 = 1.7037)
      parameter(CSAB = 27.d0 / (512.d0 * PI**2 * zeta5), Cmin2 = Cplus2 - 1.5)
      B = B12 / UNIB ! m.f. in rel.un.
      B13 = B12 / 10.
      Y = dsqrt(1.d0 + X**2)
      z = B / Y / T ! =T_B/T
      T9 = T * 5.93
      QsynB = 9.04d14 * B13**2 * T9**5
      XI = 1.5 * z * X**3
      XI23 = XI**TwoThird
      Y1 = dsqrt((1.d0 + 3172. * XI23)**TwoThird - 1.d0)**3
      Y2 = dsqrt((1.d0 + 172.2 * XI23)**TwoThird - 1.d0)**3
      Fplus = 44.01 * ((1. + 3.675d-4 * Y1) / (1. + 2.036d-4 * Y1 + 7.405d-8 * Y1**2)**2)**2  
      Fmin = 36.97 * (1. + Y2 * (1.436d-2 + Y2 * (1.024d-5 + Y2 * 7.647d-8))) / &
            (1. + 3.356d-3 * Y2 + 1.536 * Y2**2)**5
      SAB = CSAB * XI**4 * (Fplus - Cmin2 / Cplus2 * Fmin)
      SBC = dexp(-z / 2.) * (1. + z * (.4228 + z * (.1014 + z * .00624))) / &
           (1. + .4535 * z**TwoThird + z * (.03008 + z * (-.05043 + z * .004314)))
      Qsyn = QsynB * SAB * SBC
      return
      end subroutine

    subroutine SINQbr(RHO, T6, Qbr)
! Neutrino electron-bremsstrahlung: YKGH(2001):Eq.(76)
! Input: RHO - density [g/cc]
!        T6 - temperature [MK]
! Output: Qbr - neutrino emission rate in erg/(s cm^3)
    implicit double precision (A-H), double precision (O-Z)
    R = dlog10(RHO) - 12.d0 ! lg rho_{12}
    TAU = dlog10(T6) - 2.d0 ! lg T_8
    R2 = R**2
    if (R .lt. -3.) R2 = -(6. * R + 9.) ! linear extrapol.beyond lower limit
    TAU2 = TAU**2
    if (T6 .lt. 50.) then ! linear extrapol.beyond lower limit
       TAU0 = -.30103
       TAU2 = TAU0 * (2. * TAU - TAU0)
    elseif (T6 .gt. 2.d3) then ! linear extrapol.beyond higher limit
       TAU0 = 1.30103
       TAU2 = TAU0 * (2. * TAU - TAU0)
    endif
    QbrLG = 11.204 + 7.304 * TAU+ .2976 * R -.37 * TAU2 + .188 * TAU * R - .103 * R2 + &
            .0547 * TAU2 * R - 6.77 * dlog10(1.d0 + .228 * RHO / 2.8d14)
    Qbr = 10.d0**QbrLG
    return
    end subroutine
    
    ! Neutrino emission due to nucleon direct Urka process: YKGH(2001):Eq.(120)
    ! Input: nE - electron number density in fm^{-3}
    !        T9 = temperature / 10^{9}
    !        EFMp,EFMn - effective proton and neutron mass factors(mP*/mP) ,(mN*/mN), respectively 
	!        RD - superfluid reduction factor
    ! Output: QDUP - neutrino emission rate in erg/(s cm^3)
    subroutine getQDUP(nE, T9, EFMp, EFMn, RD, QDUP)
        implicit none
        real(8) :: nE, T9, EFMp, EFMn, RD
        real(8) :: QDUP
        QDUP = RD * 4.0d27 * (nE / N0)**(1.0d0 / 3.0d0) * EFMn * EFMp * M_P / M_N * T9**(6.0d0)
    end subroutine
    
    ! Neutrino emission due to nucleon modified Urka process: 
    ! YKGH(2001):Eq.(140), Eq.(142) with theta_{Mp} = 1
    ! Input: nP - proton number density in fm^{-3}
    !        T9 = temperature / 10^{9}
    !        EFMp, EFMn - effective proton and neutron mass factors(mP*/mP) ,(mN*/mN), respectively 
    !        pFE, pFP, pFN - electron, proton and neutron Fermi momentum, respectively 
	!        RMn - superfluid reduction factor for neutron branch
	!        RMp - superfluid reduction factor for proton branch
    ! Output: QMUP - neutrino emission rate in erg/(s cm^3)
    subroutine getQMUP(nP, T9, EFMp, EFMn, pFE, pFP, pFN, RMn, RMp, QMUP)
        implicit none
        real(8) :: nP, T9, EFMp, EFMn, pFE, pFP, pFN, RMn, RMp 
        real(8) :: QMUP
        real(8), parameter :: alphaN = 1.13d0, betaN = 0.68d0 
        ! Neutron branch
        QMUP = RMn * 8.1d21 * EFMn**(3.0d0) * EFMp * (nP / N0)**(1.0d0 / 3.0d0) * T9**(8.0d0) * alphaN * betaN 
        ! Append proton branch
        QMUP = QMUP * (1.0d0 + RMp * (EFMp / EFMn)**(2.0d0) * (pFE + 3.0d0 * pFP - pFN)**(2.0d0) / (8.0d0 * pFE * pFP))
    end subroutine
    
    ! Neutrino emission due to baryon-baryon collisions:
    ! YKGH(2001):Eq.(165), Eq.(166), Eq.(167)
    ! Input: nN, nP - neutron and proton number density in fm^{-3}, , respectively    
    !        T9 = temperature / 10^{9}
    !        EFMn, EFMp - effective neutron and proton mass factors(mP*/mP) ,(mN*/mN), respectively 
	!        Rnn, Rnp, Rpp - superfluid reduction factor for nn, np, pp collisions, respectively
    ! Output: QBBR - neutrino emission rate in erg/(s cm^3)
    subroutine getQBBR(nN, nP, T9, EFMn, EFMp, Rnn, Rnp, Rpp, QBBR)
        implicit none
        real(8) :: nN, nP, T9, EFMp, EFMn, Rnn, Rnp, Rpp
        real(8) :: QBBR
        real(8), parameter :: alphaNN = 0.59d0, alphaNP = 1.06d0, alphaPP = 0.11d0
        real(8), parameter :: betaNN = 0.56d0, betaNP = 0.66d0, betaPP = 0.7d0
        real(8), parameter :: nNu = 3.0d0
        ! Neutron-neutron collisions
        QBBR = Rnn * 7.5d19 * EFMn**(4.0d0) * (nN / N0)**(1.0d0 / 3.0d0) * & 
               alphaNN * betaNN * nNu * T9**(8.0d0) 
        ! Neutron-proton collisions
        QBBR = QBBR + Rnp * 1.5d20 * (EFMn * EFMp)**(2.0d0) * (nP / N0)**(1.0d0 / 3.0d0) * & 
               alphaNP * betaNP * nNu * T9**(8.0d0)
        ! Proton-proton collisions
        QBBR = QBBR + Rpp * 7.5d19 * EFMp**(4.0d0) * (nP / N0)**(1.0d0 / 3.0d0) * &
               alphaPP * betaPP * nNu * T9**(8.0d0)
    end subroutine
    
    ! Neutrino emission due to Cooper pairing of neutrons in star core:
    ! YKGH(2001): (236), (241)
    ! Input: T - temperature
    !        TCN - neutron superfluidity temperature
    !        EFMn - effective neutron mass factor (mN*/mN)   
    !        pFN - neutron Fermi momentum
    ! Output: QCPCoreNN - neutrino emission rate in erg/(s cm^3)
    subroutine getQCPCoreNN(T, TCN, EFMn, pFN, QCPCoreNN)
        implicit none
        real(8) :: T, TCN, EFMn, pFN
        real(8) :: QCPCoreNN
        
        real(8) :: nuB, FBnu
        
        if (T > TCN) then
            QCPCoreNN = 0.0d0
        else
            nuB = dimensionlessGapB(T, TCN)
            FBnu = (1.204d0 * nuB * nuB + 3.733d0 * (nuB)**(4.d0) + 0.3191d0 * nuB**(6.0d0)) / & 
                   (1.0d0 + 0.3511d0 * nuB * nuB) * &
                   (0.7591d0 + DSQRT(0.2409d0**(2.0d0) + 0.3145d0 * nuB * nuB)) * &
                   DEXP(0.4616d0 - DSQRT(4.0d0 * nuB * nuB + 0.4616d0**(2.0d0)))    
            QCPCoreNN = 1.17d21 * EFMn * pFN / (M_N * c) * (T / 1.0d9)**(7.0d0) * nNu * aNB * FBnu
            ! correction factor from PPP(2015) section 3.3
            QCPCoreNN = QCPCoreNN * 0.19d0;
        end if
        
    end subroutine
    
    ! Neutrino emission due to Cooper pairing of neutrons in star crust:
    ! YKGH(2001): (236), (241)
    ! Input: T - temperature
    !        TCN - neutron superfluidity temperature
    !        EFMn - effective neutron mass factor (mN*/mN)   
    !        pFN - neutron Fermi momentum
    !        xnuc - nuclear size parameter = R_{eff} / R_{WS}
    ! Output: QCPCrustNN - neutrino emission rate in erg/(s cm^3)
    subroutine getQCPCrustNN(T, TCN, EFMn, pFN, xnuc, QCPCrustNN)
        implicit none
        real(8) :: T, TCN, EFMn, pFN, xnuc
        real(8) :: QCPCrustNN
        
        real(8) :: nuA, FAnu 
        
        if (T > TCN) then
            QCPCrustNN = 0.0d0
        else
            nuA = dimensionlessGapA(T, TCN)
            FAnu = (0.602d0 * nuA * nuA + 0.5942d0 * nuA**(4.0d0) + 0.288d0 * nuA**(6.0d0)) * &
                   (0.5547d0 + DSQRT(0.4453**(2.0d0) + 0.0113 * nuA * nuA))**(0.5d0) * &
                   DEXP(2.245d0 - DSQRT(4.0d0 * nuA * nuA + 2.245**(2.0d0)))
            QCPCrustNN = 1.17d21 * EFMn * pFN / (M_N * c) * (T / 1.0d9)**(7.0d0) * nNu * aNA * FAnu * & 
                        (1.0d0 - xnuc**(3.0d0))
            ! correction factor from PPP(2015) section 3.3
            QCPCrustNN = QCPCrustNN * (pFN / (EFMn * M_N * c))**(2.0d0)
        end if
        
    end subroutine
    
    ! Neutrino emission due to Cooper pairing of protons in star core:
    ! YKGH(2001): (236), (238), (241)
    ! Input: T - temperature
    !        TCP - proton superfluidity temperature
    !        EFMp - effective proton mass factor (mN*/mN)   
    !        pFP - proton Fermi momentum
    ! Output: QCPCorePP - neutrino emission rate in erg/(s cm^3)
    subroutine getQCPCorePP(T, TCP, EFMp, pFP, QCPCorePP)
        implicit none
        real(8) :: T, TCP, EFMp, pFP
        real(8) :: QCPCorePP
        
        real(8) :: nuA, aPA, FAnu 
        
        if (T > TCP) then
            QCPCorePP = 0.0d0
        else
            nuA = dimensionlessGapA(T, TCP)
            ! YKGH(2001): (238)
            aPA = 0.0064d0 + 1.59d0 * (pFP / (EFMp * M_P * c))**(2.0d0) * (EFMp**(2.0d0) + 11.0d0 / 42.0d0)
            FAnu = (0.602d0 * nuA * nuA + 0.5942d0 * nuA**(4.0d0) + 0.288d0 * nuA**(6.0d0)) * &
                   (0.5547d0 + DSQRT(0.4453**(2.0d0) + 0.0113 * nuA * nuA))**(0.5d0) * &
                   DEXP(2.245d0 - DSQRT(4.0d0 * nuA * nuA + 2.245**(2.0d0)))
            QCPCorePP = 1.17d21 * EFMp * pFP / (M_P * c) * (T / 1.0d9)**(7.0d0) * nNu * aPA * FAnu
            ! correction factor from PPP(2015) section 3.3
            QCPCorePP = QCPCorePP * (pFP / (EFMp * M_P * c))**(2.0d0)
        end if
        
    end subroutine
    
    ! This function estimates dimensionless gap amplitude for ^1S_{0} superfluidity
    ! YKGH(2001): Eq.(188)
    ! Input: T - temperature in K
    !        Tc - critical temperature in K
    function dimensionlessGapA(T, Tc)
        implicit none
        real(8) :: dimensionlessGapA
        
        real(8) :: T, Tc
        
        real(8) :: tau
        
        tau = T / Tc
        dimensionlessGapA = DSQRT(1 - tau) * (1.456d0 - 0.157d0 / DSQRT(tau) + 1.764d0 / tau)        
    end function
    
    ! This function estimates dimensionless gap amplitude for ^3P_{2}(m_{J} = 0) superfluidity
    ! YKGH(2001): Eq.(188)
    ! Input: T - temperature in K
    !        Tc - critical temperature in K
    function dimensionlessGapB(T, Tc)
        implicit none
        real(8) :: dimensionlessGapB
        
        real(8) :: T, Tc
        
        real(8) :: tau
        
        tau = T / Tc
        dimensionlessGapB = DSQRT(1 - tau) * (0.7893d0 + 1.188d0 / tau)
    end function
    
    ! This function estimates dimensionless gap amplitude for ^3P_{2}(m_{J} = 2) superfluidity
    ! YKGH(2001): Eq.(188)
    ! Input: T - temperature in K
    !        Tc - critical temperature in K
    function dimensionlessGapC(T, Tc)
        implicit none
        real(8) :: dimensionlessGapC
        
        real(8) :: T, Tc
        
        real(8) :: tau
        
        tau = T / Tc
        dimensionlessGapC = (1 - tau**(4.0d0))**(0.5d0) / tau * &
                            (2.030d0 - 0.4903 * tau**(4.0d0) + 0.1727d0 * tau**(8.0d0))
    end function
    
    ! This subroutine estimates superfluidity factor R^{D} for direct Urka process
    ! R^{D}_{B} - for neutron superfluidity, R^{D}_{A} - for proton superfluidity,
    ! R^{D}_{BA} - for neutron and proton superfluidity
    ! YKGH(2001): Eq.(199), Eq.(205)
	!!
	! This subroutine estimates superfluidity factors R^{Mn}, R^{Mp}
	! for modified Urka process in neutron and proton branches, respectively
	! R_{pA}^{Mn}, R_{pA}^{Mp} - for proton superfluidity
	! YKGH(2001): Eq.(211), Eq.(212)
	! R_{nB}^{Mn}, R_{nB}^{Mp} - for neutron superfluidity
    ! YKGH(2001): Eq.(227), Eq.(215)
	! R_{BA}^{Mn}, R_{BA}^{Mp} - for neutron and proton superfluidity 
	! YKGH(2001): Eq.(226), Eq.(225)
	!! 
	! This subroutine estimates superfluidity factors R^{nn}, R^{np}, R^{pp}
	! for neutrino bremsstrahlung in Nucleon-nucleon collisions
	! R_{pA}^{np}, R_{pA}^{pp} - for proton superfluidity
	! YKGH(2001): Eq.(220), Eq.(221)
	! R_{nB}^{np}, R_{nB}^{nn} - for neutron superfluidity
	! YKGH(2001): Eq.(220), Eq.(228)
	! R_{nB}^{np}, R_{nB}^{nn}, R_{pp} - for proton and neutron superfluidity
	! YKGH(2001): Eq.(229), Eq.(228), Eq.(221)
	
	! T - temperature in K, 
	! Tcp - critical temperature for proton superfluidity
	! Tcn - critical temperature for neutron superfluidity
    subroutine superfluidReductionFactors(T, Tcp, Tcn, RD, RMn, RMp, Rnn, Rnp, Rpp)
        implicit none
        real(8) :: RD, RMn, RMp, Rnn, Rnp, Rpp ! output
        
        real(8) :: T, Tcp, Tcn !input
        
        real(8) :: nuP, nuN
		real(8) :: Rba, Rnb, Rb ! for Eq.(225), Eq.(226)
		real(8) :: Rpa, Ra      ! for Eq.(226)
        
        if (Tcp + 1.0d-10 < T) then
            ! neutron superfluidity
            if (Tcn > T + 1.0d-10) then
				nuN = dimensionlessGapB(T, Tcn)
                RD = nSFDurcaR(nuN)
				RMn = nSFMurcaNR(nuN)
				RMp = nSFMurcaPR(nuN)
				Rnn = nSFBremNNR(nuN)
				Rnp = nSFBremNPR(nuN)	
				Rpp = 1.0d0
            ! non superfluidity                            
            else
                RD = 1.0d0
				RMn = 1.0d0
				RMp = 1.0d0
            end if
        else 
            ! neutron and proton superfluidity
            if (Tcn > T + 1.0d-10) then
                nuP = dimensionlessGapA(T, Tcp)
                nuN = dimensionlessGapB(T, Tcn)
				RD = npSFDurcaR(nuP, nuN)
				! Let's estimate R^{Mp} Eq.(225)
				Rba = npSFDurcaR(nuN, 2.0d0 * nuP)
				Rb = nSFDurcaR(nuN)
				Rnb = nSFMurcaPR(nuN)
                
				if ((Rba < 1.0d-10).OR.(Rnb < 1.0d-10)) then
                    RMp = 0.0d0
                else
                    RMp = Rba * Rnb / Rb
                end if
                
				! Let's estimate R^{Mn} Eq.(226)
				Rba = npSFDurcaR(2.0d0 * nuN, nuP)
				Ra = pSFDurcaR(nuP)
				Rpa = pSFMurcaNR(nuP)
				!!!!
                !write(*,'(a)') 'Rba            Ra            Rpa'
                !write(*,'(e15.6, e15.6, e15.6)') Rba, Ra, Rpa
				!!!!
                if ((Rba < 1.0d-10).OR.(Rpa < 1.0d-10)) then
                    RMn = 0.0d0
                else
                    RMn = Rba * Rpa / Ra
                endif
                
				! Let's estimate R^{np} Eq.(229)
				Rba = npSFDurcaR(nuN, nuP)
				Ra = pSFDurcaR(nuP)
				Rpa = pSFBremNPR(nuP)
                
                if ((Rba < 1.0d-10).OR.(Rpa < 1.0d-10)) then
                    Rnp = 0.0d0
                else
                    Rnp = Rba * Rpa / Ra
                endif
				
                ! Let's estimate R^{nn}
				Rnn = nSFBremNNR(nuN)
                !!!!
                !write(*,'(a)') '       Rnn      '
                !write(*, '(e15.6)')    Rnn
                !!!!
				Rpp = pSFBremPPR(nuP)
				
            ! proton superfluidity
            else 
				nuP = dimensionlessGapA(T, Tcp)
                RD = pSFDurcaR(nuP)
				RMn = pSFMurcaNR(nuP) 
				RMp = pSFMurcaPR(nuP)
				Rnn = 1.0d0
				Rnp = pSFBremNPR(nuP)
				Rpp = pSFBremPPR(nuP)
            end if
        end if
		
    end subroutine
	
	function nSFDurcaR(nuN)
        implicit none
		real(8) :: nSFDurcaR, nuN
		nSFDurcaR = (0.2546 + DSQRT((0.7454d0)**(2.0d0) + (0.1284d0 * nuN)**(2.0d0)))**(5.0d0) * &
                    DEXP(2.701d0 - DSQRT(2.701**(2.0d0) + nuN**(2.0d0)))
	end function
	
	function npSFDurcaR(nuN, nuP)
        implicit none
        real(8) :: npSFDurcaR, nuN, nuP
        real(8) :: candidate1, candidate2
		
		if ((nuN**(2.0d0) + nuP**(2.0d0)) < 25.0d0) then
            npSFDurcaR = (1.0d4 - 2.839d0 * nuP**(4.0d0) - 5.022d0 * nuN**(4.0d0)) / &
                         (1.0d4 + 757.0d0 * nuP**(2.0d0) + 1494.0d0 * nuN**(2.0d0) + & 
                         211.1d0 * (nuP * nuN)**(2.0d0) + 0.4832d0 * (nuP * nuN)**(4.0d0))				
        else 
            candidate1 = (0.2546 + DSQRT((0.7454d0)**(2.0d0) + (0.1284d0 * nuN)**(2.0d0)))**(5.0d0) * &
                         DEXP(2.701d0 - DSQRT(2.701**(2.0d0) + nuN**(2.0d0)))
            candidate2 = (0.2312d0 + DSQRT(0.7688d0**(2.0d0) + (0.1438d0 * nuP)**(2.0d0)))**(5.5d0) * &
                         DEXP(3.427d0 - DSQRT(3.427d0**(2.0d0) + nuP**(2.0d0)))
            npSFDurcaR = MIN(candidate1, candidate2)
        end if
	end function
	
	function pSFDurcaR(nuP)
        implicit none
		real(8) :: pSFDurcaR, nuP
		pSFDurcaR = (0.2312d0 + DSQRT(0.7688d0**(2.0d0) + (0.1438d0 * nuP)**(2.0d0)))**(5.5d0) * &
                    DEXP(3.427d0 - DSQRT(3.427d0**(2.0d0) + nuP**(2.0d0)))
	end function
	
	function nSFMurcaPR(nuN)
        implicit none
		real(8) :: nSFMurcaPR, nuN
		real(8) :: a, b
		
		a = 0.1612d0 + DSQRT(0.8388d0**(2.0d0) + (0.1117d0 * nuN)**(2.0d0))
		b = 0.1612d0 + DSQRT(0.8388d0**(2.0d0) + (0.1274d0 * nuN)**(2.0d0))
		nSFMurcaPR = (a**(7.0d0) + b**(5.0d0)) / 2.0d0 * & 
					 DEXP(2.398d0 - DSQRT(2.398d0**(2.0d0) + nuN**(2.0d0)))
	end function
	
	function nSFMurcaNR(nuN)
        implicit none
		real(8) :: nSFMurcaNR, nuN
		nSFMurcaNR = (0.2414d0 + DSQRT(0.7586d0**(2.0d0) + (0.1318d0 * nuN)**(2.0d0)))**(7.0d0) * &
					 DEXP(5.339d0 - DSQRT(5.339d0**(2.0d0) + (2.0d0 * nuN)**(2.0d0)))
	end function
	
	function pSFMurcaPR(nuP)
        implicit none
		real(8) :: pSFMurcaPR, nuP
		pSFMurcaPR = (0.2414d0 + DSQRT(0.7586d0**(2.0d0) + (0.1318d0 * nuP)**(2.0d0)))**(7.0d0) * &
					 DEXP(5.339d0 - DSQRT(5.339d0**(2.0d0) + (2.0d0 * nuP)**(2.0d0)))
	end function
	
	function pSFMurcaNR(nuP)
        implicit none
		real(8) :: pSFMurcaNR, nuP
		real(8) :: a, b
		
		a = 0.1477d0 + DSQRT(0.8523**(2.0d0) + (0.1175d0 * nuP)**(2.0d0))
	    b = 0.1477d0 + DSQRT(0.8523**(2.0d0) + (0.1297d0 * nuP)**(2.0d0))
		pSFMurcaNR = (a**(7.5d0) + b**(5.5d0)) / 2.0d0 * & 
					 DEXP(3.4370d0 - DSQRT(3.470d0**(2.0d0) + nuP**(2.0d0))) 
	end function
	
	function pSFBremNPR(nuP)
        implicit none
		real(8) :: pSFBremNPR, nuP
		real(8) :: a, b
		
		a = 0.9982d0 + DSQRT(0.0018d0**(2.0d0) + (0.3815d0 * nuP)**(2.0d0))
		b = 0.3949d0 + DSQRT(0.6051d0**(2.0d0) + (0.2666 * nuP)**(2.0d0))
		pSFBremNPR = (a * DEXP(1.306d0 - DSQRT(1.306d0**(2.0d0) + nuP**(2.0d0))) + &
			   1.732d0 * b**(7.0d0) * DEXP(3.303d0 - DSQRT(3.303**(2.0d0) + (2 * nuP)**(2.0d0)))) / 2.732d0
	end function
	
	function pSFBremPPR(nuP)    
        implicit none
		real(8) :: pSFBremPPR, nuP
		REAL(8) :: c, d
		
		c = 0.1747d0 + DSQRT(0.8253d0**(2.0d0) + (0.007933d0 * nuP)**(2.0d0))
		d = 0.7333d0 + DSQRT(0.2667d0**(2.0d0) + (0.1678d0 * nuP)**(2.0d0))
		pSFBremPPR = (c**(2.0d0) * DEXP(4.228d0 - DSQRT(4.228d0**(2.0d0) + & 
		             (2.0d0 * nuP)**(2.0d0))) + &
					 d**(7.5d0) * DEXP(7.762d0 - DSQRT(7.762d0**(2.0d0) + & 
			        (3.0d0 * nuP)**(2.0d0)))) / 2.0d0
	end function

	function nSFBremNNR(nuN)
        implicit none
		real(8) :: nSFBremNNR, nuN
		nSFBremNNR = pSFBremPPR(nuN)
	end function
	
	!????????????????????????????????????
	function nSFBremNPR(nuN)
        implicit none
		real(8) :: nSFBremNPR, nuN
		nSFBremNPR = pSFBremNPR(nuN)
	end function
	
 end module