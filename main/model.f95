! Needs to be linked with: mechanicalModel.f95 mechanicalModelIO.f95 cvfull.f condBSk21.f neutrinos.f95
!                          compos21.f sfgaps.f
module model

use coreConductivity
use mechanicalModelIO
use mechanicalModel
use neutrinos
use sphericalHeatEquationIO
use sphericalHeatEquation

implicit none

    integer(4), private :: spatialKnotsNumber
       
	integer(4), private :: KEOS 
    real(8), private :: stellarRadius ! in cm
    real(8), private :: g14 
	 
    real(8), private, allocatable, dimension(:) :: spatialKnots
    real(8), private, allocatable, dimension(:) :: rhoDistribution
    real(8), private, allocatable, dimension(:) :: pressureDistribution
    real(8), private, allocatable, dimension(:) :: redshiftFactor
    real(8), private, allocatable, dimension(:) :: volumeFactor
    
    integer(4), private :: flashesNumber ! number of flashes
    integer(4), private :: flashesCounter ! ordinal number of the next flash
    
    real(8), private, allocatable, dimension(:) :: flashesEnergy 
    real(8), private, allocatable, dimension(:)	:: flashesEmissivity
    real(8), private, allocatable, dimension(:) :: flashesStartTime ! flashes start time from the star birth
    real(8), private, allocatable, dimension(:) :: flashesFinishTime ! flashes finish time from the star birth
    real(8), private, allocatable, dimension(:) :: flashesUpperRho
    integer(4), private, allocatable, dimension(:) :: flashesUpperKnots
    real(8), private, allocatable, dimension(:) :: flashesBottomRho
    integer(4), private, allocatable, dimension(:) :: flashesBottomKnots
	
    real(8), private, allocatable, dimension(:) :: neutron3P2Tc
    real(8), private, allocatable, dimension(:) :: neutron1S0Tc
    real(8), private, allocatable, dimension(:) :: proton1S0Tc
	
    real(8), public, parameter :: STEFAN_BOLTZMANN = 5.67d-5
    real(8), public, parameter :: N0 = 0.16d0 ! [fm -3] is the nucleon number density in standard nuclear matter
    real(8), public, parameter :: RHO_0 = 2.8d+14
    real(8), public, parameter :: PLANCK_BAR = 1.054571726d-27
    real(8), public, parameter :: M_E = 9.10958215d-28  !the rest electron mass 
    real(8), public, parameter :: M_P = 1.6726231d-24   !the rest proton mass 
    real(8), public, parameter :: M_N = 1.6749286d-24   !the rest neutron mass
    real(8), public, parameter :: K_BOLTAZMANN = 1.3806488d-16
    real(8), public, parameter :: XNC = .0809 ! number density of baryons in fm^{-3} at the crust/core 
                                              ! interface (Pearson et al. 2012)

contains 

	! Needs to be linked with: compos21.f bskfit.f condBSk21.f
	! This function compute the conductivity produced by electron - proton scattering in NS core
    function conductivity02(thermal, spatialKnotsNumber, timeKnot)
        integer(4) :: i, coreCrustBorder;
        integer(4) :: spatialKnotsNumber
        real(8) :: conductivity02(2 : spatialKnotsNumber), thermal(spatialKnotsNumber)
        real(8) :: timeKnot
		
        !for core conductivity
        real(8) rho, XN !for BsKNofR
		
        real(8) T
        real(8) Ye, Ymu !for FRACORE
        real(8) Yp
        real(8) nE, nMu, nP ! electron, muon and proton number density in [fm -3]
		
        real(8) EFMp, EFMn ! for EFFMASS
		
        !for crust conductivity
        real(8) Tlg, RHOlg
        real(8) SIGMA, SIGMAT, SIGMAH
        real(8) CKAPPA, CKAPPAT, CKAPPAH
		
        coreCrustBorder = 1
        do while (rhoDistribution(coreCrustBorder) > RHO_0 / 2.025d0) 
            coreCrustBorder = coreCrustBorder + 1
        end do
		
        ! Let's estimate conductivity in core
        do i = 2, (coreCrustBorder - 1) 
            rho = (rhoDistribution(i) + rhoDistribution(i - 1)) / 2.0d0
            T = (thermal(i) + thermal(i - 1)) / 2.0d0
			
            call BSkNofR(21, rho, XN)
            call FRACORE(21, XN, Ye, Ymu)
            Yp = Ye + Ymu
            nE = XN * Ye
            nMu = XN * Ymu
            nP = XN * Yp
            call EFFMASS(21, XN, Yp, EFMp, EFMn)
			
            call eMuConductivity(T, proton1S0Tc(i), nE, nMu, nP, EFMp, conductivity02(i))
            
        end do
		
        ! Let's estimate conductivity in crust
        do i = coreCrustBorder, spatialKnotsNumber	
            rho = (rhoDistribution(i) + rhoDistribution(i - 1)) / 2.0d0

            RHOlg = log10(rho) 
            Tlg = log10((thermal(i) + thermal(i - 1)) / 2.0d0)
            !write(*, *) i
            !write(*, *) Tlg, RHOlg
            
            call CONDBSk(Tlg, RHOlg, 0.0d0, 0.0d0, SIGMA, SIGMAT, SIGMAH, CKAPPA, CKAPPAT, CKAPPAH)
            conductivity02(i) = CKAPPA
        end do
    end function
    
    function sources02(thermal, spatialKnotsNumber, timeKnot)
        integer(4) :: spatialKnotsNumber
        real(8) :: sources02(spatialKnotsNumber), thermal(spatialKnotsNumber)
        real(8) :: timeKnot
        
        integer(4) :: i
        
        ! emissivity in the crust
        ! Neutrino emission due to e-e+ annihilation: YKGH(2001):Eq.(22)
        real(8) :: XN      ! baryon number density in fm^{-3}
        real(8) :: KZ, CMI, CMI1, xnuc, XND ! arguments for CRUBSK21 from compos.f
        real(8) :: Ye      ! number of electrons divided by number of nucleons
        real(8) :: nE      ! electron number density in fm^{-3}
        real(8) :: pFE     ! electron Fermi momentum 
        real(8) :: X, T    ! X - relativity parameter, T - temperature in relativistic units
        real(8) :: Qpair   ! neutrino emission rate in erg/(s cm^3)
        ! Neutrino emission due to plasmon decay: YKGH(2001):Eq.(38)
        real(8) :: Qpl     ! neutrino emission rate in erg/(s cm^3)
        ! Neutrino electron-bremsstrahlung: YKGH(2001):Eq.(76)
        real(8) :: T6      ! T6 = thermal(i) / 10^{6}  
        real(8) :: Qbr     ! neutrino emission rate in erg/(s cm^3)   
        ! Neutrino emission due to Cooper pairing of neutrons in star crust:
        ! YKGH(2001): (236), (241)
		real(8) :: QCPCrustNN ! neutrino emission rate in erg/(s cm^3)
        
        ! emissivity in the core
		! superfluidReductionFactors(T, Tcp, Tcn, RD, RMn, RMp, Rnn, Rnp, Rpp)
        ! threshold critical density for nucleon direct Urca process
        real(8), parameter :: RHO_CRIT1 = 8.2d14
		! reduction factors
		real(8) :: RD, RMn, RMp, Rnn, Rnp, Rpp
        ! Neutrino emission due to nucleon direct Urka process in nonesuperfluid case: YKGH(2001):Eq.(120)
        real(8) :: T9      ! T9 = thermal(i) / 10^{9}
        real(8) :: Ymu     ! number of muons divided by number of nucleons
        real(8) :: Yp      ! protons fraction
        real(8) :: EFMp, EFMn    ! effective proton and neutron mass factors(mP*/mP) ,(mN*/mN), respectively 
        real(8) :: QDUP     ! neutrino emission rate in erg/(s cm^3)
        ! Neutrino emission due to nucleon modified Urka process in nonesuperfluid case:
        ! YKGH(2001): Eq.(140), Eq.(142) with theta_{Mp} = 1
        real(8) :: nP, nN   ! proton and neutron number density in fm^{-3}
        real(8) :: pFP, pFN ! proton and neutron Fermi momentum, respectively
        real(8) :: QMUP     ! neutrino emission rate in erg/(s cm^3)
        ! Neutrino emission due to baryon-baryon collisions in nonesuperfluid case:
        ! YKGH(2001): Eq.(165), Eq.(166), Eq.(167)
        real(8) :: QBBR     ! neutrino emission rate in erg/(s cm^3)
        ! Neutrino emission due to Cooper pairing of neutrons in star core:
        ! YKGH(2001): Eq.(236), (241)
        real(8) :: QCPCoreNN
        ! Neutrino emission due to Cooper pairing of protons in star core:
        ! YKGH(2001): (236), (238), (241)
        real(8) :: QCPCorePP 
        
        do i = 1, spatialKnotsNumber
            sources02(i) = 0.0d0
        end do
         
        ! Let's estimate lockal heating 
        if (flashesCounter <= flashesNumber) then
            if ((timeKnot > flashesStartTime(flashesCounter)) .AND. (timeKnot < flashesFinishTime(flashesCounter))) then
                do i = flashesBottomKnots(flashesCounter), flashesUpperKnots(flashesCounter)
                    sources02(i) = flashesEmissivity(flashesCounter)
                end do
				write(*, *) 'heat!'
            end if
			if (timeKnot > flashesFinishTime(flashesCounter)) then
				flashesCounter = flashesCounter + 1
			end if
        end if
        
		! Let's estimate neutrino losses in crust
        i = spatialKnotsNumber
        do while(rhoDistribution(i) < (0.5d0 * RHO_0) )          
            T = thermal(i) / UNITEMP / 1.0d6
            !Let's find baryon number density from mass density 
            call BSkNofR(21, rhoDistribution(i), XN)
            call CRUBSk21(XN, KZ, CMI, CMI1, xnuc, XND)
            Ye = KZ / CMI1 
            nE = XN * Ye
            ! Let's find relativistic parameter X = pFE / (mE * c), 
            ! pFE = hbar * (3 * pi^{2} * nE)^{1 / 3}
            pFE = getPF(nE)
            X = pFE / (M_E * c)
            ! Let's find neutrino emission rate due to e-e+ annihilation
            call SINQpair(X, T, Qpair)
            !!!!
            if (Qpair < 0.0d0) then
                write(*, *) i, ' Qpair', Qpair
            end if
            !!!!
            sources02(i) = sources02(i) - Qpair
            !Let's find neutrino emission rate due to plasmon decay
            call SINQpl(X, T, Qpl)
            !!!!
            if (Qpl < 0.0d0) then
                write(*, *) i, ' Qpl'
            end if
            !!!!
            sources02(i) = sources02(i) - Qpl
            !Let's find neutrino emission rate due to electron-bremsstrahlung
            T6 = thermal(i) / 1.0d6
            call SINQbr(rhoDistribution(i), T6, Qbr)
            !!!!
            if (Qbr < 0.0d0) then
                write(*, *) i, ' Qbr'
            end if
            !!!!
            sources02(i) = sources02(i) - Qbr
            ! Let's find neutrino emission rate due to Cooper pairing of neutrons
            call EFFMASS(21, XN, 0.0d0, EFMp, EFMn)
            nN = (CMI1 - CMI) * XN / CMI1
            pFN = getPF(nN)
            call getQCPCrustNN(thermal(i), neutron1S0Tc(i), EFMn, pFN, xnuc, QCPCrustNN)
            !!!!
            if (QCPCrustNN < 0.0d0) then
                write(*, *) i, ' QCPCrustNN'
            end if
            !!!!
            if (nN > 1.0d-9) then
                sources02(i) = sources02(i) - QCPCrustNN
            end if
            
            i = i - 1
        end do
        
		!neutrino losses in core
        do while(i > 0)
			!Let's estimate reduction factors
			call superfluidReductionFactors(thermal(i), proton1S0Tc(i), neutron3P2Tc(i), &
				                            RD, RMn, RMp, Rnn, Rnp, Rpp)
            !Let's find neutrino emission rate due to nucleon direct Urka process
            call BSkNofR(21, rhoDistribution(i), XN)
            call FRACORE(21, XN, Ye, Ymu)
            Yp = Ye + Ymu
            nE = XN * Yp
            call EFFMASS(21, XN, Yp, EFMp, EFMn)
            T9 = thermal(i) / 1.0d9
                
            if (rhoDistribution(i) > RHO_CRIT1) then              
                call getQDUP(nE, T9, EFMp, EFMn, RD, QDUP)
                sources02(i) = sources02(i) - QDUP
            endif
                
            !Let's find neutrino emission rate due to nucleon modified Urca process
            nP = nE
            nN = XN - nP
            pFE = getPF(nE)
            pFP = getPF(nP)
            pFN = getPF(nN)
            call getQMUP(nP, T9, EFMp, EFMn, pFE, pFP, pFN, RMn, RMp, QMUP)
            
            sources02(i) = sources02(i) - QMUP
               
            ! Let's find neutrino emission rate due to baryon-baryon collisions
            call getQBBR(nN, nP, T9, EFMn, EFMp, Rnn, Rnp, Rpp, QBBR)
                
            sources02(i) = sources02(i) - QBBR
            
            ! Let's find neutrino emission rate due to Cooper pairing of neutrons
            call getQCPCoreNN(thermal(i), neutron3P2Tc(i), EFMn, pFN, QCPCoreNN)
            
            sources02(i) = sources02(i) - QCPCoreNN
            
            ! Let's find neutrino emission rate due to Cooper pairing of protons
            call getQCPCorePP(thermal(i), proton1S0Tc(i), EFMp, pFP, QCPCorePP)
            
            sources02(i) = sources02(i) - QCPCorePP
            
            i = i - 1
        end do
        
        contains
        
        !Input:  n -- number density of some particles in fm^{-3}
        !Output: the Fermi momentum in CGS pF = hbar * (3 * pi^{2} * n)^{1 / 3}
        function getPF(n) 
            real(8) :: getPF        
            real(8) :: n
            
            real(8), parameter :: M_PI = 3.141592653589793d0 
            
            getPF = PLANCK_BAR * (3 * M_PI**(2.0d0) * n * 1.0d39)**(1.0d0 / 3.0d0)
        end function

    end function
	
    ! Needs to be linked with: cvfull.f 
    function capacity02(thermal, spatialKnotsNumber, timeKnot)
        integer(4) :: spatialKnotsNumber
        real(8) :: capacity02(spatialKnotsNumber), thermal(spatialKnotsNumber)
        real(8) :: timeKnot
		
        real(8), parameter :: B12 = 0.0d0
        real(8) :: RHO, T6, XN, Zion, CMI, CMI1, Yn, Yp, Ye, Ymu, EFMp, EFMn, CVtot, CVE, CVI, CVC, CVN, CVmu
        real(8) :: xnuc, XNN, KZ, XND
        real(8) :: tauN, tauP ! T / Tcn, T / Tcp
        real(8) :: RCVN, RCVP ! neutron and proton superfluidity reduction factors
	
        integer(4) i
        
        do i = 1, spatialKnotsNumber
            T6 = thermal(i) / 1.0d+6
            RHO = rhoDistribution(i)
           !subroutine CVCRU21(XN,T6,B12,Zion,CMI,CMI1,EFMn,CVtot,CVE,CVI,CVC,CVN)
            call BSkNofR(KEOS, RHO, XN)
            
            if (XNC .lt. XN) then ! estimate capacity in core
                call FRACORE(KEOS,XN,Ye,Ymu)
                Yp = Ye + Ymu
                Yn = dim(1.d0,Yp)
                call EFFMASS(KEOS, XN, Yp, EFMp, EFMn)
                CMI = 1.008 * EFMp ! eff.mass of a proton in a.m.u.
                CMI1 = 1.d0 / Yp
                Zion = 1.d0 ! proton charge
                CVC = 0. ! neglect Coulomb nonideality

                call CVCOR21(XN, T6, B12, Ye, Ymu, EFMp, EFMn, CVtot, CVN, CVI, CVE, CVmu) ! CVI equal CVP in core

                tauN = thermal(i) / neutron3P2Tc(i)
                call REDUCV(tauN, 3, RCVN)
                tauP = thermal(i) / proton1S0Tc(i)
                call REDUCV(tauP, 2, RCVP)
                 
                CVtot = (CVE + CVmu + CVI * RCVP + CVN * RCVN)
                
            else ! estimate capacity in crust
                
                call CVFULL21(RHO,T6,B12,XN,Zion,CMI,CMI1,Yn,Yp,Ye,Ymu,EFMp,EFMn,CVtot,CVE,CVI,CVC,CVN)
            
                call CVCRU21(XN, T6, B12, Zion, CMI, CMI1, EFMn, CVtot, CVE, CVI, CVC, CVN)
                
                tauN = thermal(i) / neutron1S0Tc(i)
                call REDUCV(tauN, 1, RCVN)
                
                CVtot = CVtot - CVN
                CVtot = CVtot + CVN * RCVN
                
            end if
            capacity02(i) = CVtot * XN * 1.0d+39 * K_BOLTAZMANN ! cm^{3} == 10^{39} fm^{3}
        end do
		
    end function
    
    subroutine relativisticFactors02 (spatialKnotsNumber, currentRedshiftFactor, currentVolumeFactor)
		integer(4) i
		integer(4) :: spatialKnotsNumber
		real(8) :: currentRedshiftFactor(1 : spatialKnotsNumber), currentVolumeFactor(1 : spatialKnotsNumber)
        currentRedshiftFactor = redshiftFactor
        currentVolumeFactor = volumeFactor
	end subroutine relativisticFactors02
	
	function beta03(thermal, timeKnot)
		real(8) :: beta03, thermal, timeKnot
		
		real(8) :: FB, B12, COSLAT, TS6
		
		COSLAT = 0.0d0
		B12 = 0.0d0
		
		call TbTsFb(g14, thermal / 1.0d9, B12, COSLAT, TS6, FB)
		beta03 = -FB * stellarRadius**(2.0d0) / thermal
		
	end function
	
	! Input: beta
    !        T10g - temperature in K at the layer with rho = 10^{10}g/cc
    ! Output : surfaceTemperature in K
    !          surfaceFlux in erg/s
	subroutine getSurfaceTemperatureFromBeta03(beta, T10g, surfaceTemperature, surfaceFlux)
		real(8) :: beta, T10g
        real(8) :: surfaceTemperature, surfaceFlux
		
		real(8) B12, COSLAT, FB
		
		COSLAT = 0.0d0
		B12 = 0.0d0
		
		call TbTsFb(g14, T10g / 1.0d9, B12, COSLAT, surfaceTemperature, FB)
		surfaceTemperature = surfaceTemperature * 1.0d6
		surfaceFlux = -FB * 4.0d0 * M_PI * stellarRadius**(2.0d0) 
		
	end subroutine
	
	! Estimate beta for surface flux GPE(1983): Eq.(32)
    function beta02(thermal, timeKnot)
        real(8) :: beta02, thermal, timeKnot
		
        beta02 = -stellarRadius**(2.0d0) * STEFAN_BOLTZMANN * thermal**(1.1978d0)
        beta02 = beta02 * g14 * 10.0d0**(6.4176d0) / 1.7441
        
    end function
    
    ! Input: beta
    !        T10g - temperature in K at the layer with rho = 10^{10}g/cc
    ! Output : surfaceTemperature in K
    !          surfaceFlux in erg/s
    subroutine getSurfaceTemperatureFromBeta02(beta, T10g, surfaceTemperature, surfaceFlux)
        real(8) :: beta, T10g
        real(8) :: surfaceTemperature, surfaceFlux

        surfaceFlux = beta * T10g * 4.0d0 * M_PI
        surfaceTemperature = (-surfaceFlux / (4.0d0 * M_PI * stellarRadius**(2.0d0) * STEFAN_BOLTZMANN))**(0.25d0)
        
    end subroutine

    ! This function estimates the power of local uniform heating 
    ! in the layer of unstable nuclei is limited by the densities 
    ! rho1, rho2 (rho1 < rho2). Heat to last dt
    ! Input: totalEnergy -- total energy of the heating in [erg]
    !        rho1 -- density at the upper boundary of the layer [g]
    !        rho2 -- density at the bottom of the layer [g]
    !        dt -- duration of the heating [s]
    ! Output: QLockal -- the power of local uniform heating [erg/cc*s] 
    !         upperKnot -- spatial knot where upper limit of the layer passes   
    !         bottomKnot -- spatial knot where bottom limit of the layer passes
    subroutine estimateLockalHeating(totalEnergy, rho1, rho2, dt, QLockal, upperKnot, bottomKnot)
        real(8) :: estimateLockalHeatin
        
        real(8) :: totalEnergy, rho1, rho2, dt
        
        real(8) :: QLockal
        integer :: upperKnot, bottomKnot
        
        integer :: i
        real(8) :: r1, r2, volume
        
        i = spatialKnotsNumber
        do while(rhoDistribution(i) < rho1)
            i = i - 1
        end do
        r1 = rhoDistribution(i)
        upperKnot = i
        do while(rhoDistribution(i) < rho2)
            i = i - 1
        end do
        r2 = rhoDistribution(i + 1)
        bottomKnot = i + 1

        volume = integratingVolume(spatialKnotsNumber, upperKnot, bottomKnot, spatialKnots, redshiftFactor)
        QLockal = totalEnergy / volume / dt
        
    end subroutine
    
    subroutine setMMechanical(curSpatialKnotsNumber, curSpatialKnots, curRhoDistribution, curPressureDistribution, &
                          curRedshiftFactor, curVolumeFactor, curKEOS, curRadius, curG14)
        integer(4) :: curSpatialKnotsNumber
        real(8) :: curSpatialKnots(curSpatialKnotsNumber), curRhoDistribution(curSpatialKnotsNumber), &
                   curPressureDistribution(curSpatialKnotsNumber), curRedshiftFactor(curSpatialKnotsNumber), &
                   curVolumeFactor(curSpatialKnotsNumber)
        integer(4) :: curKEOS 
        real(8) :: curRadius, curG14
        integer(4) :: n
        
        n = curSpatialKnotsNumber
        allocate(spatialKnots(n), rhoDistribution(n), pressureDistribution(n), redshiftFactor(n), volumeFactor(n))
        
        spatialKnotsNumber = curSpatialKnotsNumber
        spatialKnots = curSpatialKnots
        rhoDistribution = curRhoDistribution
        pressureDistribution = curPressureDistribution
        redshiftFactor = curRedshiftFactor
        volumeFactor = curVolumeFactor
        KEOS = curKEOS
        stellarRadius = curRadius
        g14 = curG14
    end subroutine
    
    subroutine setFlashesNumber(curFlashesNumber)
        integer(4) :: curFlashesNumber
        flashesNumber = curFlashesNumber
        flashesCounter = 1
    end subroutine
    
    subroutine setFlashesParameters(curFlashesEnergy, curFlashesStartTime, curFlashesFinishTime, & 
                                    curFlashesUpperRho, curFlashesBottomRho, curFlashesEmissivity, & 
                                    curFlashesUpperKnots, curFlashesBottomKnots)
        real(8) :: curFlashesEnergy(flashesNumber), curFlashesStartTime(flashesNumber), & 
                   curFlashesFinishTime(flashesNumber), curFlashesUpperRho(flashesNumber), & 
                   curFlashesBottomRho(flashesNumber), curFlashesEmissivity(flashesNumber)
        integer(4) :: curFlashesUpperKnots(flashesNumber), curFlashesBottomKnots(flashesNumber)
        
        allocate(flashesEnergy(flashesNumber), flashesStartTime(flashesNumber), flashesFinishTime(flashesNumber), &
                 flashesUpperRho(flashesNumber), flashesBottomRho(flashesNumber), flashesEmissivity(flashesNumber), &
                 flashesUpperKnots(flashesNumber), flashesBottomKnots(flashesNumber))
        flashesEnergy = curFlashesEnergy
        flashesStartTime = curFlashesStartTime
        flashesFinishTime = curFlashesFinishTime
        flashesUpperRho = curFlashesUpperRho
        flashesBottomRho = curFlashesBottomRho
        flashesEmissivity = curFlashesEmissivity
        flashesUpperKnots = curFlashesUpperKnots
        flashesBottomKnots = curFlashesBottomKnots
    end subroutine
    
    subroutine setSFCurves(curNeutron3P2Tc, curNeutron1S0Tc, curProton1S0Tc)
        real(8) :: curNeutron3P2Tc(spatialKnotsNumber), curNeutron1S0Tc(spatialKnotsNumber), &
                   curProton1S0Tc(spatialKnotsNumber)
                   
        allocate(neutron3P2Tc(spatialKnotsNumber), neutron1S0Tc(spatialKnotsNumber), proton1S0Tc(spatialKnotsNumber))
        neutron3P2Tc = curNeutron3P2Tc
        neutron1S0Tc = curNeutron1S0Tc
        proton1S0Tc = curProton1S0Tc
    end subroutine
    
end module