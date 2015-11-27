module coreConductivity

implicit none

    real(8), private, parameter :: N0 = 0.16d0 ! [fm -3] is the nucleon number density in standard nuclear matter
    real(8), private, parameter :: M_P = 1.6726231d-24 ! [g] the rest proton mass 
    real(8), private, parameter :: M_N = 1.674927471d-24 ! [g] the rest neutron mass 
    real(8), private, parameter :: N_MU_FERMI = 0.0052d0 ! [fm -3] the muon Fermi number density

contains
    
    ! This subroutine estimates core eMu conductivity in npeMu-matter
    ! GY(1995): Eq.(65), Eq.(66)
    ! Input: T - temperature in K
    !        Tcp - proton critical superfluidity temperature in K
    !        nE, nMu, nP - electron muon and proton number density in [fm -3]
    !        EFMp - effective proton mass factors(mP*/mP)
    ! Output: eMuCond - thermal eMu conductivity in erg / (cm * s * K) 
    subroutine eMuConductivity(T, Tcp, nE, nMu, nP, EFMp, eMuCond)
        real(8) :: eMuCond
        
        real(8) :: T, Tcp, nE, nMu, nP, EFMp
        
        real(8) :: RP, ZP
        real(8) :: screeningMomentum
        real(8) :: nuE, nuMu
        real(8) :: crossNuEMu,   crossNuMuE
        real(8) :: tauE, tauMu
        
        real(8) eCond, muCond
        
        if (nMu > 1.0d-5) then
            call SFFactors(T, Tcp, RP, ZP)
		
            screeningMomentum = getScreeningMomentum(nE, nMu, nP, EFMp, ZP)
        
            call getEMuFrequency(screeningMomentum, nE, nMu, EFMp, T, RP, nuE, nuMu)
        
            call getCrossFrequency(screeningMomentum, nE, nMu, T, crossNuEMu, crossNuMuE)
        
            call getRelaxationTimes(nuE, nuMu, crossNuEMu, crossNuMuE, tauE, tauMu)
        
            eCond = 1.7d39 * 1.2d0 * (T / 1.0d8) * tauE * (nE / N0)**(2.0d0 / 3.0d0)
		
            muCond = eCond * nMu * tauMu / (nE * tauE)
		
            eMuCond = eCond + muCond
        else 
            if (T > Tcp) then ! Consider non-superfluid protons
                eMuCond = 6.3d32 * 1.2d0 / T * EFMp**(-0.5d0) * (nE / N0)**(7.0d0 / 6.0d0)
            else               ! Consider superfluid protons
                eMuCond = 4.3d32 * 1.2d0 / T * (nE / N0)
            end if
        end if       
    end subroutine
    
    ! This subroutine estimates superfluidity reduction factors R and Z for ^{1}S_{0} baryons 
    ! GY(1995): Eq.(45), Eq.(48)
    ! Input: T - temperature in K
    !        Tc - critical temperature in K
    ! Output: RSFFactor, ZSFFactor - superfluidity reduction factors for ^{1}S_{0} baryons 
    subroutine SFFactors(T, Tc, RSFFactor, ZSFFactor) 
        real(8) :: RSFFactor, ZSFFactor  

        real(8) :: T, Tc
        
        real(8) :: nuA
        
        if ((T + 1.0d0) > Tc) then
            RSFFactor = 1.0d0
            ZSFFactor = 1.0d0
        else
            nuA = dimensionlessGap(T, Tc)

            RSFFactor = (0.7694d0 + DSQRT(0.2306d0**(2.0d0) + (0.07207d0 * nuA)**(2.0d0)) + &
                        (27.0d0 * nuA * nuA + 0.1476d0 * nuA**(4.0d0)) * DEXP(-DSQRT(4.273d0**(2.0d0) + nuA * nuA)) + &
                        0.5051d0 * (DEXP(4.273d0 - DSQRT(4.273d0**(2.0d0) + nuA * nuA)) - 1.0d0) ) * &
                        DEXP(1.187d0 - DSQRT(1.187d0**(2.0d0) + nuA * nuA))
            
            ZSFFactor = (0.9443d0 + DSQRT(0.0557d0**(2.0d0) + (0.1886d0 * nuA)**(2.0d0)))**(0.5d0) * &
                        DEXP(1.753d0 - DSQRT(1.753d0**(2.0d0) + nuA * nuA))
        end if
        
        contains
		
        ! This internal function estimates dimensionless gap amplitude for ^1S_{0} superfluidity
        ! GY(1995): Eq.(39)
        ! Input: T - temperature in K
        !        Tc - critical temperature in K
        function dimensionlessGap(T, Tc)
            real(8) :: dimensionlessGap
        
            real(8) :: T, Tc
        
            real(8) :: tau
        
            tau = T / Tc
            dimensionlessGap = DSQRT(1 - tau) * (1.456d0 - 0.157d0 / DSQRT(tau) + 1.764d0 / tau)
        end function		
    end subroutine
    
    ! This function estimates screening momentum for protons
    ! GY(1995): Eq.(69)
    ! Input: nE, nMu, nP - electron muon and proton number density in [fm -3]
    !        EFMp - effective proton mass factors(mP*/mP)
    !        ZP - proton superfluidity screening reduction factor
    function getScreeningMomentum(nE, nMu, nP, EFMp, ZP)
        real(8) :: getScreeningMomentum
        
        real(8) :: nE, nMu, nP, EFMp, ZP
        
        getScreeningMomentum = 0.00929d0 * (1.0d0 + (nMu / nE)**(1.0d0/ 3.0d0) + &
                               2.83d0 * EFMp *(nP * N0 / nE**(2.0d0))**(1.0d0 / 3.0d0) * ZP)
    end function
    
    ! This subroutine estimates the frequency of collisions for electrons and muons in npemu matter
    ! GY(1995): Eq.(67), Eq.(68), Eq.(70), Eq.(71), Eq.(72), Eq.(73)
    ! Input: screeningMomentum - screening momentum for protons
    !        nE, nMu - electron and muon number density in [fm -3]
    !        EFMp - effective proton mass factors(mP*/mP)
    !        T - temperature in K
    !        RP - proton superfluidity reduction factor
    ! Output: nuE, nuMu - frequency of collisions for electrons and muons in npemu matter 
    subroutine getEMuFrequency(screeningMomentum, nE, nMu, EFMp, T, RP, nuE, nuMu)
        real(8) :: nuE, nuMu
        
        real(8) :: screeningMomentum, nE, nMu, EFMp, T, RP
        
        real(8) :: T8
        real(8) :: nuEP, nuMuP ! for Eq.(67), Eq.(68)
        real(8) :: nuEMu, nuMuE ! for Eq.(70), Eq.(71)
        real(8) :: nuEE, nuMuMu ! for Eq.(72), Eq.(73)
        
        T8 = T / 1.0d8
        
        nuEP = 1.15d12 * screeningMomentum**(-1.5d0) * EFMp**(2.0d0) * (N0 / nE) * T8 * T8 * RP
        nuMuP = nuEP * (nE / nMu)**(1.0d0 / 3.0d0)
        
        nuEMu = 1.43d11 * screeningMomentum**(-1.5d0) * (N0 / nE)**(1.0d0 / 3.0d0) * &
                (1.0d0 + 0.5d0 * (nMu / nE)**(2.0d0 / 3.0d0)) * T8 * T8
        nuMuE = nuEMu * (nE / nMu)**(1.0d0 / 3.0d0)
        
        nuEE = 3.58d0 * screeningMomentum**(-1.5d0) * (N0 / nE)**(1.0d0 / 3.0d0) * T8 * T8
        nuMuMu = nuEE * (nMu / nE) * (1.0d0 + 6.0d0 / 5.0d0 * (N_MU_FERMI / nMu)**(2.0d0 / 3.0d0) + &
                 2.0d0 / 5.0d0 * (N_MU_FERMI / nMu)**(1.0d0 / 3.0d0) )

        nuE = nuEP + nuEMu + nuEE
        nuMu = nuMuP + nuMuE + nuMuMu
    end subroutine
    
    ! This subroutine estimates the cross-frequencies of collisions for electrons and muons in npemu matter
    ! GY(1995): Eq.(74), Eq.(75)
    ! Input: screeningMomentum - screening momentum for protons
    !        nE, nMu - electron and muon number density in [fm -3]
    !        T - temperature in K
    ! Output: crossNuEMu, crossNuMuE - cross-frequencies of collisions for electrons and muons in npemu matter
    subroutine getCrossFrequency(screeningMomentum, nE, nMu, T, crossNuEMu, crossNuMuE)
        real(8) :: crossNuEMu, crossNuMuE
        
        real(8) :: screeningMomentum, nE, nMu, T
        
        crossNuEMu = 1.43d11 * screeningMomentum**(-1.5d0) * (N0 / nE)**(1.0d0 / 3.0d0) * &
                     (nMu / nE)**(2.0d0 / 3.0d0) * T * T / 1.0d16
        crossNuMuE = crossNuEMu * (nE / nMu)
    end subroutine
    
    ! This subroutine estimates relaxation time for electrons and muons
    ! GY(1995): Eq.(17)
    ! Input: nuE, nuMu - frequency of collisions for electrons and muons in npemu matter 
    !        crossNuEMu, crossNuMuE - cross-frequencies of collisions for electrons and muons in npemu matter
    ! Output: tauE, tauMu - relaxation time for electrons and muons, respectively 
    subroutine getRelaxationTimes(nuE, nuMu, crossNuEMu, crossNuMuE, tauE, tauMu)
        real(8) :: tauE, tauMu
        
        real(8) :: nuE, nuMu, crossNuEMu, crossNuMuE
        
        tauE = (nuMu - crossNuEMu) / (nuE * nuMu - crossNuEMu * crossNuMuE)
        tauMu = (nuE - crossNuMuE) / (nuE * nuMu - crossNuEMu * crossNuMuE)    
    end subroutine
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This subroutine estimates core n conductivity in npeMu-matter
    ! BHY(2001): Eq.(52)
    ! Input: EFMn, EFMp - effective neutron proton mass factors (mN*/mN), (mP*/mP)
    !        T - temperature in K
    !        kFn, kFp - [fm^-1] neutron and proton Fermi wave-vectors respectively
    !        Tcp, Tcn - neutron and proton critical superfluidity temperatures in K
    !        Nn - [fm^{-3}] neutron number density
    ! Output: nCond - [erg / (cm * s * K)] thermal n conductivity
    subroutine nConductivity(EFMn, EFMp, T, kFn, kFp, Tcn, Tcp, Nn, nCond)
        real(8) :: EFMn, EFMp, kFn, kFp, Tcn, Tcp, Nn
        
        real(8) :: nCond
        
        real(8) :: T8, RC, nuNN, nuNP
        
        T8 = T / 1.0d+8
        RC = RCFactor(T, Tcn)
        call getNPFrequency(EFMn, EFMp, T, kFn, kFp, Tcn, Tcp, nuNN, nuNP)
        
        nCond = 7.2d+38 * T8 * RC**2 * EFMn**(-1.0d0) / (nuNN + nuNP) * (nN / N0)   
    end subroutine
    
    ! This subroutine estimates the frequency of collisions for neutrons and proton in npemu matter
    ! BHY(2001): Eq.(53)
    ! Input: EFMn, EFMp - effective neutron proton mass factors (mN*/mN), (mP*/mP)
    !        T - temperature in K
    !        kFn, kFp - [fm^-1] neutron and proton Fermi wave-vectors respectively
    !        Tcp, Tcn - neutron and proton critical superfluidity temperatures in K
    ! Output: nuNN, nuNP - [s^{-1}] frequency of nn and np collisions respectively 
    subroutine getNPFrequency(EFMn, EFMp, T, kFn, kFp, Tcn, Tcp, nuNN, nuNP)
        real(8) :: EFMn, EFMp, kFn, kFp, Tcn, Tcp

        real(8) :: nuNN, nuNP
        
        real(8) :: T8, Rn1, Rn2, Rp1, Rp2, Kn1, Kn2, Kp1, Kp2, Sn1, Sn2, Sp1, Sp2
        real(8) :: tau, yn, yp
        
        T8 = T / 1.0d+8
        
        call RFactors(T, Tcn, Tcp, Rn1, Rn2, Rp1, Rp2)
        call KFactors(kFn, kFp, EFMn, EFMp, Kn1, Kn2, Kp1, Kp2)
        call S0Factors(kFn, kFp, Sn1, Sn2, Sp1, Sp2)
                
        nuNN = 3.48d+15 * EFMn**3 * T8**2 * &
               (Sn2 * Kn2 * Rn2 + 3 * Sn1 * Kn1 * (Rn1 - Rn2))
        nuNP = 3.48d+15 * EFMn * EFMp**2 * T8**2 * &
               (Sp2 * Kp2 * Rp2 + 0.5d0 * Kp1 * Sp1 * (3 * Rp1 - Rp2)) 
    end subroutine
    
    ! This function estimates SF reduction factor R_{c}(yn)
    ! BHY(2001): Eq.(51)
    ! Input: T - temperature in K
    !        Tcn - neutron critical superfluidity temperature in K
    function RCFactor(T, Tcn)
        real(8) :: T, Tcn
        
        real(8) :: RCFactor
        
        real(8) :: yn, tau
        if (T > Tcn) then 
            RCFactor = 1
        else
            tau = T / Tcn
            yn = DSQRT(1 - tau) * (0.7893d0 + 1.186d0 / tau) ! BHY(2001) : Eq.(39)
        
            RCFactor = (0.647d0 + DSQRT(0.353d0**2 + 0.109d0 * yn**2))**1.5d0 * &
                       DEXP(1.39d0 - DSQRT(1.39d0**2 + yn**2))
        endif
    end function
    
    ! This subroutine estimates SF factor R_{alpha}
    ! BHY(2001): Eq.(45), Eq.(46), Eq.(47), Eq.(48)
    ! Input: T - temperature in K
    !        Tcp, Tcn - neutron and proton critical superfluidity temperatures in K
    ! Output: Rn1, Rn2 - SF factors for nn collisions
    !         Rp1, Rp2 - SF factors for np collisions
    subroutine RFactors(T, Tcn, Tcp, Rn1, Rn2, Rp1, Rp2)
        real(8) :: T, Tcn, Tcp
        
        real(8) :: Rn1, Rn2, Rp1, Rp2
        
        real(8) :: yn, yp, tau
        
        if (T > Tcp) then
            if (T > Tcn) then ! without SF 
                Rn1 = 1
                Rn2 = 1
                Rp1 = 1
                Rp2 = 1
            else              ! neutron SF
                tau = T / Tcn
                yn = DSQRT(1 - tau) * (0.7893d0 + 1.186d0 / tau) ! BHY(2001) : Eq.(39)
                call nRFactors(yn, Rn1, Rn2, Rp1, Rp2)
            endif
        else
            if (T > Tcn) then ! proton SF 
                tau = T / Tcp
                yp = DSQRT(1 - tau) * (1.456d0 - 0.157d0 / DSQRT(tau) + 1.764d0 / tau) ! BHY(2001) : Eq.(38)
                call pRFactors(yp, Rn1, Rn2, Rp1, Rp2)
            else              ! neutron and proton SF
                tau = T / Tcn
                yn = DSQRT(1 - tau) * (0.7893d0 + 1.186d0 / tau) ! BHY(2001) : Eq.(39)
                tau = T / Tcp
                yp = DSQRT(1 - tau) * (1.456d0 - 0.157d0 / DSQRT(tau) + 1.764d0 / tau) ! BHY(2001) : Eq.(38)          
                call npRFactors(yn, yp, Rn1, Rn2, Rp1, Rp2)
            endif
        endif
        
        contains
        
        ! This internal subroutine estimates SF factors for neutron superfluidity in absence of 
        ! proton superfluidity
        ! BHY(2001): Eq.(45), Eq.(46)
        ! Input: yn - the neutron dimensionless gap
        ! Output: Rn1, Rn2 - SF factors for nn collisions
        !         Rp1, Rp2 - SF factors for np collisions
        subroutine nRFactors(yn, Rn1, Rn2, Rp1, Rp2)
            real(8) :: yn
            
            real(8) :: Rn1, Rn2, Rp1, Rp2
            
            Rn1 = 1.5d0 * (0.9468d0 + DSQRT(0.0532d0**2 + 0.5346d0 * yn**2))**3 * &
                  DEXP(0.377d0 - DSQRT(0.377d0**2 + 4 * yn**2)) + & 
                  (1 + 1.351d0 * yn**2)**2 / 3 * &
                  DEXP(0.169d0 - DSQRT(0.169d0**2 + 9 * yn**2))
            Rn2 = 0.5d0 * (0.6242d0 + DSQRT(0.3758d0**2 + 0.07198d0 * yn**2))**3 * &
                  DEXP(3.6724d0 - DSQRT(3.6724d0**2 + 4 * yn**2)) + &
                  0.5d0 * (1 + 0.01211d0 * yn**2)**9 * &
                  DEXP(7.5351d0 - DSQRT(7.5351d0**2 + 9 * yn**2))
            
            Rp1 = (0.4459d0 + DSQRT(0.5541d0**2 + 0.03016d0 * yn**2))**2 * &
                  DEXP(2.1178d0 - DSQRT(2.1178d0**2 + yn**2))
            Rp2 = (0.801d0 + DSQRT(0.199d0**2 + 0.04645d0 * yn**2))**2 * &
                  DEXP(2.3569d0 - DSQRT(2.3569d0**2 + yn**2))
        end subroutine
        
        ! This internal subroutine estimates SF factors for proton superfluidity in absence of 
        ! neutron superfluidity
        ! BHY(2001): Eq.(46)
        ! Input: yp - the proton dimensionless gap
        ! Output: Rn1, Rn2 - SF factors for nn collisions
        !         Rp1, Rp2 - SF factors for np collisions
        subroutine pRFactors(yp, Rn1, Rn2, Rp1, Rp2)
            real(8) :: yp
            
            real(8) :: Rn1, Rn2, Rp1, Rp2
            
            Rn1 = 1
            Rn2 = 1
            
            Rp1 = 0.5d0 * (0.3695d0 + DSQRT(0.6305d0**2 + 0.01064d0 * yp**2)) * &
                  DEXP(2.4451d0 - DSQRT(2.4451d0**2 + yp**2)) + &
                  0.5d0 * (1 + 0.1917 * yp**2)**(1.4d0) * &
                  DEXP(4.6627d0 - DSQRT(4.6627d0**2 + 4 * yp**2))
            Rp2 = 0.0436d0 * (DSQRT(4.345d0**2 + 19.55d0 * yp**2) - 3.345d0) * &
                  DEXP(2.0247 - DSQRT(2.0247d0**2 + yp**2)) + &
                  0.0654d0 * DEXP(8.992d0 - DSQRT(8.992d0**2 + 1.5d0 * yp**2)) + &
                  0.891d0 * DEXP(9.627d0 - DSQRT(9.627d0**2 + 9 * yp**2))
        end subroutine
        
        ! This internal subroutine estimates SF factors for proton and neutron superfluidity
        ! BHY(2001): Eq.(45), Eq.(47), Eq.(48)
        ! Input: yn, yp - the neutron and proton dimensionless gaps respectively
        ! Output: Rn1, Rn2 - SF factors for nn collisions
        !         Rp1, Rp2 - SF factors for np collisions
        subroutine npRFactors(yn, yp, Rn1, Rn2, Rp1, Rp2)
            real(8) :: yn, yp
            
            real(8) :: Rn1, Rn2, Rp1, Rp2
            
            real(8) :: yMin, yMax, un, up, uMin, uMax
            
            Rn1 = 1.5d0 * (0.9468d0 + DSQRT(0.0532d0**2 + 0.5346d0 * yn**2))**3 * &
                  DEXP(0.377d0 - DSQRT(0.377d0**2 + 4 * yn**2)) + & 
                  (1 + 1.351d0 * yn**2)**2 / 3 * &
                  DEXP(0.169d0 - DSQRT(0.169d0**2 + 9 * yn**2))
            Rn2 = 0.5d0 * (0.6242d0 + DSQRT(0.3758d0**2 + 0.07198d0 * yn**2))**3 * &
                  DEXP(3.6724d0 - DSQRT(3.6724d0**2 + 4 * yn**2)) + &
                  0.5d0 * (1 + 0.01211d0 * yn**2)**9 * &
                  DEXP(7.5351d0 - DSQRT(7.5351d0**2 + 9 * yn**2))
             
            yMin = MIN(yn, yp) 
            yMax = MAX(yn, yp)
            
            un = DSQRT(yn**2 + 1.485d0**2) - 1.485d0
            up = DSQRT(yp**2 + 1.485d0**2) - 1.485d0     
            uMin =  DSQRT(yMin**2 + 1.485d0**2) - 1.485d0
            uMax =  DSQRT(yMax**2 + 1.485d0**2) - 1.485d0
            
            Rp1 = DEXP(-uMax - uMin) * (0.7751d0 + 0.4823d0 * un + 0.1124d0 * up + &
                  0.04991d0 * un**2 + 0.08513d0 * un * up + 0.01284d0 * un**2 * up) + &
                  DEXP(-2 * uMax) * (0.2249d0 + 0.3539d0 * uMax - 0.2189d0 * uMin - &
                  0.6069d0 * un * uMin + 0.7362d0 * up * uMax)
                  
            un = DSQRT(yn**2 + 1.761d0**2) - 1.761d0
            up = DSQRT(yp**2 + 1.761d0**2) - 1.761d0
            uMin = DSQRT(yMin**2 + 1.761d0**2) - 1.761d0
            uMax = DSQRT(yMax**2 + 1.761d0**2) - 1.761d0
            
            Rp2 = DEXP(-uMax - uMin) * (1.1032d0 + 0.8645d0 * un + 0.2042d0 * up + &
                  0.07937d0 * un**2 + 0.1451d0 * un * up + 0.01333d0 * un**2 * up) + &
                  DEXP(-2 * uMax) * (-0.1032d0 - 0.234d0 * uMax + 0.06152d0 * un * uMax + &
                  0.7533d0 * un * uMin - 1.007d0 * up * uMax)
        end subroutine    
    end subroutine
    
    ! This subroutine estimates in-medium factors K for nn, np collisions
    ! BHY(2001): Eq.(30)
    ! Input: kFn, kFp - [fm^-1] neutron and proton Fermi wave-vectors respectively
    !        EFMn, EFMp - effective neutron proton mass factors (mN*/mN), (mP*/mP)
    ! Output: Kn1, Kn2 - in-medium factors for nn collisions
    !         Kp1, Kp2 - in-medium factors for np collisions
    subroutine KFactors(kFn, kFp, EFMn, EFMp, Kn1, Kn2, Kp1, Kp2)
        real(8) :: kFn, kFp, EFMn, EFMp
        
        real(8) :: Kn1, Kn2, Kp1, Kp2
        
        real(8) :: u
        
        u = kFn - 1.665d0 
        Kn1 = EFMn**(-2.0d0) * (0.4583d0 + 0.892d0 * u**2 - 0.5497d0 * u**3 - &
              0.06205d0 * kFp + 0.04022d0 * kFp**2 + 0.2122d0 * u * kFp)
              
        u = kFn - 1.556d0
        Kn2 = EFMn**(-2.0d0) * (0.4891d0 + 1.111d0 * u**2 - 0.2283d0 * u**3 + &
              0.01589d0 * kFp - 0.02099d0 * kFp**2 + 0.2773d0 * u * kFp)
              
        u = kFn - 2.126d0
        Kp1 = EFMp**(-2.0d0) * (0.04377d0 + 1.1d0 * u**2 + 0.118d0 * u**3 + &
              0.1626d0 * kFp + 0.3871d0 * u * kFp - 0.299d0 * u**4)
              
        u = kFn - 2.116d0 
        Kp2 = EFMp**(-2.0d0) * (0.0001313d0 + 1.248d0 * u**2 + 0.2403d0 * u**3 + &
              0.3257d0 * kFp + 0.5536d0 * u * kFp - 0.3237d0 * u**4 + &
              0.09786d0 * u**2 * kFp)
    end subroutine
    
    ! This subroutine estimates in-vacuum integrals 
    ! BHY(2001): Eq.(29)
    ! Input: kFn, kFp - [fm^-1] neutron and proton Fermi wave-vectors respectively
    ! Output: Sn1, Sn2 - in-vacuum integrals for nn collisions
    !         Sp1, Sp2 - in-vacuum integrals for np collisions
    subroutine S0Factors(kFn, kFp, Sn1, Sn2, Sp1, Sp2)
        real(8) :: kFn, kFp
        
        real(8) :: Sn1, Sn2, Sp1, Sp2

        Sn1 = 14.57d0 / kFn**1.5d0 * (1 - 0.0788d0 * kFn + 0.0883d0 * kFn**2) / &
              (1 - 0.1114d0 * kFn)
        Sn2 = 7.88d0 / kFn**2 * (1 - 0.2241d0 * kFn + 0.2006d0 * kFn**2) / &
              (1 - 0.1742d0 * kFn)
              
        Sp1 = 0.8007 * kFp / kFn**2 * (1 + 31.28d0 * kFp - 0.0004285d0 * kFp**2 + &
              26.85d0 * kFn + 0.08012d0 * kFn**2) / (1 - 0.5898d0 * kFn + &
              0.2368d0 * kFn**2 + 0.5838d0 * kFp**2 + 0.884d0 * kFn * kFp)
        Sp2 = 0.383d0 * kFp**4 / kFn**5.5d0 * (1 + 102 * kFp + 53.91d0 * kFn) / &
              (1 - 0.7087d0 * kFn + 0.2537d0 * kFn**2 + &
              9.404d0 * kFp**2 - 1.589d0 * kFn * kFp)
    end subroutine
    
end module 
