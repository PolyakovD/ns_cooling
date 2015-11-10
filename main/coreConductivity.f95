module coreConductivity

implicit none

    real(8), private, parameter :: N0 = 0.16d0 ! [fm -3] is the nucleon number density in standard nuclear matter
    real(8), private, parameter :: M_P = 1.6726231d-24   !the rest proton mass 
    real(8), private, parameter :: N_MU_FERMI = 0.0052d0 ! [fm -3] the muon Fermi number density

contains
    
    ! This subroutine estimates core conductivity in npeMu-matter
    ! GY(1995): Eq.(65), Eq.(66)
    ! Input: T - temperature in K
    !        Tcp - proton critical superfluidity temperature in K
    !        nE, nMu, nP - electron muon and proton number density in [fm -3]
    !        EFMp - effective proton mass factors(mP*/mP)
    ! Output: eMuCond - thermal conductivity in erg / (cm * s * K) 
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
    
end module 