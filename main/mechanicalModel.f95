!This module computes neutron star mechanical structure 
!The distributions are calculated from the center to the density rho_{b} = 10^{10} g/cc
module mechanicalModel

    implicit none
    
    type, public :: mModel
        
        logical, private :: isDefined = .false.
        
        integer(4), private :: KEOS ! EOS number
        integer(4), private :: isDeep 
        
        real(8), private :: epsMass ! relative accuracy in total stellar mass 
        
        real(8), private :: centralRho
        real(8), private :: radius  ! radius on the surface in cm > spatialKnots(spatialKnotsNumber)
        real(8), private :: stellarMass
        real(8), private :: g14
        
        integer(4), private :: spatialKnotsNumber
        
        real(8), private, allocatable :: spatialKnots(:)
        real(8), private, allocatable :: rhoDistribution(:)
        real(8), private, allocatable :: pressureDistribution(:)
        real(8), private, allocatable :: lagrangianMass(:)
        real(8), private, allocatable :: phiDistribution(:) 
        real(8), private, allocatable :: redshiftFactor(:)
        real(8), private, allocatable :: volumeFactor(:)
        
    end type mModel
    
    real(8), public, parameter :: sunMass = 1.989d+33
    real(8), public, parameter :: c = 3.0d+10
    real(8), public, parameter :: gG = 6.673d-8
    
contains
    
    subroutine mMConstructor(mM)
        type(mModel) :: mM
        integer(4) :: i, j, n
        
        n = mM%spatialKnotsNumber

        allocate(mM%phiDistribution(n), mM%redshiftFactor(n), mM%volumeFactor(n))
        call calculatePhi02(500 * mM%spatialKnotsNumber)
        call calculateFactors()
       
        mM%isDefined = .true.
		
        contains
        
        ! This internal subroutine solves differential equation for potential
        ! d Phi / d P = - 1 / (P + rho c^{2})
        ! with the boundary condition of the Schwarzschild
        ! Phi = ln (1 - r_{g} / R) / 2, 
        ! where r_{g} --  gravitational radius
        ! this equation equal to d Phi / d t = exp(t) / (exp(t) + rho * c^{2}),
        ! where t = log p
        ! the equation is solved by the implicit method of Adams 
        ! (А.А. Самарский, А.В. Гулин Численные методы М: Наука 1989, стр. 235-236)
        subroutine calculatePhi02(gridSize)
            integer(4) :: gridSize
        
            real(8) :: uniformLogPressure(gridSize) ! log(pressure) distribution
            real(8) :: uniformRho(gridSize) ! rho distribution by the log(pressure) scale
            real(8) :: uniformPhi(gridSize) ! solution by the log(pressure) scale
            real(8) :: rightPart(gridSize) ! right part of equation by the log(pressure) scale
            
            real(8) rG, dLogP
            
            call getUniformDistributions02(uniformLogPressure, uniformRho, gridSize)
            
            rightPart = getRightPart02(uniformLogPressure, uniformRho, gridSize)
            
            dLogP = uniformLogPressure(gridSize - 1) - uniformLogPressure(gridSize);
            rG = 2.0d0 * mM%stellarMass * gG / c**2 
            
            !Let's start calculation by the interpolation Adams method
            !lucky, we know right part, therefore it will be just
            uniformPhi(gridSize) = DLOG(1.0d0 - rG / mM%radius) / 2.0d0
			
            uniformPhi(gridSize - 1) = uniformPhi(gridSize) + dLogP * (rightPart(gridSize) + rightPart(gridSize - 1)) / 2.0d0
            
            uniformPhi(gridSize - 2) = uniformPhi(gridSize - 1) + dLogP / 12.0d0 * &
                                       (5.0d0 * rightPart(gridSize - 2) + 8.0d0 * rightPart(gridSize - 1) - rightPart(gridSize))
                                      
            uniformPhi(gridSize - 3) = uniformPhi(gridSize - 2) + dLogP / 24.0d0 * &
                                       (9.0d0 * rightPart(gridSize - 3) + 19.0d0 * rightPart(gridSize - 2) - & 
                                        5.0d0 * rightPart(gridSize - 1) + rightPart(gridSize))
            do i = (gridSize - 4), 1, -1 
                 uniformPhi(i) = uniformPhi(i + 1) + dLogP / 720.0d0 * &
                                 (251.0d0 * rightPart(i) + 646.0d0 * rightPart(i + 1) - 264.0d0 * rightPart(i + 2) + &
                                  106.0d0 * rightPart(i + 3) - 19.0d0 * rightPart(i + 4))
            end do
            
            call getPhiFromUniform02(uniformLogPressure, uniformPhi, gridSize)
        
        end subroutine
		
        ! This internal subroutine makes rho distribution by the log(pressure) scale to
        ! apply the numerical method to solve the differential equation
        subroutine getUniformDistributions02(uniformLogPressure, uniformRho, gridSize)
            integer(4) :: gridSize
            real(8) :: uniformLogPressure(gridSize), uniformRho(gridSize)
            
            real(8) :: dLogP
            
            dLogP = (DLOG(mM%pressureDistribution(1)) - DLOG(mM%pressureDistribution(n))) / gridSize
            uniformLogPressure(1) = DLOG(mM%pressureDistribution(1))
            do i = 2, gridSize
                uniformLogPressure(i) = uniformLogPressure(1) - dLogP * (i - 1)
            end do

            uniformRho(1) = mM%rhoDistribution(1)
            uniformRho(gridSize) = mM%rhoDistribution(n)
            j = 1
            do i = 2, gridSize - 1
                do while (mM%pressureDistribution(j + 1) > DEXP(uniformLogPressure(i))) 
                    j = j + 1
                end do
                uniformRho(i) = mM%rhoDistribution(j) * (DEXP(uniformLogPressure(i)) - mM%pressureDistribution(j + 1)) + &
				                mM%rhoDistribution(j + 1) * (mM%pressureDistribution(j) - DEXP(uniformLogPressure(i)))
                uniformRho(i) = uniformRho(i) / (mM%pressureDistribution(j) - mM%pressureDistribution(j + 1))
            end do
            
        end subroutine
		
        function getRightPart02(uniformLogPressure, uniformRho, gridSize)
            integer(4) :: gridSize
            real(8) :: uniformLogPressure(gridSize), uniformRho(gridSize)
        
            real(8) :: getRightPart02(gridSize)
            
            do i = 1, gridSize
                getRightPart02(i) = -DEXP(uniformLogPressure(i)) / & 
                                  (DEXP(uniformLogPressure(i)) + uniformRho(i) * c**2) 
            end do
            
        end function
        
        subroutine getPhiFromUniform02(uniformLogPressure, uniformPhi, gridSize)
            integer(4) :: gridSize
            real(8) :: uniformLogPressure(gridSize), uniformPhi(gridSize)
            
            integer(4) :: j

            mM%phiDistribution(1) = uniformPhi(1)
            mM%phiDistribution(n) = uniformPhi(gridSize)

            j = 1
            do i = 2, n - 1
                do while (DEXP(uniformLogPressure(j + 1)) > mM%pressureDistribution(i))
                    j = j + 1
                end do
                mM%phiDistribution(i) = uniformPhi(j) * (mM%pressureDistribution(i) - DEXP(uniformLogPressure(j + 1))) + &
                                        uniformPhi(j + 1) * (DEXP(uniformLogPressure(j)) - mM%pressureDistribution(i))
                mM%phiDistribution(i) = mM%phiDistribution(i) / & 
                                        (DEXP(uniformLogPressure(j)) - DEXP(uniformLogPressure(j + 1)))
            end do
            
        end subroutine 
        
        !This subroutine calculates relativistic factors from phiDistribution
        subroutine calculateFactors()
            real(8) :: multiplier
            multiplier = 2.0 * gG / c**(2.0d0)
            mM%volumeFactor(1) = (1.0d0 - multiplier * mM%lagrangianMass(1) / mM%spatialKnots(2) / 2.0d0)**(-0.5d0)
            mM%redshiftFactor(1) = EXP(mM%phiDistribution(1))
            do i = 2, n
                mM%redshiftFactor(i) = EXP(mM%phiDistribution(i))
                mM%volumeFactor(i) = (1.0d0 - multiplier * mM%lagrangianMass(i) / mM%spatialKnots(i))**(-0.5d0)
            end do
         end subroutine
         
    end subroutine
	
    ! Input: mM - the mechanical structure of the star
    ! Output: spatialKnotsNumber - the number of grid nodes in mM
    !        spatialKnot - the grid nodes
    !        rhoDistribution - the distribution of mass in the star
    !        pressureDistribution - the distribution of pressure in the star
    !        redshiftFactor - distribution of redshift relativistic factor in the star
    !        volumeFactor - distribution of volume relativistic factor in the star
    !        myKEOS - EOS number 
    !        isDeep - model depth('1' equal to 10^{10}g/cc, other equal to 10^{8}g/cc)
    !        radius - stellar radius in cm 
    !        g14 - local gravitational acceleration in 10^{14} cm/s^2
    subroutine getMMStructure(mM, spatialKnotsNumber, spatialKnots, rhoDistribution, pressureDistribution, &
                            redshiftFactor, volumeFactor, myKEOS, isDeep, radius, g14)
        type(mModel) :: mM ! input model
        
        integer(4) :: spatialKnotsNumber 
        real(8), allocatable, dimension(:) :: spatialKnots, rhoDistribution, pressureDistribution, & ! output distributions
                                              redshiftFactor, volumeFactor
        integer(4) :: myKEOS, isDeep
        real(8) :: radius, g14
        
        integer(4) :: n 
		
        spatialKnotsNumber = mM%spatialKnotsNumber
        n = mM%spatialKnotsNumber
        allocate(spatialKnots(n), rhoDistribution(n), pressureDistribution(n), redshiftFactor(n), volumeFactor(n))
        
        spatialKnots = mM%spatialKnots
        rhoDistribution = mM%rhoDistribution
        pressureDistribution = mM%pressureDistribution
        redshiftFactor = mM%redshiftFactor
        volumeFactor = mM%volumeFactor
		myKEOS = mM%KEOS
        isDeep = mM%isDeep
        radius = mM%radius
        g14 = mM%g14
        
    end subroutine
    
    function getStellarMass(mM)
        real(8) :: getStellarMass
        type(mModel) :: mM
        getStellarMass = mM%stellarMass
    end function
    
    function getLagrangianMass(mM) 
        real(8), allocatable, dimension(:) :: getLagrangianMass
        type(mModel) :: mM
        allocate(getLagrangianMass(mM%spatialKnotsNumber))
        getLagrangianMass = mM%lagrangianMass
    end function
    
    function getPhiDistribution(mM) 
        real(8), allocatable, dimension(:) :: getPhiDistribution
        type(mModel) :: mM
        allocate(getPhiDistribution(mM%spatialKnotsNumber))
        getPhiDistribution = mM%phiDistribution
    end function
    
    function getDefinedStatus(mM) 
        logical :: getDefinedStatus
        type(mModel) :: mM
        getDefinedStatus = mM%isDefined
    end function
    
    function getG14(mM)
        real(8) :: getG14
        type(mModel) :: mM
        getG14 = mM%g14
    end function
        
    function getRadius(mM)    
        real(8) :: getRadius
        type(mModel) :: mM
        getRadius = mM%radius
    end function
    
    function getSpatialKnotsNumber(mM)    
        integer(4) :: getSpatialKnotsNumber
        type(mModel) :: mM
        getSpatialKnotsNumber = mM%spatialKnotsNumber
    end function
    
    subroutine validateModel(mM)
        type(mModel) :: mM
        
        integer(4) :: i
        
        if (mM%isDefined) then
            do i = 2, mM%spatialKnotsNumber
                if (mM%spatialKnots(i) < mM%spatialKnots(i - 1)) then
                    write(*, *) 'mechanicalModel: we have problems with the monotony of spatial nodes.'
                    write(*, *) 'number from stellar center'
                    write(*, '(i12, e18.7)') i, mM%spatialKnots(i)
                    write(*, '(i12, e18.7)') (i + 1), mM%spatialKnots(i + 1)
                end if
            end do  
        else
            write(*,*) 'The mechanical model is not defined.'
        end if
        
    end subroutine
    
    subroutine setMMInputParameters(mM, KEOS, isDeep, centralRho, epsMass)
        type(mModel) :: mM
        integer :: KEOS, isDeep
        real(8) :: centralRho, epsMass
        
        if (centralRho < 5.0d+13) then
            write(*, *) 'mechanicalModel:Warning: centralRho = ', centralRho
        end if
        if (epsMass > 5.0d-6) then
            write(*, *) 'mechanicalModel:Warning: epsMass = ', epsMass
        end if
        
        mM%KEOS = KEOS
        mM%isDeep = isDeep
        mM%centralRho = centralRho
        mM%epsMass = epsMass
        
    end subroutine
    
    subroutine setIntParameters(mM, spatialKnotsNumber, radius, stellarMass, g14)
        type(mModel) :: mM
        integer(4) :: spatialKnotsNumber 
        real(8) :: radius, stellarMass, g14
        
        mM%spatialKnotsNumber = spatialKnotsNumber
        mM%radius = radius
        mM%stellarMass = stellarMass
        mM%g14 = g14
    
    end subroutine
    
    ! Input: mM - our completed mechanicalModel
    ! Output: Qpl - neutrino emission rate in erg/(s cm^3)
    subroutine setDistParameters(mM, spatialKnots, rhoDistribution, pressureDistribution, lagrangianMass)
        type(mModel) :: mM
        real(8) :: rhoDistribution(mM%spatialKnotsNumber), pressureDistribution(mM%spatialKnotsNumber), &
                   lagrangianMass(mM%spatialKnotsNumber), spatialKnots(mM%spatialKnotsNumber)
        integer(4) :: n
        
        n = mM%spatialKnotsNumber
        
        allocate(mM%spatialKnots(n), mM%rhoDistribution(n), mM%pressureDistribution(n), mM%lagrangianMass(n))
        mM%spatialKnots = spatialKnots
        mM%rhoDistribution = rhoDistribution
        mM%pressureDistribution = pressureDistribution 
        mM%lagrangianMass = lagrangianMass
        
    end subroutine
    
end module 