!This module computes neutron star mechanical structure 
!The distributions are calculated from the center to the density rho_{b} = 10^{10} g/cc
module mechanicalModel

    implicit none
    
    type, public :: mModel
        
        logical, private :: isDefined = .false.
        
        integer(4), private :: KEOS !EOS number
        
        real(8), private :: epsMass !relative accuracy in total stellar mass 
        
        real(8), private :: centralRho
        real(8), private :: radius
        real(8), private :: stellarMass
        real(8), private :: g14
        
        integer(4), private :: spatialKnotsNumber
        
        real(8), private, allocatable, dimension(:) :: spatialKnots
        real(8), private, allocatable, dimension(:) :: rhoDistribution
        real(8), private, allocatable, dimension(:) :: pressureDistribution
        real(8), private, allocatable, dimension(:) :: lagrangianMass
        real(8), private, allocatable, dimension(:) :: phiDistribution 
        real(8), private, allocatable, dimension(:) :: redshiftFactor
        real(8), private, allocatable, dimension(:) :: volumeFactor
        
    end type mModel
    
    real(8), public, parameter :: sunMass = 1.989d+33
    real(8), public, parameter :: c = 3.0d+10
    real(8), public, parameter :: gG = 6.673d-8
    
contains
    
    subroutine mMConstructor(mM)
        type(mModel) :: mM
        integer(4) :: i, n
        
        n = mM%spatialKnotsNumber

        allocate(mM%phiDistribution(n), mM%redshiftFactor(n), mM%volumeFactor(n))
        call calculatePhi(150 * mM%spatialKnotsNumber)
        call calculateFactors()
       
        mM%isDefined = .true.
		
        contains
        
        !This internal subroutine solves differential equation for potential
        ! d Phi / d P = - 1 / (P + rho c^{2})
        ! with the boundary condition of the Schwarzschild
        ! Phi = ln (1 - r_{g} / R) / 2, 
        ! where r_{g} --  gravitational radius
        subroutine calculatePhi(gridSize)
            integer(4) :: gridSize
			
            real(8) :: uniformPressure(gridSize) ! pressure scale gluing from two uniform scales 
            real(8) :: uniformRho(gridSize) ! rho distribution by the pressure scale
            real(8) :: uniformPhi(gridSize) ! solution by the pressure scale
            real(8) :: rightPart(gridSize)  ! right part of equation by the pressure scale
			
            real(8) :: rG, dP
            
            integer(4) :: bigN

            call getUniformDistributions(uniformPressure, uniformRho, gridSize)
			
            rightPart = getRightPart(uniformPressure, uniformRho, gridSize)
			
            dP = uniformPressure(gridSize - 1) - uniformPressure(gridSize);
            rG = 2.0d0 * mM%stellarMass * gG / c**(2.0d0) 
			
            !Let's start calculation by the interpolation Adams method
            !lucky, we know right part, therefore it will be just
            
            !firs we will calculate surface layer
            uniformPhi(gridSize) = LOG(1.0d0 - rG / mM%radius) / 2.0d0
			
            uniformPhi(gridSize - 1) = uniformPhi(gridSize) + dP * (rightPart(gridSize) + rightPart(gridSize - 1)) / 2.0d0
            
            uniformPhi(gridSize - 2) = uniformPhi(gridSize - 1) + dP / 12.0d0 * &
                                       (5.0d0 * rightPart(gridSize - 2) + 8.0d0 * rightPart(gridSize - 1) - rightPart(gridSize))
                                      
            uniformPhi(gridSize - 3) = uniformPhi(gridSize - 2) + dP / 24.0d0 * &
                                       (9.0d0 * rightPart(gridSize - 3) + 19.0d0 * rightPart(gridSize - 2) - & 
                                        5.0d0 * rightPart(gridSize - 1) + rightPart(gridSize))
            do i = (gridSize - 4), (gridSize - 400), -1 
                 uniformPhi(i) = uniformPhi(i + 1) + dP / 720.0d0 * &
                                 (251.0d0 * rightPart(i) + 646.0d0 * rightPart(i + 1) - 264.0d0 * rightPart(i + 2) + &
                                  106.0d0 * rightPart(i + 3) - 19.0d0 * rightPart(i + 4))
            end do
            
            !then calculate remaining part of the star
            bigN = gridSize - 396
            dP = uniformPressure(1) - uniformPressure(2)
            
            uniformPhi(bigN - 5) = uniformPhi(bigN - 4) + dP / 720.0d0 * &
                                   (251.0d0 * rightPart(bigN - 5) + 646.0d0 * rightPart(bigN - 4) - & 
                                    264.0d0 * rightPart(bigN - 4 + 100) + 106.0d0 * rightPart(bigN - 4 + 200) - &
                                    19.0d0 * rightPart(bigN - 4 + 300))
            uniformPhi(bigN - 6) = uniformPhi(bigN - 5) + dP / 720.0d0 * &
                                   (251.0d0 * rightPart(bigN - 6) + 646.0d0 * rightPart(bigN - 5) - & 
                                    264.0d0 * rightPart(bigN - 4) + 106.0d0 * rightPart(bigN - 4 + 100) - &
                                    19.0d0 * rightPart(bigN - 4 + 200))
            uniformPhi(bigN - 7) = uniformPhi(bigN - 6) + dP / 720.0d0 * &
                                   (251.0d0 * rightPart(bigN - 7) + 646.0d0 * rightPart(bigN - 6) - & 
                                    264.0d0 * rightPart(bigN - 5) + 106.0d0 * rightPart(bigN - 4) - &
                                    19.0d0 * rightPart(bigN - 4 + 100))
           
            do i = (bigN - 8), 1, -1 
                 uniformPhi(i) = uniformPhi(i + 1) + dP / 720.0d0 * &
                                 (251.0d0 * rightPart(i) + 646.0d0 * rightPart(i + 1) - 264.0d0 * rightPart(i + 2) + &
                                  106.0d0 * rightPart(i + 3) - 19.0d0 * rightPart(i + 4))
            end do
            
            call getPhiFromUniform(uniformPressure, uniformPhi, gridSize)

        end subroutine
		
        !This internal subroutine makes rho distribution by the pressure scale to
        !apply the numerical method to solve the differential equation
        subroutine getUniformDistributions(uniformPressure, uniformRho, gridSize)
		    integer(4) :: gridSize
            real(8) :: uniformPressure(gridSize), uniformRho(gridSize)
			
            real(8) :: dP

            integer(4) :: j
            integer(4) :: bigN
            
            bigN = gridSize - 396
            dP = (mM%pressureDistribution(1) - mM%pressureDistribution(n)) / (bigN - 1)
            uniformPressure(1) = mM%pressureDistribution(1)
            do i = 2, bigN - 4
                uniformPressure(i) = mM%pressureDistribution(1) - dP * (i - 1)
            end do
            
            dP = dP / 100.0d0
            
            do i = 1, 400
                uniformPressure(bigN - 4 + i) = uniformPressure(bigN - 4) - dP * i 
            end do
            
            uniformRho(1) = mM%rhoDistribution(1)
            uniformRho(gridSize) = mM%rhoDistribution(n)
            j = 1
            do i = 2, gridSize - 1
                do while (mM%pressureDistribution(j + 1) > uniformPressure(i)) 
                    j = j + 1
                end do
                uniformRho(i) = mM%rhoDistribution(j) * (uniformPressure(i) - mM%pressureDistribution(j + 1)) + &
				                mM%rhoDistribution(j + 1) * (mM%pressureDistribution(j) - uniformPressure(i))
                uniformRho(i) = uniformRho(i) / (mM%pressureDistribution(j) - mM%pressureDistribution(j + 1))
            end do
            
        end subroutine
		
        !This internal function calculates right part of differential equation
        function getRightPart(uniformPressure, uniformRho, gridSize)
            integer(4) :: gridSize
            real(8) :: uniformPressure(gridSize), uniformRho(gridSize)
            
            real(8) :: getRightPart(gridSize)
            
            do i = 1, gridSize
                getRightPart(i) = -1.0d0 / (uniformPressure(i) + uniformRho(i) * c**(2.0d0)) 
                !!!!
                if (getRightPart(i) > 0.0d0) then
                    write(*,*) 'rightPart ', i, '  ', getRightPart(i)
                end if
                !!!!
            end do
            
		end function
        
        subroutine getPhiFromUniform(uniformPressure, uniformPhi, gridSize)
            integer(4) :: gridSize
            real(8) :: uniformPressure(gridSize), uniformPhi(gridSize)
            
            integer(4) :: j

            mM%phiDistribution(1) = uniformPhi(1)
            mM%phiDistribution(n) = uniformPhi(gridSize)

            j = 1
            do i = 2, n - 1
                do while (uniformPressure(j + 1) > mM%pressureDistribution(i))
                    j = j + 1
                end do
                mM%phiDistribution(i) = uniformPhi(j) * (mM%pressureDistribution(i) - uniformPressure(j + 1)) + &
                                        uniformPhi(j + 1) * (uniformPressure(j) - mM%pressureDistribution(i))
                mM%phiDistribution(i) = mM%phiDistribution(i) / (uniformPressure(j) - uniformPressure(j + 1))
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
	
    !Input: mM - the mechanical structure of the star
    !Output: spatialKnotsNumber - the number of grid nodes in mM
    !        spatialKnot - the grid nodes
    !        rhoDistribution - the distribution of mass in the star
    !        pressureDistribution - the distribution of pressure in the star
    !        redshiftFactor - distribution of redshift relativistic factor in the star
    !        volumeFactor - distribution of volume relativistic factor in the star
    !        myKEOS - EOS number 
    !        radius - stellar radius in cm 
    !        g14 - local gravitational acceleration in 10^{14} cm/s^2
    subroutine getMMStructure(mM, spatialKnotsNumber, spatialKnots, rhoDistribution, pressureDistribution, &
                            redshiftFactor, volumeFactor, myKEOS, radius, g14)
        type(mModel) :: mM ! input model
        
        integer(4) :: spatialKnotsNumber 
        real(8), allocatable, dimension(:) :: spatialKnots, rhoDistribution, pressureDistribution, & ! output distributions
                                              redshiftFactor, volumeFactor
        integer(4) :: myKEOS 
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
    
    subroutine setMMInputParameters(mM, KEOS, centralRho, epsMass)
        type(mModel) :: mM
        integer :: KEOS
        real(8) :: centralRho, epsMass
        
        mM%KEOS = KEOS
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