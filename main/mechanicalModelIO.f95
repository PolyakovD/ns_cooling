!This code should be linked with bskfit.f95 tovbskfit.f95 mechanicalModel.f95
!This module is interface for mechanicalModel
module mechanicalModelIO

use mechanicalModel

implicit none
    
    integer(4), private, parameter :: fin = 10, fout = 20
    
    integer, private :: isDeep
    integer, private :: KEOS
    real(8), private :: centralRho, epsMass
	
	private :: TOVINTWrapper
    
contains
    
    subroutine constructMMFromCL(mM)
        type(mModel) :: mM 
        
        call setMMParametersFromCL()
        call TOVINTWrapper(mM)
        call mMConstructor(mM)
        
        contains
        
        subroutine setMMParametersFromCL()
            KEOS = 10
            write(*, *) 'Please, set the parameters of the mechanical structure'
            
            write(*, *) 'KEOS (= 19, 20 or 21)'
            read(*, *) KEOS
            do while ((KEOS < 19) .OR. (KEOS > 21))
                write(*, *) 'KEOS can be 19, 20 or 21'
                read(*, *) KEOS
            end do
            
            write(*, *) 'is our model deep? (enter "1" and model will be estimated to 10^{10} g/cc'
            write(*, *) 'enter other integer and model will be estimated to 10^{8} g/cc'
            read(*, *) isDeep
            
            write(*, *) 'Central mass density in g/cc'
            read(*, *) centralRho
            do while (centralRho < 0.0d0)
                write(*, *) 'Central mass density must be positive'
                read(*, *) centralRho
            end do
            
            write(*, *) 'Relative accuracy in total stellar mass'
            read(*, *) epsMass
            do while ((epsMass <= 0.0d0) .OR. (epsMass >= 1.0d0))
                write(*, *) 'Relative accuracy must lie in interval (0, 1)'
                read(*, *) epsMass
            end do
            
            write(*, *) 'Mechanical structure parameters is set successfully.'
            call setMMInputParameters(mM, KEOS, isDeep, centralRho, epsMass) 
            
        end subroutine
       
    end subroutine   
    
    subroutine constructMMFromFile(mM, inputfile)
        type(mModel) :: mM 
        character(*) :: inputfile
        
        call setMMParametersFromFile()
        call TOVINTWrapper(mM)
        call mMConstructor(mM)
        
        contains
        
        subroutine setMMParametersFromFile()
            open(fin, file = inputfile)
            
            read(fin, *) KEOS, isDeep
            read(fin, *) centralRho, epsMass
                
            close(fin)
            
            if ((KEOS < 19) .OR. (KEOS > 21)) then
                stop 'KEOS can be only 19, 20 or 21; check input file'
            end if
            if (centralRho < 0.0d0) then
                stop 'Central mass density must be positive; check input file'
            end if
            if ((epsMass <= 0.0d0) .OR. (epsMass >= 1.0d0)) then
                stop 'Relative accuracy must lie in interval (0, 1); check your file'
            end if
            call setMMInputParameters(mM, KEOS, isDeep, centralRho, epsMass) 
            
        end subroutine
        
    end subroutine
    
    subroutine TOVINTWrapper(mM)
    
        type(mModel) :: mM 
        
        character(len = 50) :: mechanicalTemporal = 'mechanicalData\\mechanicalTemporal.txt'
        
        real(8) ::  radius, stellarMass, g14   
        real(8) :: RLG, PLG, CHIrho, Gamma ! For bskFit
        real(8) :: GR_R ! For g14
        
        real(8), allocatable, dimension(:) :: spatialKnots, rhoDistribution, &
                                              pressureDistribution, lagrangianMass
        integer(4) :: n ! n == spatialKnotsNumber
        integer(4) :: i
        
        real(8), parameter :: M_PI = 3.141592653589793d0 
        
        call TOVINT(centralRho, KEOS, isDeep, epsMass, radius, stellarMass, n, mechanicalTemporal)
        ! code below have taken from tovbskfit.f95
        GR_R = dsqrt(1.d0 - .295423 * stellarMass / radius) ! radius (inverse vol.) corr.
        g14 = 1.32757d0 * stellarMass / radius**2 / GR_R    ! g(r) in 10^{14} cm/s^2
        ! code above have taken from tovbskfit.f95
        
        radius = radius * 1.0d+6
        stellarMass = stellarMass * sunMass
        
        ! Let's set estimated integral parameters in our model mM
        call setIntParameters(mM, n, radius, stellarMass, g14)
        
        allocate(spatialKnots(n), rhoDistribution(n), pressureDistribution(n), lagrangianMass(n))
        open(fin, file = mechanicalTemporal) 
        
        do i = 1, n
            read(fin, *) spatialKnots(i), rhoDistribution(i), lagrangianMass(i)
            lagrangianMass(i) = sunMass * lagrangianMass(i)
            
            RLG = log10(rhoDistribution(i))
            call BSKfit(KEOS, RLG, PLG, CHIrho, Gamma)
            pressureDistribution(i) = 10.0d0**PLG  
        end do
        
        close(fin)
        
        !Let's correct lagrangianMass for 2 first knots
        lagrangianMass(1) = M_PI * spatialKnots(2)**(3.0d0) * & 
                               (rhoDistribution(1) * 5.0d0 / 48.0d0 + rhoDistribution(2) / 16.0d0)
        lagrangianMass(2) = (lagrangianMass(1) * (spatialKnots(3) - spatialKnots(2)) + &
                                lagrangianMass(3) * spatialKnots(2)) / (spatialKnots(3)) / 2.0d0 
        
        call setDistParameters(mM, spatialKnots, rhoDistribution, pressureDistribution, lagrangianMass)
        
    end subroutine
    
    subroutine printMm(mM, outputFile)
        type(mModel) :: mM
        character(*) :: outputFile
        
        character(len = 2000) :: description
        
        integer(4) :: i, spatialKnotsNumber
        logical :: isDefined
        real(8) :: g14, radius
        real(8), allocatable, dimension(:) :: spatialKnots, rhoDistribution, pressureDistribution, & ! output distributions
                                              lagrangianMass, phiDistribution, redshiftFactor, volumeFactor
        
        isDefined = getDefinedStatus(mM)
        call getMMStructure(mM, spatialKnotsNumber, spatialKnots, rhoDistribution, pressureDistribution, &
                            redshiftFactor, volumeFactor, KEOS, isDeep, radius, g14)
        lagrangianMass = getLagrangianMass(mM) 
        phiDistribution = getPhiDistribution(mM)
        
        if (isDefined) then
            open(fout, file = outputFile)
            description = getMMDescription(mM)
            
            write(fout, '(a)') description
            write(fout, '(a, a)') ' spatialKnot      rho            pressure         lagrangianMass      Phi', &
                                '                    redshiftFactor              volumeFactor'
            do i = 1, spatialKnotsNumber
                write(fout, '(e12.7, 4x, e12.7, 4x, e12.7, 4x e12.7, 5x, e20.10, 5x, e20.10, 5x, e20.10)') & 
                           spatialKnots(i), rhoDistribution(i), pressureDistribution(i), lagrangianMass(i), &
                           phiDistribution(i), redshiftFactor(i), volumeFactor(i)
            end do
            
            close(fout)
        else
            write(*,*) 'mechanicalModelIO: The mechanical model is not defined.'
        end if
        
    end subroutine
    
    function getMMDescription(mM)
        type(mModel) :: mM
        character(len = 2000) :: getMMDescription
        character(len = 11) :: modelDepth
        
        real(8) :: stellarMass, g14, radius
        integer(4) :: spatialKnotsNumber 
        
        stellarMass = getStellarMass(mM)
        g14 = getG14(mM)
        radius = getRadius(mM)
        spatialKnotsNumber = getSpatialKnotsNumber(mM)
        
        
        if (isDeep .EQ. 1) then
            modelDepth = '10^{10}g/cc'
        else
            modelDepth = '10^{8}g/cc'
        end if
        
        write(getMMDescription, '(a, i3, a, a, a, e12.7, a, e12.7, a, e12.7, a, e12.3, a, e12.7, a, e12.7, a, i6)') & 
            'model: ', KEOS, &
            '\nmodel depth: ', modelDepth, &
            '\nrelative accuracy in total stellar mass SMASS ', epsMass, &
            '\ncentral grav.mass density [g/cc] ', centralRho, &
            '\nstellar mass [g] ', stellarMass, ' stellr mass [M_{sun}] ', stellarMass / sunMass, &
            '\ng14 - gravitational acceleration in 10^{14} cm/s^2 ', g14, &
            '\nradius ', radius, &
            '\nspatialKnotsNumber: ', spatialKnotsNumber   
    end function
    
end module 