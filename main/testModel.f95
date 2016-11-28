module testModel

use sphericalHeatEquationIO
use sphericalHeatEquation

implicit none

    type(heatEquation), private :: hE
    
    character(len = 50), private :: fileForHeatEquation

	integer(4), private :: spatialKnotsNumber
        
    real(8), private, allocatable, dimension(:) :: spatialKnots
    real(8), private :: stellarRadius = 0.5d0
    real(8), private, parameter :: STEFAN_BOLTZMANN = 5.67d-5
	
contains 

	subroutine conditionGenerator01(configFile)
        character(*) :: configFile
    
        integer, parameter :: n = 41, Nm = 2, fout = 101, foutSection = 107, fin = 100, m = 2000, Nout = 4
        real(8), parameter :: x0 = 0.0d0, t0 = 0.0d0, sigma = 0.5d0
        character(len = 50) :: sectionsParameters = 'inputData\\sectionsParameters.txt'
    
        integer :: i
        real(8) :: h, x
        
        open(fin, file = configFile)
        read(fin, *) fileForHeatEquation
        close(fin)
        
        open(fout, file = fileForHeatEquation)
    
        write(fout, *) sigma
        write(fout, *) n, t0
		allocate(spatialKnots(n))
		spatialKnotsNumber = n
        h = stellarRadius / real(n - 1)
        do i = 0, n - 1
            x = x0 + h * real(i)
			spatialKnots(i + 1) = x
            write(fout, '(f12.7, f12.7)') x, 1.0d0 - x * x
        end do
        
        close(fout)
        
        open(foutSection, file = sectionsParameters)
        
        write(foutSection, *) Nm
        write(foutSection, '(f12.7, f12.7, f12.7, i8)') 0.0d0, 0.0d0, 1.0d0, 2000
        write(foutSection, '(f12.7, f12.7, f12.7, i8)') 0.0d0, 0.0d0, 1.5d0, 2000
        write(foutSection, *) Nout
        write(foutSection, *) 1
        write(foutSection, *) 1001
        write(foutSection, *) 2001
        write(foutSection, *) 4001
        
        close(foutSection)
        
    end subroutine
    
    subroutine makeSolution()
        ! for Windows
        character(len = 50) :: heatOutputfile = 'outputData\\heatResult.txt'
        character(len = 50) :: outputfileForFlux = 'outputData\\resultFlux.txt'
        ! for Unix
        !character(len = 50) :: heatOutputfile = 'outputData//heatResult.txt'
        
        call constructSHEFromFile(hE, fileForHeatEquation)
        call hESolution(hE, conductivity03, sources03, capacity03, beta03, getSurfaceTemperatureFromBeta02, & 
                        relativisticFactors03, heatOutputfile, outputfileForFlux)
        
    end subroutine

    function conductivity01(thermal, spatialKnotsNumber, timeKnot)
        integer(4) i
		integer(4) :: spatialKnotsNumber
		real(8) :: conductivity01(2 : spatialKnotsNumber), thermal(spatialKnotsNumber)
		real(8) :: timeKnot
		do i = 2, spatialKnotsNumber
			conductivity01(i) = 1.0d0
		end do
	end function
	
	function sources01(thermal, spatialKnotsNumber, timeKnot)
        integer(4) i
		integer(4) :: spatialKnotsNumber
        real(8) :: sources01(spatialKnotsNumber), thermal(spatialKnotsNumber)
        real(8) :: timeKnot
        real(8) :: expT
        expT = exp(-timeKnot)
        do i = 1, spatialKnotsNumber
			sources01(i) = (5.0d0 + spatialKnots(i)**(2.0d0) ) * expT 
        end do
	end function
	
	function capacity01(thermal, spatialKnotsNumber, timeKnot)
        integer(4) i
		integer(4) :: spatialKnotsNumber
		real(8) :: capacity01(spatialKnotsNumber), thermal(spatialKnotsNumber)
		real(8) :: timeKnot
        do i = 1, spatialKnotsNumber
			capacity01(i) = 1.0d0
		end do
	end function
	
	function betaForTest(thermal, timeKnot)
        real(8) :: betaForTest, thermal, timeKnot
        betaForTest = (-1.0d0) / 3.0d0
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
    
	subroutine relativisticFactors01 (spatialKnotsNumber, redshiftFactor, volumeFactor)
		integer(4) i
		integer(4) :: spatialKnotsNumber
		real(8) :: redshiftFactor(1 : spatialKnotsNumber), volumeFactor(1 : spatialKnotsNumber)
		do i = 1, spatialKnotsNumber
			redshiftFactor(i) = 1.0d0
			volumeFactor(i) = 1.0d0
		end do
	end subroutine relativisticFactors01
	
	function conductivity03(thermal, spatialKnotsNumber, timeKnot)
        integer(4) i
		integer(4) :: spatialKnotsNumber
		real(8) :: conductivity03(2 : spatialKnotsNumber), thermal(spatialKnotsNumber)
		real(8) :: timeKnot
		do i = 2, spatialKnotsNumber
			conductivity03(i) = 0.75d0
		end do
	end function
	
	function sources03(thermal, spatialKnotsNumber, timeKnot)
        integer(4) i
		integer(4) :: spatialKnotsNumber
        real(8) :: sources03(spatialKnotsNumber), thermal(spatialKnotsNumber)
        real(8) :: timeKnot
        real(8) :: expT
        expT = exp(-timeKnot)
        do i = 1, spatialKnotsNumber
			sources03(i) = (4.0d0 + spatialKnots(i) * (10.5d0 + spatialKnots(i) * 6.5d0) ) * expT 
        end do
	end function
	
	function capacity03(thermal, spatialKnotsNumber, timeKnot)
        integer(4) i
		integer(4) :: spatialKnotsNumber
		real(8) :: capacity03(spatialKnotsNumber), thermal(spatialKnotsNumber)
		real(8) :: timeKnot
        do i = 1, spatialKnotsNumber
			capacity03(i) = 0.5d0
		end do
	end function
	
	function beta03(thermal, timeKnot)
        real(8) :: beta03, thermal, timeKnot
        beta03 = (-3.0d0) / 8.0d0
    end function
	
	subroutine relativisticFactors03 (spatialKnotsNumber, redshiftFactor, volumeFactor)
		integer(4) i
		integer(4) :: spatialKnotsNumber
		real(8) :: redshiftFactor(1 : spatialKnotsNumber), volumeFactor(1 : spatialKnotsNumber)
		do i = 1, spatialKnotsNumber
			redshiftFactor(i) = 1.0d0
			volumeFactor(i) = 1.0d0 / (1.0d0 + spatialKnots(i))
		end do
	end subroutine relativisticFactors03
	
end module