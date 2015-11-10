!This module implements some simple calculations algorithms
module calculations 

implicit none
    
    real(8), public, parameter :: M_PI = 3.141592653589793d0    
    
contains
 
    ! estimate int_{bottomR}^{upperR} 4 PI r^{2} integrand dr
    ! upperKnot > bottomKnot
    function integratingVolume(spatialKnotsNumber, upperKnot, bottomKnot, spatialKnots, integrand)
        real(8) :: integratingVolume
        
        integer(4) :: spatialKnotsNumber, upperKnot, bottomKnot
        real(8) :: spatialKnots(1 : spatialKnotsNumber), integrand(1 : spatialKnotsNumber)
        
        integer(4) :: i
        real(8) :: upperR, bottomR
        
        if (upperKnot <= bottomKnot) then
            stop 'Calculations: upperKnot and bottomKnot must satisfy: upperKnot > bottomKnot'
        end if
        
        integratingVolume = 0.0d0
        if (bottomKnot /= 1) then        
            upperR = spatialKnots(bottomKnot)
            bottomR = (spatialKnots(bottomKnot) + spatialKnots(bottomKnot - 1)) / 2.0d0
            integratingVolume = simpleIntegral() * &
                (3.0d0 / 4.0d0 * integrand(bottomKnot) + &
                 1.0d0 / 4.0d0 * integrand(bottomKnot - 1))
        end if
        
        do i = bottomKnot, upperKnot - 1, 1
            upperR = spatialKnots(i + 1)
            bottomR = spatialKnots(i)
            integratingVolume = integratingVolume + simpleIntegral() * &
                (integrand(i) + integrand(i + 1)) / 2.0d0
        end do
        
        if (upperKnot /= spatialKnotsNumber) then
            upperR = (spatialKnots(upperKnot + 1) + spatialKnots(upperKnot)) / 2.0d0
            bottomR = spatialKnots(upperKnot)
            integratingVolume = integratingVolume + simpleIntegral() * &
                (3.0d0 / 4.0d0 * integrand(upperKnot) + &
                 1.0d0 / 4.0d0 * integrand(upperKnot + 1))
        end if
        
        contains 
        ! estimate int_{bottomR}^{upperR} 4 * PI r^{2} dr
        ! upperR >= bottomR
        function simpleIntegral()
            real(8) :: simpleIntegral 
            if (upperR < bottomR) then
                stop 'Calculations: upperR and bottomR must satisfy: upperR >= bottomR'
            end if
            simpleIntegral = 4.0d0 / 3.0d0 * M_PI * (upperR**(3.0d0) - bottomR**(3.0d0)) 
        end function
            
    end function
	
end module calculations

	