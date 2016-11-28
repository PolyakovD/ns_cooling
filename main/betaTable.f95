! Needs to be linked with: tbtsfb.f

module betaTable
	
implicit none

    integer, private :: m = 0 ! massKnotsNumber
    integer, private :: n = 0 ! lgTKnotsNumber

    real(8), private, allocatable :: masses(:)
    real(8), private, allocatable :: lgTb(:)

    real(8), private, allocatable :: betaValues(:,:) ! surface flux from rho_b surface divided by (4 * Pi)
    real(8), private, allocatable :: lgTsValues(:,:)
    
    real(8), private, allocatable :: buffer(:) ! for input trick
    
    real(8), private :: rb ! the radius on the border of the grid
    real(8), private :: stellarMass ! stellarMass
    real(8), private :: stellarRadius ! radius on the surface in cm > spatialKnots(spatialKnotsNumber)
    real(8), private :: g14 
    
    integer, private :: buttomKnot, upperKnot ! masses(buttomKnot) < stellarMass < masses(upperKnot)
    real(8), private :: bMWeight, uMWeight 
    
    real(8), private, parameter :: SUN_MASS = 1.989d+33
    real(8), private, parameter :: STEFAN_BOLTZMANN = 5.67d-5
    real(8), private, parameter :: M_PI = 3.141592653589793d0 

contains

    subroutine constructGridFromFile(input, curRb, curStellarMass, curStellarRadius, curG14)
        character(*) :: input
        real(8) :: curRb, curStellarMass, curStellarRadius, curG14
        integer, parameter :: fin = 101
        
        character :: absorbDescription = 'a'
        integer :: i
        real(8) :: betaFactor
        
        ! check block
        if ((curRb < 9.0d+5) .OR. (curRb > 1.4d+6)) then
            write(*,*) 'betaTable:Warning: detected an unexpected value of the parameter rb = ', curRb
        end if
        if ((curStellarMass < 1.989d+32) .OR. (curStellarMass > 4.376d+33)) then
            write(*,*) 'betaTable:Warning: detected an unexpected value of the parameter stellarMass = ', curStellarMass
        end if
        if ((curStellarRadius < 9.0d+5) .OR. (curStellarRadius > 1.4d+6)) then
            write(*,*) 'betaTable:Warning: detected an unexpected value of the parameter stellarRadius = ', curStellarRadius
        end if 
        if ((curG14 < 1.0d-4) .OR. (curG14 > 10.0)) then
            write(*,*) 'betaTable:Warning: detected an unexpected value of the parameter g14 = ', curG14
        end if
        
        rb = curRb
        stellarMass = curStellarMass
        stellarRadius = curStellarRadius
        g14 = curG14
        
        open(fin, file = input)
        do while (absorbDescription .NE. '#')
            read (fin, *) absorbDescription
        end do
        
        read(fin, *) n, m 
        
        allocate(masses(2), lgTb(n), buffer(3 * m), & 
                 betaValues(n, 2), lgTsValues(n, 2))
         
        read(fin, *) masses
        do i = 1, m
            masses(i) = masses(i) * SUN_MASS 
        end do
        
        if ((masses(1) > curStellarMass) .OR. (masses(m) < curStellarMass)) then
            stop 'model: stellarMass is out of approximation range'
        end if
        do i = 1, m - 1 
            if (masses(i + 1) > curStellarMass) then
                buttomKnot = i
                upperKnot = i + 1
                exit
            end if
        end do
        bMWeight = masses(upperKnot) - curStellarMass
        uMWeight = curStellarMass - masses(buttomKnot)

        read(fin, *) absorbDescription
        
        betaFactor = 1.0d0 / 12.5663706144d0 ! 1 / (4 * Pi)
        do i = 1, n
            read(fin, *) lgTb(i), buffer
            lgTsValues(i, 1) = buffer(3 * buttomKnot - 1)
            lgTsValues(i, 2) = buffer(3 * upperKnot - 1)
            betaValues(i, 1) = buffer(3 * buttomKnot) * betaFactor  
            betaValues(i, 2) = buffer(3 * upperKnot) * betaFactor
        end do
        close(fin)
    end subroutine
    
    ! Input: beta
    !        T10g - temperature in K at the layer with rho = 10^{8}g/cc
    ! Output : surfaceTemperature in K
    !          surfaceFlux in erg/s from rho_b surface
	subroutine getSurfaceTemperatureFromBeta04(beta, T8g, surfaceTemperature, surfaceFlux)
        real(8) :: beta, T8g
        real(8) :: surfaceTemperature, surfaceFlux
		
        integer :: i
        real (8) :: upperWeight, buttomWeight
        
        real(8) :: lgPoint
        
        lgPoint = DLOG10(T8g)
        if (lgPoint > lgTb(n)) then
            stop 'model: stellar border temperature is out of approximation range'
        end if
        if (lgPoint < lgTb(1)) then 
            surfaceTemperature = T8g 
            surfaceFlux = -surfaceTemperature**4 * STEFAN_BOLTZMANN * 12.5663706144d0 * stellarRadius**2 
        else 
            i = 1
            do while(lgPoint > lgTb(i))
                i = i + 1
            end do
            buttomWeight = lgTb(i) - lgPoint
            upperWeight = lgPoint - lgTb(i - 1)
            surfaceTemperature = (lgTsValues(i, 1) * upperWeight + lgTsValues(i - 1, 1) * buttomWeight) / & 
                    (upperWeight + buttomWeight) * bMWeight + &
                    (lgTsValues(i, 2) * upperWeight + lgTsValues(i - 1, 2) * buttomWeight) / & 
                    (upperWeight + buttomWeight) * uMWeight
            surfaceTemperature = 10.0d0**(surfaceTemperature / (bMWeight + uMWeight))
            
            surfaceFlux = (betaValues(i, 1) * upperWeight + betaValues(i - 1, 1) * buttomWeight) / & 
                    (upperWeight + buttomWeight) * bMWeight + &
                    (betaValues(i, 2) * upperWeight + betaValues(i - 1, 2) * buttomWeight) / & 
                    (upperWeight + buttomWeight) * uMWeight
            surfaceFlux = surfaceFlux / (bMWeight + uMWeight)
            surfaceFlux = -surfaceFlux * 12.5663706144d0 ! beta04 * 4 * Pi
        end if   
		
	end subroutine
    
    ! estimate beta in shallow model (rho_{b} = 10^{8}g/cc)
    ! Input: thermal - temperature in rho_{b} = 10^{8}g/cc in K
    !        timeKnot - the time from the beginning of the simulation      
    function beta04(thermal, timeKnot)
        real(8) :: beta04
        real(8) :: thermal, timeKnot
        integer :: i
        real (8) :: upperWeight, buttomWeight
        
        real(8) :: lgPoint
        
        lgPoint = DLOG10(thermal)
        
        if (lgPoint > lgTb(n)) then
            stop 'model: stellar border temperature is out of approximation range'
        end if
        if (lgPoint < lgTb(1)) then 
            beta04 = -thermal**3 * STEFAN_BOLTZMANN * rb**2
        else 
            i = 1
            do while (lgPoint > lgTb(i))
                i = i + 1
            end do
            buttomWeight = lgTb(i) - lgPoint
            upperWeight = lgPoint - lgTb(i - 1)
            beta04 = (betaValues(i, 1) * upperWeight + betaValues(i - 1, 1) * buttomWeight) / & 
                    (upperWeight + buttomWeight) * bMWeight + &
                    (betaValues(i, 2) * upperWeight + betaValues(i - 1, 2) * buttomWeight) / & 
                    (upperWeight + buttomWeight) * uMWeight
            beta04 = beta04 / (bMWeight + uMWeight)
            beta04 = -beta04 / thermal
        end if
    end function
    
    function beta03(thermal, timeKnot)
        real(8) :: beta03, thermal, timeKnot
		
        real(8) :: FB, B12, COSLAT, TS6
		
        COSLAT = 0.0d0
        B12 = 0.0d0
		
        call TbTsFb(g14, thermal / 1.0d9, B12, COSLAT, TS6, FB) ! from tbtsfb
        beta03 = -FB * stellarRadius**2 / thermal
		
    end function
	
	! Input: beta
    !        T10g - temperature in K at the layer with rho = 10^{10}g/cc
    ! Output : surfaceTemperature in K
    !          surfaceFlux in erg/s from rho_b surface
	subroutine getSurfaceTemperatureFromBeta03(beta, T10g, surfaceTemperature, surfaceFlux)
        real(8) :: beta, T10g
        real(8) :: surfaceTemperature, surfaceFlux
		
		real(8) B12, COSLAT, FB
		
		COSLAT = 0.0d0
		B12 = 0.0d0
		
		call TbTsFb(g14, T10g / 1.0d9, B12, COSLAT, surfaceTemperature, FB)
		surfaceTemperature = surfaceTemperature * 1.0d6
		surfaceFlux = -FB * 4.0d0 * M_PI * stellarRadius**2 
		
	end subroutine
    
    ! Deprecated
	! Estimate beta for surface flux GPE(1983): Eq.(32)
    function beta02(thermal, timeKnot)
        real(8) :: beta02, thermal, timeKnot
		
        beta02 = -stellarRadius**2 * STEFAN_BOLTZMANN * thermal**(1.1978d0)
        beta02 = beta02 * g14 * 10.0d0**(6.4176d0) / 1.7441
        
    end function
    
    ! Deprecated
    ! Input: beta
    !        T10g - temperature in K at the layer with rho = 10^{10}g/cc
    ! Output : surfaceTemperature in K
    !          surfaceFlux in erg/s
    subroutine getSurfaceTemperatureFromBeta02(beta, T10g, surfaceTemperature, surfaceFlux)
        real(8) :: beta, T10g
        real(8) :: surfaceTemperature, surfaceFlux

        surfaceFlux = beta * T10g * 4.0d0 * M_PI
        surfaceTemperature = (-surfaceFlux / (4.0d0 * M_PI * stellarRadius**2 * STEFAN_BOLTZMANN))**(0.25d0)
        
    end subroutine
end module