! Needs to be linked with: tridiagonal.f95 calculations.f95
module sphericalHeatEquation

    use calculations
    use tridiagonal
    
    implicit none
    
    type, public :: heatEquation
        integer, private :: fout = 20, ferror = 30
        integer, private :: foutBeta = 40
        
        logical, private :: isSolved
        
        real(8), private :: sigma
        
        real(8), private :: t0
        
        integer(4), private :: stepNumbers 
        integer(4), private :: timePeriodsNumber
        real(8), private, allocatable, dimension(:) :: timePeriods
        integer(4), private, allocatable, dimension(:) :: numbersOfPeriods
        
        integer(4), private :: spatialKnotsNumber
        real(8), private, allocatable, dimension(:) :: spatialKnots
        real(8), private, allocatable, dimension(:) :: initialCondition
        
        integer(4), private :: sectionsNumber
        integer(4), private, allocatable, dimension(:) :: sections
        real(8), private, allocatable, dimension(:) :: timeOfSections
        real(8), private, allocatable, dimension(:,:) :: outputSections 
        
    end type heatEquation
     
    real(8), public, parameter :: DAY_IN_SECONDS = 86400.0d0
    real(8), public, parameter :: YEAR_IN_SECONDS = 31556926.0d0
    real(8), public, parameter :: MEGAYEAR_IN_SECONDS = 0.3155692600000d14 
	
contains 
	
    subroutine hESolution(hE, getConductivity, getSources, getCapacity, getBeta, getSurfaceTemperatureFromBeta, & 
                          getRelativisticFactors, outputFile)
        interface 
		
			subroutine getRelativisticFactors (spatialKnotsNumber, redshiftFactor, volumeFactor)
                implicit none
                integer(4) :: spatialKnotsNumber
                real(8) :: redshiftFactor(1 : spatialKnotsNumber), volumeFactor(1 : spatialKnotsNumber)
			end subroutine getRelativisticFactors
		
            function getConductivity (thermal, spatialKnotsNumber, timeKnot)
                implicit none
                integer(4) :: spatialKnotsNumber
                real(8) :: getConductivity(2 : spatialKnotsNumber), thermal(spatialKnotsNumber)
                real(8) :: timeKnot
            end function getConductivity
           
            function getSources (thermal, spatialKnotsNumber, timeKnot)
            implicit none
                integer(4) :: spatialKnotsNumber
                real(8) :: getSources(spatialKnotsNumber), thermal(spatialKnotsNumber)
                real(8) :: timeKnot
            end function getSources
            
            function getCapacity (thermal, spatialKnotsNumber, timeKnot)
            implicit none
                integer(4) :: spatialKnotsNumber
                real(8) :: getCapacity(spatialKnotsNumber), thermal(spatialKnotsNumber)
                real(8) :: timeKnot
            end function getCapacity
            
            function getBeta (thermal, timeKnot) 
            implicit none
                real(8) :: getBeta, thermal, timeKnot
            end function getBeta
            
            subroutine getSurfaceTemperatureFromBeta(beta, T10g, surfaceTemperature, surfaceFlux)
                implicit none
                real(8) :: beta, T10g, surfaceTemperature, surfaceFlux
            end subroutine
        end interface
    
        type(heatEquation) :: hE  
        
        character(*) :: outputFile
        !for Windows 
        character(len = 50) :: surfaceFluxFile = 'outputData\\surfaceFlux.txt'
        character(len = 50) :: errorFile = 'outputData\\error.txt'
        !for Unix 
        ! character(len = 50) :: surfaceFluxFile = 'outputData//surfaceFlux.txt'
        ! character(len = 50) :: neutrinoLosses = 'outputData//neutrinoLosses.txt'
    
        integer :: i, j, k, l, m, n, counter, outputCounter, iterationCounter
        real(8) :: currentTime, currentTimeInYear
        real(8) :: tau
        real(8), allocatable, dimension(:,:) :: desk 
        real(8), allocatable, dimension(:) :: rightPart
        real(8), allocatable, dimension(:) :: mConductivity, mSources, mCapacity, sources, & 
                                              pConductivity, pSources, pCapacity, &
                                              meanSquare, mRedshiftFactor, mVolumeFactor,  h, hbar
        real(8), allocatable, dimension(:) :: previousLayer, currentLayer
        
        real(8), dimension(:) :: integrand(hE%spatialKnotsNumber) !for neutrino losses
        
        real(8) :: beta0, beta1, surfaceTemperature, surfaceFlux, prevSurfaceFlux
        !!!!
        real(8) :: prevCurrentSurfFlux, ratio, volumeHeatLoss, surfaceLoss, initialRight
        !!!!
        real(8) :: RdV     ! R_{0} / V_{0} + R_{1} / V_{1}
        real(8) :: meanRV  ! 0.75 * R_{0}V_{0} + 0.25 * R_{1}V_{1}
        real(8) :: waveR2  ! (r_{n}^{2} + (r_{n} - h_{n} / 2)^{2}) / 2
        
        logical :: isAccurate
        
        real(8), parameter :: FOUR_PI = 12.56637061d0
		
        type(tridiagonalMatrix) :: tM
        
        open(hE%ferror, file = errorFile)
        
        if (hE%isSolved .EQV. .false.) then
        
            m = hE%timePeriodsNumber
            n = hE%spatialKnotsNumber
            counter = 1
            outputCounter = 1
            currentTime = hE%t0

            allocate(h(2 : n), hbar(3 : n), meanSquare(2 : n), mConductivity(2 : n), mSources(n), mCapacity(n), &
                    pConductivity(2 : n), pSources(n), pCapacity(n), sources(n), &
                    mRedshiftFactor(n), mVolumeFactor(n), desk(3, n), rightPart(n), previousLayer(n), currentLayer(n))
            
			call getRelativisticFactors(n, mRedshiftFactor, mVolumeFactor)
			
            do i = 2, n
                h(i)          = hE%spatialKnots(i) - hE%spatialKnots(i - 1)
                meanSquare(i) = (hE%spatialKnots(i)**(2.0d0) * mRedshiftFactor(i) / mVolumeFactor(i) + & 
				                 hE%spatialKnots(i - 1)**(2.0d0) * mRedshiftFactor(i - 1) / mVolumeFactor(i - 1)) / 2.0d0
		    end do
            
            do i = 3, n
                hbar(i) = (h(i) + h(i - 1) ) / 2.0d0
            end do
            
            previousLayer = hE%initialCondition
            currentLayer = hE%initialCondition
            
            open(hE%fout, file = outputFile)
            open(hE%foutBeta, file = surfaceFluxFile)
            
            write(hE%fout, '(i8)') hE%sectionsNumber
            
            write(hE%foutBeta, '(a)') '   time                timeInYear       surfaceFlux    Tg = 10d10 g/cc       Ts     &
                                        neutrinoLosses'
            write(hE%foutBeta, '(i8)') hE%stepNumbers
            
            if ((outputCounter <= hE%sectionsNumber) .AND. (hE%sections(outputCounter) == counter)) then
                hE%timeOfSections(outputCounter) = currentTime
                hE%outputSections(1: n, outputCounter) = currentLayer
                outputCounter = outputCounter + 1
                
                ! Output section
                pConductivity = getConductivity(previousLayer, n, currentTime)
                pSources      = getSources     (previousLayer, n, currentTime)
                pCapacity     = getCapacity    (previousLayer, n, currentTime)
                currentTimeInYear = currentTime / YEAR_IN_SECONDS
                
                do i = 1, n
                    write(hE%fout, '(e22.13, e15.6, e15.6, e18.8, e15.6, e15.6, e15.6)') currentTime, currentTimeInYear, & 
                    hE%spatialKnots(i), hE%outputSections(i, outputCounter - 1), pCapacity(i), pConductivity(i), pSources(i)
                end do
                
                ! Let's write flux from star`s surface 
                beta1 = getBeta(currentLayer(n), currentTime)
                call getSurfaceTemperatureFromBeta(beta1, currentLayer(n), surfaceTemperature, surfaceFlux)
                
                !Let's write neutrino losses 
                do i = 1, hE%spatialKnotsNumber
                    integrand(i) = mRedshiftFactor(i) * pSources(i)
                end do
                
                write(hE%foutBeta, '(e22.13, e15.6, e15.6, e15.6, e15.6, e15.6)') currentTime, currentTimeInYear, &
                                    -surfaceFlux, currentLayer(n), surfaceTemperature, -estimateNeutrinoLosses()
            endif
            !for specialStep 
            initialRight = 20000.0d0
            
            isAccurate = .true.
			
            do j = 1, m
                tau = hE%timePeriods(j)
                do k = 1, hE%numbersOfPeriods(j)
                
                    counter = counter + 1
                    currentTime = currentTime + tau
                    previousLayer = currentLayer
                    
                    pConductivity = getConductivity(previousLayer, n, currentTime - tau)
                    pSources      = getSources     (previousLayer, n, currentTime - tau)
                    pCapacity     = getCapacity    (previousLayer, n, currentTime - tau)
                    
                    if (isAccurate) then
                        call step()
                        
                        mSources      = getSources     (currentLayer, n, currentTime)
                        mCapacity     = getCapacity    (currentLayer, n, currentTime)
                        
                        call estimateLosses()
                        isAccurate = validateAccuracy()
                        if (isAccurate .eqv. .false.) then
                            write(*, *) 'we are here' 
                            call monotonicity()
                        end if
                       
                        !!!!
                        if (MOD(counter, 100) == 1) then 
                            write(*, '(i8, e18.8)') counter, (volumeHeatLoss / surfaceLoss)
                        endif
                        write(hE%ferror, '(e18.8)') (volumeHeatLoss / surfaceLoss)
                        !!!!
                        
                    else
                        call specialStep()
                        !!!!   
                        if (MOD(counter, 100) == 1) then                        
                            call estimateLosses()
                            write(*, '(i8, e18.8)') counter, (volumeHeatLoss / surfaceLoss)
                        endif
                        !!!!                        
                    end if
                                 
                    !!!!
                    do i = 1, n
                        if (currentLayer(i) < 0.0d0) then
                            write(*, '(i8, e18.7, e18.7, e18.7)') counter, currentTime, hE%spatialKnots(i), currentLayer(i)
                            write(hE%ferror, '(a, a)') '      time             radius            thermal', &
                            '          capacity         conductivity       sources'                                
                            do l = 1, n 
                                write(hE%ferror, '(e18.7, e18.7, e18.7, e18.7, e18.7, e18.7)') currentTime - tau, & 
                                    hE%spatialKnots(l), previousLayer(l), pCapacity(l), pConductivity(l), pSources(l)
                            end do
                            do l = 1, n 
                                write(hE%ferror, '(e18.7, e18.7, e18.7, e18.7, e18.7, e18.7, e18.7)') currentTime, & 
                                    hE%spatialKnots(l), currentLayer(l),  mCapacity(l), mConductivity(l), mSources(l), sources(l)
                            end do
                        end if
                    end do
                    !!!!
                    
                    ! Output section
                    currentTimeInYear = currentTime / YEAR_IN_SECONDS
                    ! Let's write flux from star`s surface 
                    beta1 = getBeta(currentLayer(n), currentTime)
                    call getSurfaceTemperatureFromBeta(beta1, currentLayer(n), surfaceTemperature, surfaceFlux)
                    
                    !Let's write neutrino losses 
                    do i = 1, hE%spatialKnotsNumber
                        integrand(i) = mRedshiftFactor(i) * mSources(i)
                    end do
                 
                    write(hE%foutBeta, '(e22.13, e15.6, e15.6, e15.6, e15.6, e15.6)') currentTime, currentTimeInYear, &
                                        -surfaceFlux, currentLayer(n), surfaceTemperature, -estimateNeutrinoLosses()
                    prevSurfaceFlux = surfaceFlux
                    !Let's write selected sections
                    if ((outputCounter <= hE%sectionsNumber) .AND. (hE%sections(outputCounter) == counter)) then
                        hE%timeOfSections(outputCounter) = currentTime
                        hE%outputSections(1: n, outputCounter) = currentLayer
                        outputCounter = outputCounter + 1
                        
                        do i = 1, n
                            write(hE%fout, '(e22.13, e15.6, e15.6, e18.8, e15.6, e15.6, e15.6)') currentTime, & 
                            currentTimeInYear, hE%spatialKnots(i), hE%outputSections(i, outputCounter - 1), & 
                            mCapacity(i), mConductivity(i), mSources(i)
                        end do
                    endif
                
                end do
            
            end do
            hE%isSolved = .true.
            !!!!
            close(hE%ferror)
            !!!!
            close(hE%fout)
            close(hE%foutBeta)
            
        end if
        
        contains
        
        subroutine step()       
            mConductivity = getConductivity(currentLayer, n, currentTime)
            mSources      = getSources     (currentLayer, n, currentTime)
            mCapacity     = getCapacity    (currentLayer, n, currentTime)
			
            sources = (1 - hE%sigma) * pSources + hE%sigma * mSources
                
            meanRV = 0.75d0 * (mRedshiftFactor(1) * mVolumeFactor(1)) + & 
                     0.25d0 * (mRedshiftFactor(2) * mVolumeFactor(2))
            RdV = mRedshiftFactor(1) / mVolumeFactor(1) + mRedshiftFactor(2) / mVolumeFactor(2)

            desk(1, 1) = 0.0d0
            desk(2, 1) = 3.0d0 * hE%sigma * mConductivity(2) * mRedshiftFactor(1) * RdV / (h(2) * h(2)) + & 
                         mCapacity(1) * meanRV / tau
            desk(3, 1) = (-3.0) * hE%sigma * mConductivity(2) * mRedshiftFactor(2) * RdV / h(2) / h(2)
            rightPart(1) = 3.0 * (1.0 - hE%sigma) * pConductivity(2) * RdV * &
                           (mRedshiftFactor(2) * previousLayer(2) - mRedshiftFactor(1) * previousLayer(1) ) / h (2) / h(2)
            rightPart(1) = rightPart(1) + (0.75d0 * mRedshiftFactor(1)**(2.0d0) * mVolumeFactor(1) + &
        				                   0.25d0 * mRedshiftFactor(2)**(2.0d0) * mVolumeFactor(2)) * sources(1) + &
                                           meanRV * previousLayer(1) * pCapacity(1) / tau
                
            do i = 2, n - 1
                desk(1, i) = (-hE%sigma) * mRedshiftFactor(i - 1) * mConductivity(i) * meanSquare(i) / &
                             h(i) / hbar(i + 1)
                desk(2, i) = hE%sigma / hbar(i + 1) * (mConductivity(i) * meanSquare(i) / h(i) + &
                             mConductivity(i + 1) * meanSquare(i + 1) / h(i + 1) ) + & 
                             mCapacity(i) * mVolumeFactor(i) * hE%spatialKnots(i)**(2.0) / tau
                desk(2, i) = desk(2, i) * mRedshiftFactor(i)
                desk(3, i) = mRedshiftFactor(i + 1) * (-hE%sigma) * meanSquare(i + 1) * mConductivity(i + 1) / &
					         h(i + 1) / hbar(i + 1)
                rightPart(i) = (1.0 - hE%sigma) / hbar(i + 1) *                           &
                               ((mRedshiftFactor(i + 1) * previousLayer(i + 1) -         &
                                 mRedshiftFactor(i) * previousLayer(i) ) / h(i + 1) *    &
                                 pConductivity(i + 1) * meanSquare(i + 1) -                 &
                                (mRedshiftFactor(i) * previousLayer(i) -                  &
						         mRedshiftFactor(i - 1) * previousLayer(i - 1) ) / h(i) * &
                                 pConductivity(i) * meanSquare(i) )
                rightPart(i) = rightPart(i) + ((sources(i) * mRedshiftFactor(i) + & 
					                            previousLayer(i) * pCapacity(i) / tau) * &
                               hE%spatialKnots(i)**(2.0) * mRedshiftFactor(i) * mVolumeFactor(i) )
            end do
                
            beta0 = getBeta(previousLayer(n), currentTime - tau)
            beta1 = getBeta(currentLayer(n), currentTime)
            meanRV = (0.75d0 * mRedshiftFactor(n) * mVolumeFactor(n) + & 
                      0.25d0 * mRedshiftFactor(n - 1) * mVolumeFactor(n - 1))
            waveR2 = (hE%spatialKnots(n)**(2.0d0) + (hE%spatialKnots(n) - (h(n) / 2.0d0) )**(2.0d0) ) / 2.0d0
				
            desk(1, n) = (-hE%sigma) * mRedshiftFactor(n - 1) * mConductivity(n) * meanSquare(n) / h(n) / h(n)
            desk(2, n) = hE%sigma / h(n) * &
                         (mRedshiftFactor(n) * mConductivity(n) * meanSquare(n) / h(n) - beta1)
            desk(2, n) = desk(2, n) + (mCapacity(n) / 2.0d0 / tau * waveR2 * meanRV)
            desk(3, n) = 0.0
            rightPart(n) = (1 - hE%sigma) / h(n) * (beta0 * previousLayer(n) - &
                           pConductivity(n) / h(n) * meanSquare(n) * &
                           (mRedshiftFactor(n) * previousLayer(n) - mRedshiftFactor(n - 1) * previousLayer(n - 1) ) ) 
            rightPart(n) = rightPart(n) + (waveR2 / 2.0d0) * (meanRV * pCapacity(n) * previousLayer(n) / tau + &
				           (0.75d0 * mRedshiftFactor(n)**(2.0d0) * mVolumeFactor(n) + &
                            0.25d0 * mRedshiftFactor(n - 1)**(2.0d0) * mVolumeFactor(n - 1) ) * sources(n) )
                
            call tMConstructor(tM, desk, rightPart, hE%spatialKnotsNumber)
            currentLayer = tMsolution(tM)
            
        end subroutine
        
        subroutine specialStep()
            real(8) :: surfaceDTLeft, surfaceDTRight, surfaceDT
            real(8) :: rightRatio, leftRatio
            real(8) :: gap
            integer(4) :: gapCounter
            
            surfaceDTLeft = 0.0d0
            surfaceDTRight = initialRight
            
            ! Let's estimate surfaceDTRight and rightRatio
            do i = 1, n
                currentLayer(i) = previousLayer(i) * (1.0 - surfaceDTRight / previousLayer(n))
            end do
            
            mCapacity = getCapacity(currentLayer, n, currentTime)
            rightRatio = estimateLossesWithoutNeutrino()
            
            do while (rightRatio > 4.0d0)
                surfaceDTRight = surfaceDTRight / 2.0d0
                do i = 1, n
                    currentLayer(i) = previousLayer(i) * (1.0 - surfaceDTRight / previousLayer(n))
                end do
            
                mCapacity = getCapacity(currentLayer, n, currentTime)
                rightRatio = estimateLossesWithoutNeutrino()
            end do
            initialRight = surfaceDTRight
            
            do while (rightRatio < 1.0d0)
                surfaceDTRight = surfaceDTRight * 2.0d0
                do i = 1, n
                    currentLayer(i) = previousLayer(i) * (1.0d0 - surfaceDTRight / previousLayer(n))
                end do
                mCapacity = getCapacity(currentLayer, n, currentTime)
                rightRatio = estimateLossesWithoutNeutrino()
            end do
            
            ! Let's estimate leftRatio 
            leftRatio = 0.0d0
                     
            gap = 0.0003d0
            gapCounter = 0
            ratio = 2.0d0
            do while ((ratio > 1.0d0 + gap) .OR. (ratio < 1.0d0 - gap))
                gapCounter = gapCounter + 1
                !!!!
                if (MOD(counter, 30) == 1) then 
                    write(*, *) counter, 'ratio: ', ratio
                end if
                !!!!
                surfaceDT = surfaceDTLeft + (surfaceDTRight - surfaceDTLeft) * (1.0d0 - leftRatio) / &
                            (rightRatio - leftRatio)
                do i = 1, n
                    currentLayer(i) = previousLayer(i) * (1.0 - surfaceDT / previousLayer(n))
                end do
                mCapacity = getCapacity(currentLayer, n, currentTime)
                ratio = estimateLossesWithoutNeutrino()
                
                if (ratio > 1.0d0) then
                    surfaceDTRight = surfaceDT
                    rightRatio = ratio
                else
                    surfaceDTLeft = surfaceDT
                    leftRatio = ratio
                end if
                
                if (gapCounter == 50) then
                    gapCounter = 0
                    gap = gap * 10.0d0
                end if
                    
            end do
            
            mConductivity = getConductivity(currentLayer, n, currentTime)
            mSources      = getSources     (currentLayer, n, currentTime)
            
        end subroutine
        
        subroutine estimateLosses()
            ! numerical losses
            do i = 1, n
                integrand(i) = (currentLayer(i) - previousLayer(i)) * (mCapacity(i) + pCapacity(i)) / 2.0d0 * & 
                                mRedshiftFactor(i) * mVolumeFactor(i)
            end do
            volumeHeatLoss = integratingVolume(n, n, 1, hE%spatialKnots, integrand)
                         
            beta1 = getBeta(currentLayer(n), currentTime)
            call getSurfaceTemperatureFromBeta(beta1, currentLayer(n), surfaceTemperature, surfaceFlux)
            beta1 = getBeta(previousLayer(n), currentTime)
            call getSurfaceTemperatureFromBeta(beta1, previousLayer(n), surfaceTemperature, prevSurfaceFlux)
            surfaceLoss = (prevSurfaceFlux + surfaceFlux) / 2.0d0 * tau
                        
            do i = 1, n
                integrand(i) = mRedshiftFactor(i)**(2.0d0) * mVolumeFactor(i) * (mSources(i) + pSources(i)) / 2.0d0
            end do
            surfaceLoss = surfaceLoss + integratingVolume(n, n, 1, hE%spatialKnots, integrand) * tau
            
        end subroutine
        
        function estimateLossesWithoutNeutrino()
            real(8) :: estimateLossesWithoutNeutrino
            ! physical losses
            beta1 = getBeta(currentLayer(n), currentTime)
            call getSurfaceTemperatureFromBeta(beta1, currentLayer(n), surfaceTemperature, surfaceFlux)
            beta1 = getBeta(previousLayer(n), currentTime)
            call getSurfaceTemperatureFromBeta(beta1, previousLayer(n), surfaceTemperature, prevSurfaceFlux)
            surfaceLoss = (prevSurfaceFlux + surfaceFlux) / 2.0d0 * tau
            do i = 1, n
                integrand(i) = mRedshiftFactor(i)**(2.0d0) * mVolumeFactor(i) * (mSources(i) + pSources(i)) / 2.0d0
            end do
            surfaceLoss = surfaceLoss + integratingVolume(n, n, 1, hE%spatialKnots, integrand) * tau
                    
            ! numerical losses
            do i = 1, n
                integrand(i) = (currentLayer(i) - previousLayer(i)) * (mCapacity(i) + pCapacity(i)) / 2.0d0 * & 
                                mRedshiftFactor(i) * mVolumeFactor(i)
            end do
            volumeHeatLoss = integratingVolume(n, n, 1, hE%spatialKnots, integrand)
                    
            estimateLossesWithoutNeutrino = volumeHeatLoss / surfaceLoss
            
        end function
        
        function validateAccuracy()
            logical :: validateAccuracy
            if (currentTime < 0.2 * MEGAYEAR_IN_SECONDS) then
                validateAccuracy = .true.
            else if (surfaceFlux / prevSurfaceFlux > 1.001) then
                validateAccuracy = .false.
            end if
        end function
        
        function estimateNeutrinoLosses()
            real(8) :: estimateNeutrinoLosses
        
            estimateNeutrinoLosses = integratingVolume(n, n, 1, hE%spatialKnots, integrand)
        
        end function
        
        subroutine monotonicity()
            real(8) :: key
            
            integer(4) :: p
            
            do i = 2, n
                key = currentLayer(i)
                p = i
                do while((p > 1) .AND. (currentLayer(p - 1) < key))
                    currentLayer(p) = currentLayer(p - 1)
                    p = p - 1
                end do
                currentLayer(p) = key 
            end do
            
        end subroutine
        
        function monotonicityChecking()
            logical monotonicityChecking
            
            monotonicityChecking = .FALSE.
            do i = hE%spatialKnotsNumber, 2, -1
                if (currentLayer(i) > currentLayer(i - 1)) then
                    monotonicityChecking = .TRUE.
                end if
            end do
        end function
       
    end subroutine  
    
    function getSHETimeOfSections(hE)
        type(heatEquation) :: hE
        real(8) :: getSHETimeOfSections(hE%sectionsNumber)
        getSHETimeOfSections = hE%timeOfSections
    end function
    
    function getSHEOutputSections(hE)
    type(heatEquation) :: hE
        real(8) :: getSHEOutputSections(hE%spatialKnotsNumber, hE%sectionsNumber)
        getSHEOutputSections = hE%outputSections
    end function
    
    subroutine setSHEInitialConditions(hE, spatialKnotsNumber, sigma, t0, spatialKnots, initialCondition) 
        type(heatEquation) :: hE
        
        integer(4) :: spatialKnotsNumber
        real(8) :: sigma, t0
        real(8) :: spatialKnots(spatialKnotsNumber), initialCondition(spatialKnotsNumber)
        
        hE%spatialKnotsNumber = spatialKnotsNumber
        hE%sigma = sigma
        hE%t0 = t0
        
        allocate(hE%spatialKnots(spatialKnotsNumber), hE%initialCondition(spatialKnotsNumber))
        hE%spatialKnots = spatialKnots
        hE%initialCondition = initialCondition
    end subroutine
	
    subroutine setSHEStepConditions(hE, stepNumbers, timePeriodsNumber, timePeriods, numbersOfPeriods)
        type(heatEquation) :: hE
        
        integer(4) :: stepNumbers, timePeriodsNumber
    
        real(8) :: timePeriods(timePeriodsNumber)
        integer(4) :: numbersOfPeriods(timePeriodsNumber)
        
        hE%stepNumbers = stepNumbers
        hE%timePeriodsNumber = timePeriodsNumber
        
        allocate(hE%timePeriods(timePeriodsNumber), hE%numbersOfPeriods(timePeriodsNumber))
        hE%timePeriods = timePeriods
        hE%numbersOfPeriods = numbersOfPeriods
    end subroutine
    
    subroutine setSHEOutputConditions(hE, sectionsNumber, sections)
        type(heatEquation) :: hE
        
        integer(4) :: sectionsNumber
        integer(4) :: sections(sectionsNumber)
        
        hE%sectionsNumber = sectionsNumber
        allocate(hE%sections(sectionsNumber), hE%timeOfSections(sectionsNumber))
        allocate(hE%outputSections(hE%spatialKnotsNumber, sectionsNumber))
        hE%sections = sections
    end subroutine
    
end module