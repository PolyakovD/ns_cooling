! Needs to be linked with: mechanicalModel.f95 mechanicalModelIO.f95 sfgaps.f betaTable.f95 
!                          model.f95 sphericalHeatEquation.f95 sphericalHeatEquation.f95
! This module is interface for model
module modelIO

use mechanicalModelIO 
use model
use betaTable
use sphericalHeatEquationIO
use sphericalHeatEquation

implicit none
    
    type(mModel), private :: mM
    type(heatEquation), private :: hE
    
    real(8), private :: sigma, initialTau
    
    integer(4), private :: spatialKnotsNumber
    
    integer(4), private :: KEOS, isDeep
    
    real(8), private, allocatable :: rhoDistribution(:) ! for SF curves
    real(8), private, allocatable :: spatialKnots(:) ! for SHE conditions
	character(len = 2000) :: mMDescription ! for knit subroutines
    character(len = 200) :: SFDescription  ! for knit subroutines
    integer(4) :: neutronTripletSFModel, neutronSingletSFModel, protonSFModel !for knit subroutines
    
    ! input files
    character(len = 50), private :: configFile = 'inputData\\configFile.txt'                ! file with (for \sigma, \tau_{0}, frozenStepKey and lockalHeatingKey)
    character(len = 50), private :: mechanicalInput = 'mechanicalData\\mechanicalInput.txt' ! file with KEOS, rhoC, EPSMASS
    character(len = 50), private :: borderConditionInput = 'tb_ts_rho8_fixed.dat'           ! config file for border conditions 
    character(len = 50), private :: superfluidFile = 'inputData\\superfluid.txt'            ! file with numbers of superfluid model for nt, ns and ps respectively
    character(len = 50), private :: frozenStepFile = 'inputData\\frozenStep.txt'            ! file with frozen conditions
    character(len = 50), private :: lockalHeatingFile = 'inputData\\lockalHeating.txt'      ! file with lockal heating configuration
    ! output files
    character(len = 50), private :: mechanicalOutput = 'outputData\\mechanicalOutput.txt'   ! file for mechanical structure
    character(len = 50), private :: fileForHeatEquation = 'inputData\\conditions.txt'       ! file for SHE
    character(len = 50), private :: outputfile = 'outputData\\result.txt'              ! We will write here the solution with rhoDistribution
    !!!! only for testModel
    character(len = 50), private :: heatOutputfile = 'outputData\\heatResult.txt'            ! file with thermal solution U(r)
    !!!!
    !character(len = 50), private :: fluxFile = 'outputData\\surfaceFlux.txt'                ! file with thermal flux and neutrino losses
    character(len = 50), private :: outputfileForFlux = 'outputData\\resultFlux.txt'         ! we will write here flux and losses with description
       
    integer(4), private, parameter :: fout = 101, fin = 201
    integer(4) :: i, frozenKey, lockalKey, SFKey

contains

    subroutine conditionGeneratorFromFiles()
        !Let's estimate mechanical model and write it to file
        call estimateMM(.TRUE.)
        
        ! Let's read configuration
        open(fin, file = configFile)
        read(fin, *) sigma, initialTau, frozenKey, lockalKey, SFKey
        close(fin)
        
        ! Let's create the conditions for SHE 
        call createConditionsForSHE()
        
        !Let's read flashes parameters
        call readFlashesParameters()
        
        !Let's estimate critical temperature for superfluidity
        call estimateSFCriticalTemperature(.TRUE.)
    end subroutine
    
    subroutine conditionGeneratorFromCL()
        character :: answer
        !Let's estimate mechanical model and write it to file
        call estimateMM(.FALSE.)
        
        ! Let's read configuration
        write(*, *) 'Please, set the parameters of problem'
        
        write(*, *) 'Enter the weight sigma'
        read(*, *) sigma
        do while ((sigma < 0.0d0) .OR. (sigma > 1.0d0))
            write(*, *) 'sigma must lie in interval (0, 1)'
            read(*, *) sigma         
        end do
        
        write(*, *) 'Enter the intial dt'
        read(*, *) initialTau
        do while (initialTau < 0.0d0)
            write(*, *) 'initialTau must be positive'
            read(*, *) initialTau       
        end do
        
        write(*, *) 'Write "y" if you have frozen initial step in "frozenStep.txt" or "n" if you haven`t'
        read(*, *) answer
        do while ((answer .NE. 'y') .AND. (answer .NE. 'n'))
            write(*, *) 'Try again'
            read(*, *) answer
        end do
        if (answer .EQ. 'y') then
            frozenKey = 1
        else 
            frozenKey = 0
        end if
        
        write(*, *) 'Write "y" if you have source of heating in "lockalHeating.txt" or "n" if you haven`t'
        read(*, *) answer
        do while ((answer .NE. 'y') .AND. (answer .NE. 'n'))
            write(*, *) 'Try again'
            read(*, *) answer
        end do
        if (answer .EQ. 'y') then
            lockalKey = 1
        else 
            lockalKey = 0
        end if
        
        write(*, *) 'Write "y" if you want to take into account the superfluidity, "n" if yoy don`t want'
        read(*, *) answer
        do while ((answer .NE. 'y') .AND. (answer .NE. 'n'))
            write(*, *) 'Try again'
            read(*, *) answer
        end do
        if (answer .EQ. 'y') then
            SFKey = 1
        else 
            SFKey = 0
        end if
        
        ! Let's create the conditions for SHE 
        call createConditionsForSHE()
        
        !Let's read flashes parameters
        call readFlashesParameters()
        
        !Let's estimate critical temperature for superfluidity
        call estimateSFCriticalTemperature(.FALSE.)
    end subroutine
    
    ! prepares difference scheme ans runs the solution
    subroutine makeSolution()        
        real(8) :: stellarMass, g14, stellarRadius ! for border conditions
        
        stellarMass = getStellarMass(mM) ! for borderCondition grid
        g14 = getG14(mM) ! for old beta approximation
        stellarRadius = getRadius(mM) ! for old beta approximation
       
        ! construct border condition
        call constructGridFromFile(borderConditionInput, spatialKnots(spatialKnotsNumber), &
                                   stellarMass, stellarRadius, g14)
        
        call constructSHEFromFile(hE, fileForHeatEquation)
        call preapareFiles() ! write mechanical and SF parametera to result file
        if (isDeep .EQ. 1) then
            call hESolution(hE, conductivity02, sources02, capacity02, beta03, getSurfaceTemperatureFromBeta03, &
                            relativisticFactors02, outputfile, outputfileForFlux)
        else 
            call hESolution(hE, conductivity02, sources02, capacity02, beta04, getSurfaceTemperatureFromBeta04, &
                            relativisticFactors02, outputfile, outputfileForFlux)
        end if
        !call knitSolutions()
        !call knitFlux()     
    end subroutine
    
    ! estimate mechanical model, set it to model.f95 and write to file
    subroutine estimateMM(fileIsSource)
        logical :: fileIsSource
        real(8), allocatable, dimension(:) :: pressureDistribution, &
                                              redshiftFactor, volumeFactor ! output distributions
        real(8) :: radius, g14
        
        if (fileIsSource) then
            call constructMMFromFile(mM, mechanicalInput)
        else
            call constructMMFromCL(mM)
        end if
        
        call validateModel(mM)
		mMDescription = getMMDescription(mM)
        call getMMStructure(mM, spatialKnotsNumber, spatialKnots, rhoDistribution, pressureDistribution, &
                          redshiftFactor, volumeFactor, KEOS, isDeep, radius, g14)
        
        call setMMechanical(spatialKnotsNumber, spatialKnots, rhoDistribution, pressureDistribution, &
                          redshiftFactor, volumeFactor, KEOS, radius, g14)
        call printMechanicalModel()        
    end subroutine
    
    ! create conditions for SHE
    subroutine createConditionsForSHE()
        integer(4), parameter :: frozenFin = 301
        integer(4) :: frozenSpatialKnotsNumber
        real(8) :: time, timeInYear, radius, rho, thermal, capacity, conductivity, sources, gridFluxes
        
        open(fout, file = fileForHeatEquation)
        write(fout, *) sigma, initialTau
        
        if (frozenKey == 1) then
            open(frozenFin, file = frozenStepFile)
            read(frozenFin, *) frozenSpatialKnotsNumber
            
            if (frozenSpatialKnotsNumber /= spatialKnotsNumber) then
                write(*, *) 'modelIO: Mechanical model have ', spatialKnotsNumber, ' knots, but initial conditions have', &
                             frozenSpatialKnotsNumber, ' knots.'
                stop 'model: Initial conditions inconsistent with mechanical model'
            end if
            
            read(frozenFin, *) time, timeInYear, radius, rho, thermal, capacity, conductivity, sources, gridFluxes
            write(fout, '(i8, 5x, e19.8)') spatialKnotsNumber, time
            write(fout, '(e18.7, 5x, e25.12, 5x, e25.12)') radius, rho,  thermal
            do i = 2, spatialKnotsNumber
                read(frozenFin, *) time, timeInYear, radius, rho, thermal, capacity, conductivity, sources, gridFluxes
                write(fout, '(e18.7, 5x, e25.12, 5x, e25.12)') radius, rho, thermal
            end do
            
            close(frozenFin)
        else 
            write(fout, '(i8, 5x, e12.7)') spatialKnotsNumber, 0.0d0
            !initial distribution
            do i = 1, spatialKnotsNumber
                write(fout, '(e18.7, 5x, e25.12, 5x, e25.12)') spatialKnots(i), rhoDistribution(i), 1.0d10
            end do
        end if
        close(fout)
    end subroutine
    
    ! read flashes parameters and set it in model.f95
    subroutine readFlashesParameters()
        integer(4), parameter :: lockalFin = 401
        
        integer(4) :: flashesNumber ! number of flashes
        
        real(8), allocatable, dimension(:) :: flashesEnergy 
        real(8), allocatable, dimension(:) :: flashesEmissivity
        real(8), allocatable, dimension(:) :: flashesStartTime ! flashes start time from the star birth
        real(8), allocatable, dimension(:) :: flashesFinishTime ! flashes finish time from the star birth
        real(8), allocatable, dimension(:) :: flashesUpperRho
        integer(4), allocatable, dimension(:) :: flashesUpperKnots
        real(8), allocatable, dimension(:) :: flashesBottomRho
        integer(4), allocatable, dimension(:) :: flashesBottomKnots
        
        real(8) :: years, days, seconds 
        real(8) :: dt, QLockal
        integer(4) :: upperI, bottomI
        
        if (lockalKey == 1) then
            open(lockalFin, file = lockalHeatingFile)
            
            read(lockalFin, *) flashesNumber
            
            allocate(flashesEnergy(flashesNumber), flashesStartTime(flashesNumber), flashesFinishTime(flashesNumber), &
                     flashesUpperRho(flashesNumber), flashesBottomRho(flashesNumber), flashesEmissivity(flashesNumber), &
                     flashesUpperKnots(flashesNumber), flashesBottomKnots(flashesNumber))
            
            do i = 1, flashesNumber
                read(lockalFin, *) flashesEnergy(i)
                read(lockalFin, *) years, days, seconds, dt
                flashesStartTime(i) = years * YEAR_IN_SECONDS + days * DAY_IN_SECONDS + seconds
                flashesFinishTime(i) = flashesStartTime(i) + dt
                read(lockalFin, *) flashesUpperRho(i), flashesBottomRho(i)
                
                call estimateLockalHeating(flashesEnergy(i), flashesUpperRho(i), flashesBottomRho(i), &
                                           dt, QLockal, upperI, bottomI)
                flashesEmissivity(i) = QLockal
                flashesUpperKnots(i) = upperI
                flashesBottomKnots(i) = bottomI               
            end do
            
            close(lockalFin)
            
            call setFlashesNumber(flashesNumber)
            call setFlashesParameters(flashesEnergy, flashesStartTime, flashesFinishTime, flashesUpperRho, flashesBottomRho, &
                                      flashesEmissivity, flashesUpperKnots, flashesBottomKnots)
        else
            call setFlashesNumber(0)
        end if
    end subroutine
    
    ! estimate SF critical temperatures from model numbers
    ! Input: fileIsSource 
    !        fileIsSource == TRUE --> read model numbers from superfluidFile
    !        fileIsSource == FALSE --> read model numbers from command line
    subroutine estimateSFCriticalTemperature(fileIsSource)
        logical :: fileIsSource
        
        integer(4), parameter :: SFfin = 501
        
        real(8) :: neutron3P2Tc(spatialKnotsNumber), neutron1S0Tc(spatialKnotsNumber), &
                   proton1S0Tc(spatialKnotsNumber)
        
        if (SFKey == 1) then
            if (fileIsSource) then 
                open(SFfin, file = superfluidFile)
                read(SFfin, *) neutronTripletSFModel, neutronSingletSFModel, protonSFModel
                close(SFfin)
            else
                call setCurvesTypesFromCL()
            endif
		
            call SFCURV(KEOS, neutronTripletSFModel, spatialKnotsNumber, rhoDistribution, neutron3P2Tc) ! called from sfgaps.f
            call SFCURV(KEOS, neutronSingletSFModel, spatialKnotsNumber, rhoDistribution, neutron1S0Tc)
            call SFCURV(KEOS, protonSFModel, spatialKnotsNumber, rhoDistribution, proton1S0Tc)
         
        else 
            do i = 1, spatialKnotsNumber
                neutron3P2Tc(i) = 0.0d0
                neutron1S0Tc(i) = 0.0d0
                proton1S0Tc(i) = 0.0d0
            end do
        end if
        
        SFDescription = getMSFDescription()
        call setSFCurves(neutron3P2Tc, neutron1S0Tc, proton1S0Tc)    

        contains
        
        subroutine setCurvesTypesFromCL()
            write(*, *) 'Please, enter the SF type`s numbers for n3P2, n1S0, p1S0 in single line'
            write(*, *) 'n3P2 = 19 ... 26  n1S0 = 1 ... 9  p1S0 = 10 ... 18' 
            write(*, *) 'See gaps.d for more details'
            read(*, *) neutronTripletSFModel, neutronSingletSFModel, protonSFModel
            do while((neutronTripletSFModel < 19) .OR. (neutronTripletSFModel > 26) .OR. &
                  (neutronSingletSFModel < 1) .OR. (neutronSingletSFModel > 9) .OR. &
                  (protonSFModel < 10) .OR. (protonSFModel > 18))
                write(*, *) 'n3P2 = 19 ... 26  n1S0 = 1 ... 9  p1S0 = 10 ... 18'
                read(*, *) neutronTripletSFModel, neutronSingletSFModel, protonSFModel
            end do
            write(*, *) 'Problem parameters is set successfully'
        end subroutine
    end subroutine
    
    ! returns mechanical model description
    function getMSFDescription()
        character(len = 200) :: getMSFDescription
        if (SFKey == 1) then
            write(getMSFDescription, '(a, a, i3, a, i3, a, i3)') & 
                ' SF model`s numbers: ', &
                '\n n3P2 = ', neutronTripletSFModel, &
                '\n n1S0 = ', neutronSingletSFModel, &
                '\n p1S0 = ', protonSFModel
        else
            write(getMSFDescription, '(a)') 'The model does not account for superfluidity'
        end if
    end function
    
    ! prints mechanical model to mechanicalOutput
    subroutine printMechanicalModel()
        call printMm(mM, mechanicalOutput)
    end subroutine
	
    ! write descriptions for file with solutions
    subroutine preapareFiles()
        open(fout, file = outputfile)
        
        write(fout, '(a)') mMDescription
        write(fout, '(a, e12.6)') ' sigma = ', sigma
        write(fout, '(a)') SFDescription
        write(fout, '(a)') '#' 
        write(fout, '(i8)') spatialKnotsNumber
        
        close(fout)
        
        open(fout, file = outputfileForFlux)
        
        write(fout, '(a)') mMDescription
        write(fout, '(a, e12.6)') ' sigma = ', sigma
        write(fout, '(a)') SFDescription
        write(fout, '(a)') '#' 
        
        close(fout)
    end subroutine
    
    ! knits solution on the grid with description
	!subroutine knitSolutions()
    !    integer :: Nout
    !    integer :: j
    !    real(8) :: time, timeInYear, radius, thermal, capacity, conductivity, sources
        
    !    open(fin, file = heatOutputfile)
    !    open(fout, file = finalOutputfile)
        
    !    write(fout, '(a)') mMDescription
    !    write(fout, '(a, e12.6)') ' sigma = ', sigma
    !    write(fout, '(a)') SFDescription
    !    write(fout, '(a)') '#' 
    !    write(fout, '(i8)') spatialKnotsNumber 
        
    !    read(fin, *) Nout
 
    !    do i = 1, Nout
    !        write(fout, '(a, a)') '   time           timeInYear       radius         rho            temperature     capacity', &
    !               '       conductivity       sources'
    !        do j = 1, spatialKnotsNumber
    !            read(fin, *) time, timeInYear, radius, thermal, capacity, conductivity, sources
    !            write(fout, '(e15.6, e15.6, e15.6, e15.6, e18.8, e15.6, e15.6, e15.6)') time, timeInYear, radius, & 
    !                  rhoDistribution(j), thermal, capacity, conductivity, sources
    !        end do
    !    end do
        
    !    close(fout)
    !    close(fin)
        
    !end subroutine 
    
    ! knits heat flux from surface with description 
	!subroutine knitFlux()
    !    character :: columns
    !    integer(4) :: stepNumbers
    !    real(8) :: time, timeInYear, surfaceFlux, Tg, Ts, neutrinoLosses  
        
    !    open(fin, file = fluxFile)
    !    open(fout, file = finalFluxFile)
        
    !    write(fout, '(a)') mMDescription
    !    write(fout, '(a, e12.6)') ' sigma = ', sigma
    !    write(fout, '(a)') SFDescription
    !    write(fout, '(a)') '   time                timeInYear       surfaceFlux    Tg = 10d10 g/cc       Ts    &
    !                                    neutrinoLosses'
        
    !    read(fin, *) columns
    !    read(fin, *) stepNumbers

    !    do i = 1, stepNumbers
    !        read(fin, *) time, timeInYear, surfaceFlux, Tg, Ts, neutrinoLosses
    !        write(fout, '(e22.13, e15.6, e15.6, e15.6, e15.6, e15.6)') time, timeInYear, &
    !                                surfaceFlux, Tg, Ts, neutrinoLosses
    !    end do
        
    !    close(fout)
    !    close(fin)
        
    !end subroutine
	
end module