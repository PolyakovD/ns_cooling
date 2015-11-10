module sphericalHeatEquationIO

use sphericalHeatEquation

implicit none
    
    integer, private, parameter :: fin = 10, finSection = 25
    integer, private, parameter :: fout = 20
    
    ! initial conditions
    integer(4), private :: spatialKnotsNumber
    
    real(8), private :: sigma, t0

    real(8), private, allocatable, dimension(:) :: spatialKnots
    real(8), private, allocatable, dimension(:) :: initialCondition
    
    ! step conditions
    integer(4), private :: stepNumbers 
    integer(4), private :: timePeriodsNumber
    
    real(8), private, allocatable, dimension(:) :: timePeriods
    integer(4), private, allocatable, dimension(:) :: numbersOfPeriods
    
    ! output conditions
    integer(4), private :: sectionsNumber
    integer(4), private, allocatable, dimension(:) :: sections
    
    private :: setSHEParametersFromFile
    
contains 

    subroutine constructSHEFromFile(hE, inputfile)
        type(heatEquation) :: hE
        character(*) :: inputfile ! This file contains weight -- sigma, T_{init} -- initial temperature of star 
                                  ! and model grid (number of knots and their coordinates)
        
        call setSHEParametersFromFile(hE, inputfile)
        
    end subroutine
            
    subroutine setSHEParametersFromFile(hE, inputfile)
        type(heatEquation) :: hE
        character(*) :: inputfile ! This file contains weight -- sigma, t_{o} -- initial time 
                                  ! and model grid (number of knots and their coordinates and temperature)
        
        integer(4) :: spatialKnotsNumber   
        character(len = 50) :: outputParameters = 'inputData\\outputParameters.txt'
        
        integer(4) :: i
        real(8) :: prevTimePoint, timePoint, years, days, seconds
        
        open(fin, file = inputfile)
	
        read(fin, *) sigma
        read(fin, *) spatialKnotsNumber, t0
        
        ! Let's read initial condition
        allocate(spatialKnots(spatialKnotsNumber))
        allocate(initialCondition(spatialKnotsNumber))
        
        do i = 1, spatialKnotsNumber
            read(fin, *) spatialKnots(i), initialCondition(i)
        end do
        
        close(fin)
        
        ! Let's read grid parameters 
        open(finSection, file = outputParameters)
        
        read(finSection, *) timePeriodsNumber
        allocate(timePeriods(timePeriodsNumber))
        allocate(numbersOfPeriods(timePeriodsNumber))
        
        sectionsNumber = timePeriodsNumber + 1
        allocate(sections(sectionsNumber))
        
        stepNumbers = 1
        sections(1) = stepNumbers
		prevTimePoint = t0
        do i = 1, timePeriodsNumber
            read(finSection, *) years, days, seconds, numbersOfPeriods(i)
    
            stepNumbers = stepNumbers + numbersOfPeriods(i)
            sections(i + 1) = stepNumbers
			timePoint = seconds + days * DAY_IN_SECONDS + years * YEAR_IN_SECONDS 
            timePeriods(i) = (timePoint - prevTimePoint) / numbersOfPeriods(i)
            prevTimePoint = timePoint           
        end do
        
        close(finSection)
        
        call setSHEInitialConditions(hE, spatialKnotsNumber, sigma, t0, spatialKnots, initialCondition) 
        call setSHEStepConditions(hE, stepNumbers, timePeriodsNumber, timePeriods, numbersOfPeriods)
        call setSHEOutputConditions(hE, sectionsNumber, sections)
    
    end subroutine
	
	subroutine printHESolution(hE, outputFile, isAppend)
        type(heatEquation) :: hE
        character(*) :: outputFile
        logical :: isAppend  
        
        real(8) :: timeOfSections(sectionsNumber), outputSections(spatialKnotsNumber, sectionsNumber)
        
        integer(4) :: i, j
        
        if (isAppend) then
            open(fout, file = outputFile, status = "old", position = "append")
        else
            open(fout, file = outputFile)
        end if
        
        timeOfSections = getSHETimeOfSections(hE)
        outputSections = getSHEOutputSections(hE)
        
        write(fout, '(i8)') sectionsNumber
        do i = 1, sectionsNumber
            do j = 1, spatialKnotsNumber
                write(fout, '(e19.7, e19.8, e19.8)') timeOfSections(i), spatialKnots(j), outputSections(j, i)
            end do
        end do
        close(fout)
    end subroutine
    
end module