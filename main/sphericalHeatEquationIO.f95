module sphericalHeatEquationIO

use sphericalHeatEquation

implicit none
    
    integer, private, parameter :: fin = 10, finSection = 25
    integer, private, parameter :: fout = 20
    
    ! initial conditions
    integer(4), private :: spatialKnotsNumber
    
    real(8), private :: sigma, t0, initialTau

    real(8), private, allocatable :: spatialKnots(:)
    real(8), private, allocatable :: rhoDistribution(:)
    real(8), private, allocatable :: initialCondition(:)
    
    ! output conditions
    integer(4), private :: outputSectionsNumber ! number of output sections
    
    real(8), allocatable :: outputSectionsTime(:) ! time of output sections
    
contains 
    
    ! prepares difference scheme from file
    subroutine constructSHEFromFile(hE, inputfile)
        type(heatEquation) :: hE
        character(*) :: inputfile ! This file contains weight -- sigma, T_{init} -- initial temperature of star 
                                  ! and model grid (number of knots and their coordinates and rhoDistribution)
        
        call setSHEParametersFromFile()
        
        contains
        
        subroutine setSHEParametersFromFile()            
            integer(4) :: spatialKnotsNumber   
            character(len = 50) :: outputParameters = 'inputData\\outputParameters.txt'
            
            integer(4) :: i
            real(8) :: prevTimePoint, timePoint, years, days, seconds
            
            open(fin, file = inputfile)
        
            read(fin, *) sigma, initialTau
            read(fin, *) spatialKnotsNumber, t0
            
            ! Let's read initial condition
            allocate(spatialKnots(spatialKnotsNumber), initialCondition(spatialKnotsNumber), & 
                     rhoDistribution(spatialKnotsNumber))
            
            do i = 1, spatialKnotsNumber
                read(fin, *) spatialKnots(i), rhoDistribution(i), initialCondition(i)
            end do
            
            close(fin)
            
            ! Let's read output parameters 
            open(finSection, file = outputParameters)
            
            read(finSection, *) outputSectionsNumber
            allocate(outputSectionsTime(outputSectionsNumber))
            
            do i = 1, outputSectionsNumber
                read(finSection, *) years, days, seconds
                outputSectionsTime(i) = seconds + days * DAY_IN_SECONDS + years * YEAR_IN_SECONDS 
            end do
            
            close(finSection)
            
            call setSHEInitialConditions(hE, spatialKnotsNumber, sigma, t0, initialTau, spatialKnots, & 
                                         rhoDistribution, initialCondition) 
            ! call setSHEStepConditions(hE, stepNumbers, timePeriodsNumber, timePeriods, numbersOfPeriods)
            call setSHEOutputConditions(hE, outputSectionsNumber, outputSectionsTime)
        
        end subroutine
        
    end subroutine
	
    ! prints solution 
    ! Input: hE - our model
    !        outputFile - path of outputFile
    !        isAppend == TRUE --> appends result to file
    !        isAppend == FALSE --> rewrite file
	subroutine printHESolution(hE, outputFile, isAppend)
        type(heatEquation) :: hE
        character(*) :: outputFile
        logical :: isAppend  
        
        real(8) :: outputSections(spatialKnotsNumber, outputSectionsNumber)
        
        integer(4) :: i, j
        
        if (isAppend) then
            open(fout, file = outputFile, status = "old", position = "append")
        else
            open(fout, file = outputFile)
        end if
        
        outputSections = getSHEOutputSections(hE)
        
        write(fout, '(i8)') outputSectionsNumber
        do i = 1, outputSectionsNumber
            do j = 1, spatialKnotsNumber
                write(fout, '(e19.7, e19.8, e19.8)') outputSectionsTime(i), spatialKnots(j), outputSections(j, i)
            end do
        end do
        close(fout)
    end subroutine
    
end module