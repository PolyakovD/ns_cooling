! Needs to be linked with: model.f95 sphericalHeatEquation.f95 tridiagonal.f95
program program1
    !use testModel and model separately to avoid naming collisions. 
    !use testModel
    use model
    use modelIO
	
	implicit none
   
    !call conditionGeneratorFromFiles()
    call conditionGeneratorFromCL()
    call makeSolution()
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Let's solve the second problem from sphericalHeatEquation2.pdf
    !character(len = 50) :: configFile = 'inputData\\data.txt' 
    
    !call conditionGenerator01(configFile)
    !call makeSolution()
	
end program program1
