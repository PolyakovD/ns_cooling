objects = main.o sphericalHeatEquationIO.o sphericalHeatEquation.o tridiagonal.o calculations.o modelIO.o model.o testModel.o mechanicalModelIO.o mechanicalModel.o bskfit.o tovbskfit.o cvfull.o compos21.o eos14.o eosmag14.o condBSk21.o neutrinos.o sfgaps.o tbtsfb.o coreConductivity.o
          
heatEquationSolution: main.exe
	main.exe
main.exe: $(objects)
	gfortran $(objects) -o main
main.o: main.f95 sphericalHeatEquation.o tridiagonal.o model.o modelIO.o testModel.o
	gfortran -c main.f95
sphericalHeatEquationIO.o: sphericalHeatEquationIO.f95 sphericalHeatEquation.o
	gfortran -c sphericalHeatEquationIO.f95
sphericalHeatEquation.o: sphericalHeatEquation.f95 tridiagonal.o calculations.o
	gfortran -c sphericalHeatEquation.f95
tridiagonal.o: tridiagonal.f95
	gfortran -c tridiagonal.f95
calculations.o: calculations.f95
	gfortran -c calculations.f95

modelIO.o: modelIO.f95 model.o	mechanicalModelIO.o sfgaps.o
	gfortran -c -fbackslash modelIO.f95
model.o: model.f95 mechanicalModelIO.o mechanicalModel.o sphericalHeatEquationIO.o sphericalHeatEquation.o tovbskfit.o cvfull.o condBSk21.o neutrinos.o sfgaps.o tbtsfb.o coreConductivity.o
	gfortran -c -fbackslash model.f95   
mechanicalModelIO.o: mechanicalModelIO.f95 tovbskfit.o mechanicalModel.o
	gfortran -c -fbackslash mechanicalModelIO.f95
mechanicalModel.o: mechanicalModel.f95
	gfortran -c mechanicalModel.f95
tovbskfit.o: tovbskfit.f95 bskfit.o
	gfortran -c tovbskfit.f95
sfgaps.o: sfgaps.f bskfit.o gaps.d
	gfortran -c sfgaps.f
bskfit.o: bskfit.f 
	gfortran -c bskfit.f
tbtsfb.o: tbtsfb.f
	gfortran -c tbtsfb.f
coreConductivity.o: coreConductivity.f95
	gfortran -c coreConductivity.f95
    
testModel.o: testModel.f95 sphericalHeatEquation.o sphericalHeatEquationIO.o
	gfortran -c testModel.f95
	
cvfull.o: cvfull.f compos21.o eos14.o eosmag14.o bskfit.o
	gfortran -c cvfull.f
compos21.o: compos21.f
	gfortran -c compos21.f
eos14.o: eos14.f 
	gfortran -c eos14.f
eosmag14.o: eosmag14.f
	gfortran -c eosmag14.f
	
condBSk21.o: condBSk21.f
	gfortran -c condBSk21.f
    
neutrinos.o : neutrinos.f95
	gfortran -c neutrinos.f95
clear: $(objects)
	rm *.o *.mod *.exe
