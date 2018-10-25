
Diff = DC * (alpha1 + alpha2) ; 
// volScalarField phiH = ("phiH", Diff*alpha1 ) ; 
        {
            fvScalarMatrix CEqn
            (
	     fvm::ddt(C)
		+ fvm::div(phi, C)
		- fvm::laplacian(DC, C)
		- fvm::div(phiHsurf,C)
		//              - fvm::laplacian(DC, alpha1)
             ==
                fvOptions(C)
            );

            CEqn.relax();
            fvOptions.constrain(CEqn);
            CEqn.solve();
            fvOptions.correct(C);
        }

