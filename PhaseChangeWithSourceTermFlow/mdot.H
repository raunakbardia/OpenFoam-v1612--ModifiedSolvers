// Update: November 1, 2017
/*  This file updates the mdot values based on the
    interface area calculated in the included header
    file using the library functionalities of the
    isoAdvector solver.
*/

#include "InterfaceDetails.H"


forAll (surfCells,i)
    mdot[surfCells[i]] = PhaseChangeFlux[surfCells[i]] * mag(interfaceArea[surfCells[i]])/cellvol[surfCells[i]];
//    AlphaSourceTerm[surfCells[i]] = PhaseChangeFlux[surfCells[i]] * mag(interfaceArea[surfCells[i]])/cellvol[surfCells[i]];

//AlphaSourceTerm = mdot/rho1;

const dictionary dict(mesh.solverDict(alpha1.name()));
const scalar extrapolationIterations(dict.lookupOrDefault<scalar>("extrapolationIterations", 1000));
const scalar extrapolationTol(dict.lookupOrDefault<scalar>("extrapolationTol", 1e-8));

for(subCycle<volScalarField> mdotSC(mdot,extrapolationIterations);!(++mdotSC).end();)
{
	fvScalarMatrix mdotEqn
	(
    		fvm::ddt(mdot)
    		+ fvm::laplacian(liquidorvapor,mdot)
	);
	mdotEqn.solve();
}
NormalVector = fvc::grad(alpha1)/(mag(fvc::grad(alpha1))+deltaN); 

interfaceArea = vector::zero;	
interfaceCenter = vector::zero;	
surfCells.clear();	     
