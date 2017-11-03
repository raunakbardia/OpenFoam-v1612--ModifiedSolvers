// Update: November 1, 2017
/*  This file updates the mdot values based on the
    interface area calculated in the included header
    file using the library functionalities of the
    isoAdvector solver.
*/

#include "InterfaceDetails.H"

forAll (surfCells,i)
	AlphaSourceTerm[surfCells[i]] = PhaseChangeFlux[surfCells[i]] * mag(interfaceArea[surfCells[i]])/cellvol[surfCells[i]];

interfaceArea = vector::zero;	
interfaceCenter = vector::zero;	
surfCells.clear();	     
