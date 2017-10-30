// Update: October 27, 2017
/*  This file updates the mdot values based on the
    interface area calculated in the included header
    file using the library functionalities of the
    isoAdvector solver.
*/

#include "InterfaceArea.H"

forAll (surfCells,i)
{
	mdot[surfCells[i]] = Minterface.value() * mag(interfaceArea[i])/cellvol[surfCells[i]];
	mdot_source[surfCells[i]] = mdot[surfCells[i]]/rho1.value();
	Info << surfCells[i] <<'\n';
}
