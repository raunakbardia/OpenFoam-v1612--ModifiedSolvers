/* This file obtains the cellID of interfacial cells and 
   evaluates the interface area and interface center from
   the sub-cell interface reconstruction provided in the
   isoAdvector solver.
*/

// Creating Dictionary to read if we want to write out interfacial cells
surfCells = advector.interfacecells();

// Controlling the writeouts to a file
// Check if we need to write the interface details
const dictionary dict(mesh.solverDict(alpha1.name()));
const bool InterfaceWrite(dict.lookupOrDefault<bool>("writeInterfaceDetails",false));

forAll (surfCells,i)
{
	vector Area = advector.AreaCalc(surfCells[i]);
	xdot[surfCells[i]] = - Area/mag(Area) * massFlux[surfCells[i]]/rho1.value(); 
	interface_upstream_source[surfCells[i]] = (massFlux[surfCells[i]]/rho1.value() * mag(Area)/cellVol[surfCells[i]]) * mesh.time().deltaTValue();
	if(mesh.time().writeTime() && InterfaceWrite)
	{
		interfaceArea[surfCells[i]] = Area;
		interfaceCenter[surfCells[i]] = advector.CentreCalc(surfCells[i]);
	}
}

if(mesh.time().writeTime() && InterfaceWrite)
{
    const fileName dirName
    (
        mesh.time().path()/"InterfaceDetails"
    );
    #include "writeInterface.H"
}
