/* This file obtains the cellID of interfacial cells and 
   evaluates the interface area and interface center from
   the sub-cell interface reconstruction provided in the
   isoAdvector solver.
*/

// Creating Dictionary to read if we want to write out interfacial cells
surfCells = advector.interfacecells();

forAll (surfCells,i)
{
	interfaceArea[surfCells[i]] = advector.AreaCalc(surfCells[i]);
	interfaceCenter[surfCells[i]] = advector.CentreCalc(surfCells[i]);
}

