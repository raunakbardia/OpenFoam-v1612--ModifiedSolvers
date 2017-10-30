
/* This file obtains the cellID of interfacial cells and 
   evaluates the interface area and interface center from
   the sub-cell interface reconstruction provided in the
   isoAdvector solver.
*/

surfCells = advector.interfacecells();

forAll (surfCells,i)
{
	interfaceArea[i] = advector.AreaCalc(surfCells[i]);
	interfaceCenter[i] = advector.CentreCalc(surfCells[i]);
	Info << surfCells[i] << '\n';
}
