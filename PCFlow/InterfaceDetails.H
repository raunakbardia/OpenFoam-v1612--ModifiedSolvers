/* This file obtains the cellID of interfacial cells and 
   evaluates the interface area and interface center from
   the sub-cell interface reconstruction provided in the
   isoAdvector solver.
*/

surfCells = advector.interfacecells();

forAll (surfCells,i)
{
	interfaceArea[surfCells[i]] = advector.AreaCalc(surfCells[i]);
	interfaceCenter[surfCells[i]] = advector.CentreCalc(surfCells[i]);
}

if(mesh.time().writeTime() && InterfaceWrite){
    const fileName dirName
    (
        mesh.time().path()/"InterfaceDetails"
    );
    const string fName1
    (
        "InterfaceAreaLog_" + Foam::name(mesh.time().timeIndex())
    );
    const string fName2
    (
        "InterfaceCenterLog_" + Foam::name(mesh.time().timeIndex())
    );
    mkDir(dirName);
    OFstream InterfaceAreaFile(dirName/fName1);
    OFstream InterfaceCenterFile(dirName/fName2);
    forAll (alpha1.internalField(),nCells)
    {
	InterfaceCenterFile << interfaceCenter[nCells].component(0) << "," << interfaceCenter[nCells].component(1) << "," << interfaceCenter[nCells].component(2) << '\n';
	InterfaceAreaFile << interfaceArea[nCells].component(0) << "," << interfaceArea[nCells].component(1) << "," << interfaceArea[nCells].component(2) << '\n';
    }
}
