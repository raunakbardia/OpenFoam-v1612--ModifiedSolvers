
const bool InterfaceWrite(dict.lookupOrDefault<bool>("writeInterfaceDetails",false));

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
