
/* This file obtains the cellID of interfacial cells and 
   evaluates the interface area and interface center from
   the sub-cell interface reconstruction provided in the
   isoAdvector solver.
*/

// Obtain the Cell ID of interfacial cells
surfCells = advector.interfacecells();

// Obtain the Area and Center of the reconstructed interface in each cell
forAll (surfCells,i)
{
	interfaceArea[surfCells[i]] = advector.AreaCalc(surfCells[i]);
	interfaceCenter[surfCells[i]] = advector.CentreCalc(surfCells[i]);
}

// Controlling the writeouts to a file
// Check if we need to write the interface details
const dictionary dict(mesh.solverDict(alpha1.name()));
bool InterfaceWrite(dict.lookupOrDefault<bool>("writeInterfaceDetails",false));

if(mesh.time().writeTime() && InterfaceWrite){

    mixture.correct();	// This will update the kappa values if momentum coupling is switched off
    kappa = mixture.K();	// This updates the kappa field, which will be written out at each write time

    // Creating directory Interface Details: Assumption is that openfoam construct accounts for a pre-existing directory
    // with the same name, i.e., it will not create the directory again
    const fileName dirName
    (
        mesh.time().path()/"InterfaceDetails"
    );

    // Identifying the area and center logs by time index
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

gradalpha = fvc::grad(alpha1);	
interfaceArea = vector::zero;	
interfaceCenter = vector::zero;	
surfCells.clear();	     
