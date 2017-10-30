volScalarField ycomp = alpha1.mesh().C().component(1);
scalarField Vcells = alpha1.mesh().V();
Vsatellite.value() = 0.0;
Vmain.value() = 0.0;

forAll (alpha1.internalField(), nCells)
{
        if(ycomp[nCells] > 0.6e-4 && ycomp[nCells] < 1.1e-4)
                Vsatellite.value() = Vsatellite.value() + alpha1.internalField()[nCells] * Vcells[nCells];
        else
                Vmain.value() = Vmain.value() + alpha1.internalField()[nCells] * Vcells[nCells];

}

reduce(Vsatellite,sumOp<dimensionedScalar>());
reduce(Vmain,sumOp<dimensionedScalar>());

VolFile << Vsatellite.value() << '\t' << Vmain.value() <<nl;
