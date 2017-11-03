#include "createRDeltaT.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();


// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();


// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);


// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture)
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

// MULES flux from previous time-step
surfaceScalarField alphaPhi
(
    IOobject
    (
        "alphaPhi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    phi*fvc::interpolate(alpha1)
);

// MULES compressed flux is registered in case scalarTransport FO needs it.
surfaceScalarField alphaPhiUn
(
    IOobject
    (
        "alphaPhiUn",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", phi.dimensions(), 0.0)
);

// MULES Correction
tmp<surfaceScalarField> talphaPhiCorr0;

#include "createMRF.H"


// Output details for Interface
// Additions by Raunak Bardia - Nov 1, 2017

const scalarField cellvol = mesh.V();	// This limits the implementation to static mesh only 
const vectorField cellcenter = mesh.C();	// This limits the implementation to static mesh only 
const label meshsize = cellcenter.size();

DynamicList<label> surfCells;	// Stores Cell Id of all interfacial cells

// To get Interface Area Vectors for the entire mesh
// It is zero for non-interfacial cells - easier to read in MATLAB
List<vector> interfaceArea(meshsize);

// To get Interface Center Position for the entire mesh
// It is (-100,-100,-100) for non-interfacial cells - easier to read in MATLAB
List<point> interfaceCenter(meshsize);

/* Calculation of mdot value for phase change: Corresponds with Eq2 and Eq 20
	of Sato & Niceno (2013), JCP */

volScalarField PhaseChangeFlux
(
    IOobject
    (
        "mdot",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(1,-2,-1,0,0,0,0), -0.1)
);
volScalarField AlphaSourceTerm
(
    IOobject
    (
        "mdotsource",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0,0,-1,0,0,0,0), 0)
);

// Instantiating isoAdvection object for advancing alpha1
isoAdvection advector(alpha1,AlphaSourceTerm,phi,U);

// Creating Dictionary to read if we want to write out interfacial cells
const dictionary dict(mesh.solverDict(alpha1.name()));
bool InterfaceWrite(dict.lookupOrDefault<bool>("writeInterfaceDetails",false));

// * * * * * * September 25, 2017 * * * * * * //
/*
// * * * * * * September 13, 2017 * * * * * * //

// Temperature
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

dimensionedScalar Cvl(PhaseChangeProperties.lookup("Cv_l"));
dimensionedScalar Cpv(PhaseChangeProperties.lookup("Cp_v"));
dimensionedScalar Cvv(PhaseChangeProperties.lookup("Cv_v"));
dimensionedScalar kl(PhaseChangeProperties.lookup("k_l"));
dimensionedScalar kv(PhaseChangeProperties.lookup("k_v"));
dimensionedScalar hfg(PhaseChangeProperties.lookup("h_fg"));
dimensionedScalar Tsat(PhaseChangeProperties.lookup("T_sat"));

// Calculate mdot at the interface
volScalarField mdotGamma
(
    IOobject
    (
        "modtGamma",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    mesh,
    dimensionedScalar("zero",dimensionSet(1,-2,-1,0,0,0,0),0.0)	//mdot calculated at area
);

// Calculate interface area for all interface cells
volScalarField intA
(
    IOobject
    (
        "intA",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    mesh,
    dimensionedScalar("zero",dimensionSet(1,-2,-1,0,0,0,0),0.0)	//interface area will be calculated
);

// Volumetric form of mdot
volScalarField mdotV
(
    IOobject
    (
        "modtV",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    mesh,
    dimensionedScalar("zero",dimensionSet(1,-3,-1,0,0,0,0),0.0)	//mdot calculated at area
);

// Need to store Thermal Conductivity (k)
volScalarField k
(
    IOobject
    (
        "k",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*kl + alpha2*kv
);

// Need to store rhoCv for ddt(rhoCv,T)
volScalarField rhoCv
(
    IOobject
    (
        "rhoCv",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1*Cvl + alpha2*rho2*Cvv
);

// Need to store rhoCp for energy convection
volScalarField rhoCp
(
    IOobject
    (
        "rhoCp",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1*Cpl + alpha2*rho2*Cpv
);

// Energy Flux at Cell Boundaries
surfaceScalarField HPhi
(
    IOobject
    (
        "HPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rhoCp)*phi
);

dimensionedVector COM_y = fvc::domainIntegrate(alpha1*mesh.C()) / fvc::domainIntegrate(alpha1);
OFstream COMFile
(
    runTime.path()/("COM_Log")
);*/
