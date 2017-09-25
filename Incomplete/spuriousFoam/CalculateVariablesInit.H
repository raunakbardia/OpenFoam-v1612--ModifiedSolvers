volVectorField centres = alpha1.mesh().C();
volScalarField xcomp = centres.component(0);
volScalarField ycomp = centres.component(1);
volScalarField zcomp = centres.component(2);

//Calculating COM
COM_x_init = fvc::domainIntegrate(alpha1*xcomp)/fvc::domainIntegrate(alpha1);
COM_y_init = fvc::domainIntegrate(alpha1*ycomp)/fvc::domainIntegrate(alpha1);
COM_z_init = fvc::domainIntegrate(alpha1*zcomp)/fvc::domainIntegrate(alpha1);

//Calculating Kinetic Energy over the whole domain
KE = 0.5*fvc::domainIntegrate(rho*(U & U));

COMFile

	<< COM_x_init.value() << '\t' << COM_y_init.value() << '\t' << COM_z_init.value() << nl;

KEFile

	<< KE.value() << nl;

RunTimeFile

        << runTime.timeName() << nl;

