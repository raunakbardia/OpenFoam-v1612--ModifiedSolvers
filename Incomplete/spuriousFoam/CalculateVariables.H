volVectorField centres = alpha1.mesh().C();
volScalarField xcomp = centres.component(0);
volScalarField ycomp = centres.component(1);
volScalarField zcomp = centres.component(2);

//Calculating COM
COM_x = fvc::domainIntegrate(alpha1*xcomp)/fvc::domainIntegrate(alpha1);
COM_y = fvc::domainIntegrate(alpha1*ycomp)/fvc::domainIntegrate(alpha1);
COM_z = fvc::domainIntegrate(alpha1*zcomp)/fvc::domainIntegrate(alpha1);

//Calculating Kinetic Energy over the whole domain
KE = 0.5*fvc::domainIntegrate(rho*(U & U));

ShiftCOM = sqrt(pow(COM_x-COM_x_init,2.0)+pow(COM_y-COM_y_init,2.0)+pow(COM_z-COM_z_init,2.0));

COMFile

	<< COM_x.value() << '\t' << COM_y.value() << '\t' << COM_z.value() << '\t' << ShiftCOM.value() << nl;

KEFile

	<< KE.value() << nl;

RunTimeFile

        << runTime.timeName() << nl;

